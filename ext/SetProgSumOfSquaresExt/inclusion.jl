
# S-procedure: Q ⊆ P <=> q(x) ≤ 1 => p(x) ≤ 1 <=> p(x) ≤ q(x) <= q - p is SOS
function JuMP.build_constraint(
    _error::Function,
    subset::Union{Sets.PolySet,Sets.ConvexPolySet},
    sup_powerset::PowerSet{<:Union{Sets.PolySet,Sets.ConvexPolySet}};
    S_procedure_scaling = nothing,
    margin = 0,
)
    supset = sup_powerset.set
    if subset isa Sets.PolySet &&
       supset isa Sets.PolySet &&
       (
           MP.mindegree(MP.polynomial(subset.p)) != subset.degree ||
           MP.mindegree(MP.polynomial(supset.p)) != supset.degree
       )
        return InclusionConstraint(subset, supset, (; S_procedure_scaling, margin))
    end
    @assert iszero(margin)
    @assert S_procedure_scaling === nothing || isone(S_procedure_scaling)
    q = subset.p
    p = sup_powerset.set.p
    JuMP.build_constraint(_error, q - p, SOSCone())
end

# S-procedure: Q ⊆ P <=> q - p is SOS
function s_procedure(model, subset, supset; S_procedure_scaling = nothing)
    q = subset.p
    p = supset.p
    if S_procedure_scaling === nothing
        S_procedure_scaling = @variable(model)
        # We want to avoid creating a non-convex problem. If one of `p` and `q`
        # is a polynomial with constant coefficients, we multiply the variable
        # by this one
        if MP.coefficient_type(q) <: Number
            s = S_procedure_scaling * q - p
        else
            s = q - S_procedure_scaling * p
        end
    else
        s = q - S_procedure_scaling * p
    end
    return JuMP.build_constraint(error, s, SOSCone())
end


#     [p(x) ≤ 1] ⊆ [⟨a, x⟩ ≤ β]
# <= λ(1 - p(x)) ≤ β - ⟨a, x⟩ (necessary if p is SOS-convex)
#              0 ≤ λ(p(x) - 1) - ⟨a, x⟩ + β
# Set `x = 0`: 0 ≤ λ(p(0) - 1) + β
# hence if p(0) ≤ 1 (i.e. 0 ∈ S), then λ ≤ β / (1 - p(0))
# Homogeneous case: λ ≤ β
# Use build_constraint when SumOfSquares#66 if λ = β (e.g. homogeneous)
function JuMP.add_constraint(
    model::JuMP.Model,
    constraint::InclusionConstraint{
        <:Union{Sets.ConvexPolySet,Sets.ConvexPolynomialSet},
        <:HalfSpace,
    },
    name::String = "",
)
    p = Sets.gauge1(constraint.subset)
    h = constraint.supset
    x = Sets.space_variables(constraint.subset)
    hs = dot(h.a, x) - h.β
    if MP.coefficient_type(p) <: Number
        λ = @variable(model, lower_bound=0.0, base_name = "λ")
        cref = @constraint(model, λ * (p - 1) - hs in SOSCone())
    elseif Polyhedra.coefficient_type(h) <: Number
        λ = @variable(model, lower_bound=0.0, base_name = "λ")
        cref = @constraint(model, (p - 1) - λ * hs in SOSCone())
    else
        β = _unstable_constantify(h.β)
        if β isa Number
            # TODO what is a good value for β ???
            λ = β/2
        else
            λ = @variable(model, lower_bound=0.0, base_name = "λ")
        end
        cref = @constraint(model, λ * (p - 1) - hs in SOSCone())
    end
    return cref
end
