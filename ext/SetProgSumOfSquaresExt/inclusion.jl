
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
           !(subset.degree isa Int) ||
           !(supset.degree isa Int) ||
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

# These certificates use the defining homogeneous polynomials directly, so
# polynomial set variables are represented in the primal space.
const _PolynomialLpSource =
    Union{SetProg.SetVariableRef{<:Any,<:SetProg.PolySet},Sets.PolySet,Sets.ConvexPolySet}
function SetProg.set_space(
    space::SetProg.Space,
    ::InclusionConstraint{<:_PolynomialLpSource,<:Sets.PolarLpSum},
)
    return SetProg.set_space(space, SetProg.PrimalSpace)
end

function _integer_polynomial_power(poly, ::Val{true}, exponent)
    return poly^Int(exponent)
end
function _homogeneous_polynomial_power(poly, degree::Int, p)
    exponent = p / degree
    return _integer_polynomial_power(poly, Val(isinteger(exponent)), exponent)
end
function _powered_gauge_polynomial(set::Union{Sets.PolySet,Sets.ConvexPolySet}, x, p)
    poly = MP.polynomial(set.p)
    return _homogeneous_polynomial_power(poly(MP.variables(poly) => x), set.degree, p)
end
function _powered_gauge_polynomial(set::Sets.Ellipsoid, x, p)
    return _homogeneous_polynomial_power(dot(x, set.Q * x), 2, p)
end
function _powered_gauge_polynomial(set::Sets.LinearPreImage, x, p)
    return _powered_gauge_polynomial(set.set, set.A * x, p)
end
function _powered_gauge_polynomial(set::Sets.PolarLpSum, x, p)
    # A nested sum must use the same exponent to remain polynomial without
    # taking roots of sums. Other cases need a different certificate.
    return _nested_powered_gauge_polynomial(Val(set.p == p), set, x, p)
end
function _nested_powered_gauge_polynomial(::Val{true}, set, x, p)
    return _powered_gauge_polynomial(set.first, x, p) +
           _powered_gauge_polynomial(set.second, x, p)
end

function JuMP.build_constraint(
    _error::Function,
    subset::Union{Sets.PolySet,Sets.ConvexPolySet},
    superset::PowerSet{<:Sets.PolarLpSum};
    margin::Real = 0,
)
    margin >= 0 || throw(ArgumentError("margin must be nonnegative"))
    if SetProg.need_variablify(superset.set)
        return InclusionConstraint(subset, superset.set, (; margin))
    end
    x = Sets.space_variables(subset)
    p = superset.set.p
    residual =
        _powered_gauge_polynomial(subset, x, p) -
        _powered_gauge_polynomial(superset.set, x, p)
    if !iszero(margin)
        residual -= margin * _homogeneous_polynomial_power(sum(x .^ 2), 2, p)
    end
    return JuMP.build_constraint(_error, residual, SOSCone())
end
