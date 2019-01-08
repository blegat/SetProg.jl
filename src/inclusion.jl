### InclusionConstraint ###
struct InclusionConstraint{SubSetType, SupSetType, KWT} <: SetConstraint
    subset::SubSetType
    supset::SupSetType
    kws::KWT
end
function need_variablify(c::InclusionConstraint)
    return need_variablify(c.subset) || need_variablify(c.supset)
end
function variablify(c::InclusionConstraint)
    return JuMP.build_constraint(error, variablify(c.subset),
                                 PowerSet(variablify(c.supset)); c.kws...)
end
function clear_spaces(c::InclusionConstraint)
    clear_spaces(c.subset)
    clear_spaces(c.supset)
end
function Sets.perspective_variable(c::InclusionConstraint)
    return synchronize_perspective(Sets.perspective_variable(c.subset),
                                   Sets.perspective_variable(c.supset))
end
function create_spaces(c::InclusionConstraint, spaces::Spaces)
    sub = create_spaces(c.subset, spaces)
    sup = create_spaces(c.supset, spaces)
    merge_spaces(spaces, sub, sup)
end

JuMP.function_string(print_mode, c::InclusionConstraint) = string(c.subset)
function JuMP.in_set_string(print_mode, c::InclusionConstraint)
    string(print_mode == JuMP.IJuliaMode ? "\\subseteq" : "⊆", " ", c.supset)
end
struct PowerSet{S}
    set::S
end

# Fallback, might be because `subset` or `sup_powerset` is a `VariableRef` or
# a `Polyhedron` (which is handled by `JuMP.add_constraint`).
function JuMP.build_constraint(_error::Function, subset, sup_powerset::PowerSet;
                               kws...)
    InclusionConstraint(subset, sup_powerset.set, kws)
end

### InclusionConstraint for sets ###

## Set in Set ##
function set_space(space::Space, ::InclusionConstraint{<:LinearImage,
                                                       <:LinearImage})
    return set_space(space, DualSpace)
end
# We can always transform  an ellipsoid to primal or dual space so we can handle
# any space
function set_space(space::Space,
                   ::InclusionConstraint{<:VariableRef,
                                         <:Sets.AbstractEllipsoid{T}}) where T<:Number
    return space
end

# S-procedure: Q ⊆ P <=> xQx ≤ 1 => xPx ≤ 1 <=> xPx ≤ xQx <=> Q - P is PSD
function JuMP.build_constraint(_error::Function,
                               subset::Sets.EllipsoidAtOrigin,
                               sup_powerset::PowerSet{<:Sets.EllipsoidAtOrigin})
    Q = subset.Q
    P = sup_powerset.set.Q
    JuMP.build_constraint(_error, Symmetric(Q - P), PSDCone())
end

# S-procedure: Q ⊆ P <=> q(x) ≤ 1 => p(x) ≤ 1 <=> p(x) ≤ q(x) <= q - p is SOS
function JuMP.build_constraint(_error::Function,
                               subset::Sets.ConvexPolynomialSublevelSetAtOrigin,
                               sup_powerset::PowerSet{<:Sets.ConvexPolynomialSublevelSetAtOrigin})
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
        if MultivariatePolynomials.coefficienttype(q) <: Number
            s = S_procedure_scaling * q - p
        else
            s = q - S_procedure_scaling * p
        end
    else
        s = q - S_procedure_scaling * p
    end
    return JuMP.build_constraint(error, s, SOSCone())
end

# S-procedure: Q ⊆ P <=> q - p is SOS
function JuMP.add_constraint(model::JuMP.Model,
                             constraint::InclusionConstraint{<:Union{Sets.PerspectiveEllipsoid,
                                                                     Sets.PerspectivePolynomialSet},
                                                             <:Union{Sets.PerspectiveEllipsoid,
                                                                     Sets.PerspectivePolynomialSet}},
                             name::String = "")
    return JuMP.add_constraint(model,
                               s_procedure(model, constraint.subset,
                                           constraint.supset;
                                           constraint.kws...), name)
end

# S ⊆ T <=> T* ⊇ S*
function JuMP.build_constraint(_error::Function,
                               subset::Sets.PerspectiveDualOf{<:Union{Sets.PerspectiveEllipsoid,
                                                                      Sets.PerspectivePolynomialSet}},
                               sup_powerset::PowerSet{<:Sets.PerspectiveDualOf{<:Union{Sets.PerspectiveEllipsoid,
                                                                                                Sets.PerspectivePolynomialSet}}};
                               kws...)
    S = subset
    T = sup_powerset.set
    JuMP.build_constraint(_error, Sets.perspective_dual(T), PowerSet(Sets.perspective_dual(S)); kws...)
end

# S ⊆ T <=> polar(T) ⊆ polar(S)
function JuMP.build_constraint(_error::Function,
                               subset::Sets.Polar,
                               sup_powerset::PowerSet{<:Sets.Polar})
    S = subset
    T = sup_powerset.set
    JuMP.build_constraint(_error, Sets.polar(T), PowerSet(Sets.polar(S)))
end

# See [LTJ18]
function JuMP.build_constraint(_error::Function,
                               subset::LinearImage{<:Sets.AbstractSet},
                               sup_powerset::PowerSet{<:LinearImage{<:Sets.AbstractSet}};
                               kws...)
    JuMP.build_constraint(_error, apply_map(subset),
                          PowerSet(apply_map(sup_powerset.set)); kws...)
end

## Set in Polyhedron ##
function set_space(space::Space, ::InclusionConstraint{<:VariableRef,
                                                       <:Polyhedra.Rep})
    return set_space(space, DualSpace)
end

# Ellipsoid #
function JuMP.add_constraint(model::JuMP.Model,
                             constraint::InclusionConstraint{<:Sets.AbstractSet{JuMP.VariableRef},
                                                             <:Polyhedra.Rep},
                             name::String = "")
    ◯ = constraint.subset
    □ = constraint.supset
    for hp in hyperplanes(□)
        @constraint(model, ◯ ⊆ hp)
    end
    for hs in halfspaces(□)
        @constraint(model, ◯ ⊆ hs)
    end
end

#  S∘ ⊆ [⟨a,   x⟩ ≤ β]
#  S∘ ⊆ [⟨a/β, x⟩ ≤ 1]
# a/β ∈ S
function JuMP.build_constraint(_error::Function, subset::Sets.Polar,
                               sup_powerset::PowerSet{<:Polyhedra.HyperPlane})
    @assert iszero(sup_powerset.set.β) # Otherwise it is not symmetric around the origin
    JuMP.build_constraint(_error, Line(sup_powerset.set.a),
                          Sets.polar(subset))
end
function JuMP.build_constraint(_error::Function, subset::Sets.Polar,
                               sup_powerset::PowerSet{<:Polyhedra.HalfSpace})
    JuMP.build_constraint(_error, ScaledPoint(sup_powerset.set.a, sup_powerset.set.β), Sets.polar(subset))
end

# τ^{-1}(τ(S)*) ⊆ [⟨a, x⟩ ≤ β]
#         τ(S)* ⊆ [⟨(β, -a), (z, x)⟩ ≥ 0]
#       (β, -a) ∈ τ(S)
function JuMP.build_constraint(_error::Function, subset::Sets.PerspectiveDual,
                               sup_powerset::PowerSet{<:Polyhedra.HyperPlane})
    JuMP.build_constraint(_error, SymScaledPoint(-sup_powerset.set.a, sup_powerset.set.β), Sets.polar(subset))
end
function JuMP.build_constraint(_error::Function, subset::Sets.PerspectiveDual,
                               sup_powerset::PowerSet{<:Polyhedra.HalfSpace})
    JuMP.build_constraint(_error, ScaledPoint(-sup_powerset.set.a, sup_powerset.set.β), Sets.perspective_dual(subset))
end

## Polyhedron in Set ##
function set_space(space::Space, ::InclusionConstraint{<:Polyhedra.Rep,
                                                       <:VariableRef})
    return set_space(space, PrimalSpace)
end

function JuMP.add_constraint(model::JuMP.Model,
                             constraint::InclusionConstraint{<:Polyhedra.Rep,
                                                             <:Sets.AbstractSet{JuMP.VariableRef}},
                             name::String = "")
    □ = constraint.subset
    ◯ = constraint.supset
    for line in lines(□)
        @constraint(model, line in ◯)
    end
    for ray in rays(□)
        @constraint(model, ray in ◯)
    end
    for point in points(□)
        @constraint(model, point in ◯)
    end
end

function JuMP.build_constraint(_error::Function,
                               subset::Sets.EllipsoidAtOrigin,
                               sup_powerset::PowerSet{<:HalfSpace})
    hs = sup_powerset.set
    P = [hs.β hs.a'
         hs.a subset.Q]
    JuMP.build_constraint(_error, Symmetric(P), PSDCone())
end

#     [p(x) ≤ 1] ⊆ [⟨a, x⟩ ≤ β]
# <= λ(1 - p(x)) ≤ β - ⟨a, x⟩ (necessary if p is SOS-convex)
#              0 ≤ λ(p(x) - 1) - ⟨a, x⟩ + β
# Set `x = 0`: 0 ≤ λ(p(0) - 1) + β
# hence if p(0) ≤ 1 (i.e. 0 ∈ S), then λ = β / (1 - p(0))
# Homogeneous case: λ = β
# Use build_constraint when SumOfSquares#66 if λ = β (e.g. homogeneous)
function JuMP.add_constraint(model::JuMP.Model,
                             constraint::InclusionConstraint{<:Union{Sets.ConvexPolynomialSublevelSetAtOrigin,
                                                                     Sets.PerspectiveConvexPolynomialSet},
                                                             <:HalfSpace},
                            name::String = "")
    p = Sets.gauge1(constraint.subset)
    λ = @variable(model, lower_bound=0.0)
    x = Sets.space_variables(constraint.subset)
    hs = dot(constraint.supset.a, x) - constraint.supset.β
    @constraint(model, λ * (p - 1) - hs in SOSCone())
end
