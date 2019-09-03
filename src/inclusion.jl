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
    return merge_spaces(spaces, sub, sup)
end

JuMP.function_string(print_mode, c::InclusionConstraint) = string(c.subset)
function JuMP.in_set_string(print_mode, c::InclusionConstraint)
    string(print_mode == JuMP.IJuliaMode ? "\\subseteq" : "⊆", " ", c.supset)
end
struct PowerSet{S}
    set::S
end

# Fallback, might be because `subset` or `sup_powerset` is a `SetVariableRef` or
# a `Polyhedron` (which is handled by `JuMP.add_constraint`).
function JuMP.build_constraint(_error::Function, subset, sup_powerset::PowerSet;
                               kws...)
    InclusionConstraint(subset, sup_powerset.set, kws)
end

### InclusionConstraint for sets ###

## Set in Set ##
# See [LTJ18] B. Legat, P. Tabuada and R. M. Jungers.
# *Computing controlled invariant sets for hybrid systems with applications to model-predictive control*.
# 6th IFAC Conference on Analysis and Design of Hybrid Systems ADHS 2018, **2018**.
function set_space(space::Space, ::InclusionConstraint{<:Sets.LinearImage,
                                                       <:Sets.LinearImage})
    return set_space(space, DualSpace)
end
# AS ⊆ S <=> S ⊆ A^{-1}S so PrimalSpace work
# AS ⊆ S <=> A^{-T}S∘ ⊆ S∘ so DualSpace work
function set_space(space::Space, ::InclusionConstraint{<:Sets.LinearImage,
                                                       <:SetVariableRef})
    return space
end
# We can always transform  an ellipsoid to primal or dual space so we can handle
# any space
function set_space(space::Space,
                   ::InclusionConstraint{<:SetVariableRef,
                                         <:Sets.AbstractEllipsoid{T}}) where T<:Number
    return space
end

# S-procedure: Q ⊆ P <=> xQx ≤ 1 => xPx ≤ 1 <=> xPx ≤ xQx <=> Q - P is PSD
function JuMP.build_constraint(_error::Function,
                               subset::Sets.EllipsoidAtOrigin,
                               sup_powerset::PowerSet{<:Sets.EllipsoidAtOrigin};
                               S_procedure_scaling = nothing)
    @assert S_procedure_scaling === nothing || isone(S_procedure_scaling)
    Q = subset.Q
    P = sup_powerset.set.Q
    return psd_constraint(Symmetric(Q - P))
end

# S-procedure: Q ⊆ P <=> q(x) ≤ 1 => p(x) ≤ 1 <=> p(x) ≤ q(x) <= q - p is SOS
function JuMP.build_constraint(_error::Function,
                               subset::Union{Sets.PolynomialSublevelSetAtOrigin,
                                             Sets.ConvexPolynomialSublevelSetAtOrigin},
                               sup_powerset::PowerSet{<:Union{Sets.PolynomialSublevelSetAtOrigin,
                                                              Sets.ConvexPolynomialSublevelSetAtOrigin}};
                               S_procedure_scaling = nothing)
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
                             constraint::InclusionConstraint{<:Sets.Householder,
                                                             <:Sets.Householder},
                             name::String = "")
    return JuMP.add_constraint(model,
                               s_procedure(model, constraint.subset,
                                           constraint.supset;
                                           constraint.kws...), name)
end

# S ⊆ T <=> T* ⊇ S*
function JuMP.build_constraint(_error::Function,
                               subset::Sets.HouseDualOf,
                               sup_powerset::PowerSet{<:Sets.HouseDualOf};
                               kws...)
    S = subset
    T = sup_powerset.set
    JuMP.build_constraint(_error, Sets.perspective_dual(T), PowerSet(Sets.perspective_dual(S)); kws...)
end

# S ⊆ T <=> polar(T) ⊆ polar(S)
function JuMP.build_constraint(_error::Function,
                               subset::Sets.Polar,
                               sup_powerset::PowerSet{<:Sets.Polar}; kws...)
    S = subset
    T = sup_powerset.set
    JuMP.build_constraint(_error, Sets.polar(T), PowerSet(Sets.polar(S));
                          kws...)
end

# See [LTJ18]
function JuMP.build_constraint(_error::Function,
                               subset::Sets.LinearImage{<:Union{Sets.Polar, Sets.PerspectiveDual}},
                               sup_powerset::PowerSet{<:Sets.LinearImage{<:Union{Sets.Polar, Sets.PerspectiveDual}}};
                               kws...)
    dim = Sets.dimension(subset)
    @polyvar x[1:dim]
    JuMP.build_constraint(_error, apply_map(subset, x),
                          PowerSet(apply_map(sup_powerset.set, x)); kws...)
end
function JuMP.build_constraint(_error::Function,
                               subset::Sets.AbstractSet,
                               sup_powerset::PowerSet{<:Sets.LinearPreImage{<:Sets.AbstractSet}};
                               kws...)
    x = Sets.space_variables(subset)
    JuMP.build_constraint(_error, subset,
                          PowerSet(apply_map(sup_powerset.set, x)); kws...)
end
function JuMP.build_constraint(_error::Function,
                               subset::Sets.LinearImage{<:Sets.AbstractSet},
                               sup_powerset::PowerSet{<:Sets.AbstractSet};
                               kws...)
    JuMP.build_constraint(
        _error, subset.set,
        PowerSet(Sets.LinearPreImage(sup_powerset.set, subset.A)); kws...)
end
function JuMP.build_constraint(_error::Function,
                               subset::Sets.LinearImage{<:Union{Sets.Polar, Sets.PerspectiveDual}},
                               sup_powerset::PowerSet{<:Sets.AbstractSet};
                               kws...)
    x = Sets.space_variables(sup_powerset.set)
    JuMP.build_constraint(_error, apply_map(subset, x), sup_powerset; kws...)
end

## Set in Polyhedron ##
function set_space(space::Space, ::InclusionConstraint{<:SetVariableRef,
                                                       <:Polyhedra.Rep})
    return set_space(space, DualSpace)
end

# Ellipsoid #
function JuMP.add_constraint(model::JuMP.Model,
                             constraint::InclusionConstraint{<:Sets.AbstractSet{<:JuMP.AbstractJuMPScalar},
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
                                                       <:SetVariableRef})
    return set_space(space, PrimalSpace)
end

function JuMP.add_constraint(model::JuMP.Model,
                             constraint::InclusionConstraint{<:Polyhedra.Rep,
                                                             <:Sets.AbstractSet{<:JuMP.AbstractJuMPScalar}},
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
    return psd_constraint(Symmetric(P))
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
                                                                     Sets.ConvexPolynomialSet},
                                                             <:HalfSpace},
                            name::String = "")
    p = Sets.gauge1(constraint.subset)
    λ = @variable(model, lower_bound=0.0, base_name = "λ")
    x = Sets.space_variables(constraint.subset)
    hs = dot(constraint.supset.a, x) - constraint.supset.β
    cref = @constraint(model, λ * (p - 1) - hs in SOSCone())
    return cref
end
function JuMP.build_constraint(_error::Function,
                               subset::Sets.Householder,
                               sup_powerset::PowerSet{<:HalfSpace})
    # 0 ≤ βz + ⟨-a, x⟩
    x = [sup_powerset.set.β; -sup_powerset.set.a]
    H = Sets._householder(subset.h)
    # The householder transformation is symmetric and orthogonal so no need to
    # worry about whether we should invert or transpose H
    y = H * x
    JuMP.build_constraint(_error, subset.set,
                          PowerSet(HalfSpace(y[2:end], -y[1])))
end
