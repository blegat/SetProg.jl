### InclusionConstraint ###
struct InclusionConstraint{SubSetType, SupSetType} <: SetConstraint
    subset::SubSetType
    supset::SupSetType
end
function need_variablify(c::InclusionConstraint)
    return need_variablify(c.subset) || need_variablify(c.supset)
end
function variablify(c::InclusionConstraint)
    return InclusionConstraint(variablify(c.subset), variablify(c.supset))
end

JuMP.function_string(print_mode, c::InclusionConstraint) = string(c.subset)
function JuMP.in_set_string(print_mode, c::InclusionConstraint)
    string(print_mode == JuMP.IJuliaMode ? "\\subseteq" : "⊆", " ", c.supset)
end
struct PowerSet{S}
    set::S
end
function JuMP.build_constraint(_error::Function, subset,
                               supset_powerset::PowerSet)
    InclusionConstraint(subset, supset_powerset.set)
end

# Primal:
#   set : x^T Q x ≤ 1
#   volume proportional to 1/det(Q)
#   t ≤ det(Q)^(1/n) <=> 1/t^n ≥ 1/det(Q)
#   volume proportional to 1/t^n
#   For t ≤ det(Q)^(1/n) to be tight we need to maximize `t`
#   hence we need to minimize the volume
function set_space(space::Space, ::InclusionConstraint{<:VariableRef, <:Polyhedra.Rep})
    return set_space(space, DualSpace)
end
function set_space(space::Space, ::InclusionConstraint{<:Polyhedra.Rep, <:VariableRef})
    return set_space(space, PrimalSpace)
end

### InclusionConstraint for sets ###

## Set in Set ##


## Set in Polyhedron ##
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
# TODO if a is not constant, use Schur Lemma
function quad_form(Q::Symmetric{<:JuMP.AbstractJuMPScalar},
                   a::AbstractVector{<:AbstractMonomialLike})
    n = length(a)
    @assert n == LinearAlgebra.checksquare(Q)
    return sum((i == j ? 1 : 2) * a[i] * Q[i, j] * a[j] for j in 1:n for i in 1:j)
end
function quad_form(Q::Symmetric{JuMP.VariableRef}, a::AbstractVector{<:Real})
    n = length(a)
    @assert n == LinearAlgebra.checksquare(Q)
    aff = zero(JuMP.GenericAffExpr{eltype(a), JuMP.VariableRef})
    for j in 1:n
        for i in 1:n
            JuMP.add_to_expression!(aff, a[i] * a[j], Q[i, j])
        end
    end
    return aff
end
function poly_eval(p::AbstractPolynomial{JuMP.AffExpr},
                   a::AbstractVector{Float64})
    vars = variables(p)
    aff = zero(JuMP.AffExpr)
    for term in terms(p)
        mono = monomial(term)
        aff = JuMP.destructive_add!(aff, mono(vars => a), coefficient(term))
    end
    return aff
end
function sublevel_eval(ell::Union{Sets.EllipsoidAtOrigin,
                                  Sets.PolarEllipsoidAtOrigin},
                       a::AbstractVector)
    return quad_form(ell.Q, a)
end
function sublevel_eval(set::Union{Sets.PolynomialSublevelSetAtOrigin,
                                  Sets.PolarPolynomialSublevelSetAtOrigin},
                       a::AbstractVector)
    return poly_eval(polynomial(set.p), a)
end
function sublevel_eval(model, set::Sets.DualQuadCone, a::AbstractVector, β)
    d = data(model)
    x = d.polyvars[1:dimension(constraint.subset)]
    z = d.perspective_polyvar
    return set.p(z => -β, x => a)
end
function JuMP.add_constraint(model::JuMP.Model,
                             constraint::InclusionConstraint{<:Union{Sets.PolarEllipsoidAtOrigin{JuMP.VariableRef},
                                                                     Sets.PolarPolynomialSublevelSetAtOrigin{JuMP.VariableRef}},
                                                             <:Polyhedra.HyperPlane},
                             name::String = "")
    @assert iszero(h.set.β) # Otherwise it is not symmetric around the origin
    @constraint(model, sublevel_eval(constraint.subset, constraint.supset.a) in MOI.EqualTo(constraint.supset.β^2))
end
function JuMP.add_constraint(model::JuMP.Model,
                             constraint::InclusionConstraint{<:Sets.DualQuadCone,
                                                             <:Polyhedra.HyperPlane},
                             name::String = "")
    val = sublevel_eval(model, constraint.subset, constraint.supset.a,
                        constraint.supset.β)
    @constraint(model, val in MOI.EqualTo(0.0))
end
function JuMP.add_constraint(model::JuMP.Model,
                             constraint::InclusionConstraint{<:Union{Sets.PolarEllipsoidAtOrigin{JuMP.VariableRef},
                                                                     Sets.PolarPolynomialSublevelSetAtOrigin{JuMP.VariableRef}},
                                                             <:Polyhedra.HalfSpace},
                             name::String = "")
    @constraint(model, sublevel_eval(constraint.subset, constraint.supset.a) in MOI.LessThan(constraint.supset.β^2))
end
function JuMP.add_constraint(model::JuMP.Model,
                             constraint::InclusionConstraint{<:Sets.DualQuadCone,
                                                             <:Polyhedra.HalfSpace},
                             name::String = "")
    val = sublevel_eval(model, constraint.subset, constraint.supset.a,
                        constraint.supset.β)
    @constraint(model, val in MOI.LessThan(0.0))
end

## Polyhedron in Set ##
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
