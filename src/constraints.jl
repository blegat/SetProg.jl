abstract type SetConstraint  <: JuMP.AbstractConstraint end

### MembershipConstraint ###
struct MembershipConstraint{P, S}
    member::P
    set::S
end
JuMP.function_string(print_mode, c::MembershipConstraint) = string(c.member)
JuMP.in_set_string(print_mode, c::MembershipConstraint) = string(JuMP.math_symbol(print_mode, :in), c.set)
function JuMP.build_constraint(_error::Function, member,
                               set::Sets.AbstractSet)
    MembershipConstraint(member, set)
end

### InclusionConstraint ###
struct InclusionConstraint{SubSetType, SupSetType} <: SetConstraint
    subset::SubSetType
    supset::SupSetType
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
function JuMP.add_constraint(model::JuMP.Model,
                             constraint::InclusionConstraint{<:Union{Sets.PolarEllipsoidAtOrigin{JuMP.VariableRef},
                                                                     Sets.PolarPolynomialSublevelSetAtOrigin{JuMP.VariableRef}},
                                                             <:Polyhedra.HyperPlane},
                             name::String = "")
    @assert iszero(h.set.β) # Otherwise it is not symmetric around the origin
    @constraint(model, sublevel_eval(constraint.subset, constraint.supset.a) in MOI.EqualTo(constraint.supset.β^2))
end
function JuMP.add_constraint(model::JuMP.Model,
                             constraint::InclusionConstraint{<:Union{Sets.PolarEllipsoidAtOrigin{JuMP.VariableRef},
                                                                     Sets.PolarPolynomialSublevelSetAtOrigin{JuMP.VariableRef}},
                                                             <:Polyhedra.HalfSpace},
                             name::String = "")
    @constraint(model, sublevel_eval(constraint.subset, constraint.supset.a) in MOI.LessThan(constraint.supset.β^2))
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
function JuMP.add_constraint(model::JuMP.Model,
                             constraint::MembershipConstraint{<:AbstractVector,
                                                              <:Union{Sets.EllipsoidAtOrigin,
                                                                      Sets.PolynomialSublevelSetAtOrigin}},
                             name::String = "")
    @constraint(model, sublevel_eval(constraint.set, constraint.member) in MOI.LessThan(1.0))
end
function JuMP.add_constraint(model::JuMP.Model,
                             constraint::MembershipConstraint{<:Polyhedra.Line,
                                                              <:Union{Sets.EllipsoidAtOrigin,
                                                                      Sets.PolynomialSublevelSetAtOrigin}},
                             name::String = "")
    # We must have (λl)^T Q (λl) ≤ 1 for all λ hence we must have l^T Q l ≤ 0
    # As Q is positive definite, it means l^T Q l = 0
    l = Polyhedra.coord(constraint.member)
    @constraint(model, sublevel_eval(constraint.set, l) in MOI.EqualTo(0.0))
end
function JuMP.add_constraint(model::JuMP.Model,
                             constraint::MembershipConstraint{<:Polyhedra.Ray,
                                                              <:Union{Sets.EllipsoidAtOrigin,
                                                                      Sets.PolynomialSublevelSetAtOrigin}},
                             name::String = "")
    # We must have (λl)^T Q (λl) ≤ 1 for all λ > 0 hence we must have l^T Q l ≤ 0
    # As Q is positive definite, it means l^T Q l = 0
    r = Polyhedra.coord(constraint.member)
    @constraint(model, sublevel_eval(constraint.set, r) in MOI.EqualTo(0.0))
end

### InclusionConstraint for variable sets  ###

## Loading  ##

function load(model::JuMP.Model, constraint::SetConstraint)
    return JuMP.add_constraint(model, variablify(constraint))
end

variablify(p::Polyhedra.Rep) = p
variablify(v::VariableRef) = v.variable

function variablify(c::InclusionConstraint)
    return InclusionConstraint(variablify(c.subset), variablify(c.supset))
end

## Adding ##

struct SetShape <: JuMP.AbstractShape end
struct ConstraintIndex
    value::Int
end

function JuMP.add_constraint(model::JuMP.Model, constraint::InclusionConstraint,
                             name::String)
    @assert constraint.subset isa VariableRef || constraint.supset isa VariableRef
    container = data(model)
    index = ConstraintIndex(container.last_index += 1)
    container.constraints[index] = constraint
    container.names[index] = name
    return JuMP.ConstraintRef(model, index, SetShape())
end

### SetConstraintRef ###

const SetConstraintRef{M} = JuMP.ConstraintRef{M, ConstraintIndex, SetShape}
function JuMP.name(cref::SetConstraintRef)
    return data(cref.model).names[cref.index]
end
function JuMP.constraint_object(cref::SetConstraintRef)
    return data(cref.model).constraints[cref.index]
end
