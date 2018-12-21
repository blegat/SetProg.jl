include("eval.jl")

### MembershipConstraint ###
struct ScaledPoint{T, S, P<:AbstractVector{T}}
    coord::P
    scaling::S
end
Base.copy(p::ScaledPoint) = ScaledPoint(p.coord, p.scaling)

const Point{T} = Union{AbstractVector{T}, ScaledPoint{T}}
Polyhedra.coord(p::ScaledPoint) = p.coord
scaling(::AbstractVector) = 1.0
scaling(p::ScaledPoint) = p.scaling

"""
    struct MembershipConstraint{P, S} <: SetConstraint
        member::P
        set::S
    end

Constrain `member in set`.
"""
struct MembershipConstraint{P, S} <: SetConstraint
    member::P
    set::S
end
need_variablify(c::MembershipConstraint) = need_variablify(c.set)
function variablify(c::MembershipConstraint)
    return MembershipConstraint(c.member, variablify(c.set))
end

JuMP.function_string(print_mode, c::MembershipConstraint) = string(c.member)
JuMP.in_set_string(print_mode, c::MembershipConstraint) = string(JuMP.math_symbol(print_mode, :in), c.set)
function JuMP.build_constraint(_error::Function, member,
                               set::Sets.AbstractSet)
    MembershipConstraint(member, set)
end

function JuMP.add_constraint(model::JuMP.Model,
                             constraint::MembershipConstraint{<:Point,
                                                              <:Sets.PolarEllipsoidAtOrigin},
                             name::String = "")
    p = constraint.member
    P = [scaling(p) coord(p)'
         coord(p)   constraint.set.Q]
    @constraint(model, Symmetric(P) in PSDCone())
end
function JuMP.add_constraint(model::JuMP.Model,
                             constraint::MembershipConstraint{<:Point{T},
                                                              <:Union{Sets.EllipsoidAtOrigin,
                                                                      Sets.ConvexPolynomialSublevelSetAtOrigin}},
                             name::String = "") where {T <: Number}
    p = constraint.member
    @constraint(model, sublevel_eval(constraint.set, coord(p)) <= scaling(p)^2)
end
function JuMP.add_constraint(model::JuMP.Model,
                             constraint::MembershipConstraint{<:Polyhedra.Line,
                                                              <:Union{Sets.EllipsoidAtOrigin,
                                                                      Sets.ConvexPolynomialSublevelSetAtOrigin}},
                             name::String = "")
    # We must have (λl)^T Q (λl) ≤ 1 for all λ hence we must have l^T Q l ≤ 0
    # As Q is positive definite, it means l^T Q l = 0
    l = Polyhedra.coord(constraint.member)
    @constraint(model, sublevel_eval(constraint.set, l) in MOI.EqualTo(0.0))
end
function JuMP.add_constraint(model::JuMP.Model,
                             constraint::MembershipConstraint{<:Polyhedra.Ray,
                                                              <:Union{Sets.EllipsoidAtOrigin,
                                                                      Sets.ConvexPolynomialSublevelSetAtOrigin}},
                             name::String = "")
    # We must have (λl)^T Q (λl) ≤ 1 for all λ > 0 hence we must have l^T Q l ≤ 0
    # As Q is positive definite, it means l^T Q l = 0
    r = Polyhedra.coord(constraint.member)
    @constraint(model, sublevel_eval(constraint.set, r) in MOI.EqualTo(0.0))
end

