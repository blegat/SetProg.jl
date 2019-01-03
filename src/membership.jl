include("eval.jl")

### MembershipConstraint ###
struct ScaledPoint{T, S, P<:AbstractVector{T}}
    coord::P
    scaling::S
end
Base.copy(p::ScaledPoint) = ScaledPoint(p.coord, p.scaling)
struct SymScaledPoint{T, S, P<:AbstractVector{T}}
    coord::P
    scaling::S
end
Base.copy(p::SymScaledPoint) = SymScaledPoint(p.coord, p.scaling)

const Point{T} = Union{AbstractVector{T}, ScaledPoint{T}}
Polyhedra.coord(p::Union{ScaledPoint, SymScaledPoint}) = p.coord
scaling(::AbstractVector) = 1.0
scaling(p::Union{ScaledPoint, SymScaledPoint}) = p.scaling

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
                                                              <:Sets.PolarOf{<:Sets.EllipsoidAtOrigin}},
                             name::String = "")
    p = constraint.member
    P = [scaling(p) coord(p)'
         coord(p)   constraint.set.set.Q]
    @constraint(model, Symmetric(P) in PSDCone())
end
function JuMP.add_constraint(model::JuMP.Model,
                             constraint::MembershipConstraint{<:Point,
                                                              Sets.EllipsoidAtOrigin{T}},
                             name::String = "") where {T <: Number}
    # The eltype of Point is an expression of JuMP variables so we cannot compute
    # x' * Q * x <= 1, we need to do
    # [ 1 x'     ]
    # [ x Q^{-1} ] ⪰ 0 so we switch to the polar representation
    @constraint(model, constraint.member in Sets.polar_representation(constraint.set))
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

function JuMP.add_constraint(model::JuMP.Model,
                             constraint::MembershipConstraint{<:Point{T},
                                                              <:Union{Sets.PerspectiveEllipsoid,
                                                                      Sets.PerspectiveConvexPolynomialSet}},
                             name::String = "") where {T <: Number}
    p = constraint.member
    val = sublevel_eval(model, constraint.set, coord(p), scaling(p))
    @constraint(model, val in MOI.LessThan(0.0))
end
function JuMP.add_constraint(model::JuMP.Model,
                             constraint::MembershipConstraint{<:SymScaledPoint{T},
                                                              <:Union{Sets.PerspectiveEllipsoid,
                                                                      Sets.PerspectiveConvexPolynomialSet}},
                             name::String = "") where {T <: Number}
    p = constraint.member
    val = sublevel_eval(model, constraint.set, coord(p), scaling(p))
    @constraint(model, val in MOI.EqualThan(0.0))
end
