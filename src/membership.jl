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
    return JuMP.build_constraint(error, c.member, variablify(c.set))
end

JuMP.function_string(print_mode, c::MembershipConstraint) = string(c.member)
JuMP.in_set_string(print_mode, c::MembershipConstraint) = string(JuMP.math_symbol(print_mode, :in), c.set)
function JuMP.build_constraint(_error::Function, member,
                               set::Sets.AbstractSet)
    MembershipConstraint(member, set)
end

function JuMP.build_constraint(_error::Function,
                               member::Point,
                               set::Sets.PolarOf{<:Sets.EllipsoidAtOrigin})
    p = member
    P = [scaling(p) coord(p)'
         coord(p)   set.set.Q]
    JuMP.build_constraint(_error, Symmetric(P), PSDCone())
end
function JuMP.build_constraint(_error::Function,
                               member::Point,
                               set::Sets.EllipsoidAtOrigin{T}) where T <: Number
    # The eltype of Point is an expression of JuMP variables so we cannot compute
    # x' * Q * x <= 1, we need to do
    # [ 1 x'     ]
    # [ x Q^{-1} ] ⪰ 0 so we switch to the polar representation
    JuMP.build_constraint(_error, member, Sets.polar_representation(set))
end
function JuMP.build_constraint(_error::Function,
                               member::Point{T},
                               set::Union{Sets.EllipsoidAtOrigin,
                                          Sets.ConvexPolynomialSublevelSetAtOrigin}) where T <: Number
    p = member
    JuMP.build_constraint(_error, sublevel_eval(set, coord(p)), MOI.LessThan(scaling(p)^2))
end
function JuMP.build_constraint(_error::Function,
                               member::Polyhedra.Line,
                               set::Union{Sets.EllipsoidAtOrigin,
                                          Sets.ConvexPolynomialSublevelSetAtOrigin})
    # We must have (λl)^T Q (λl) ≤ 1 for all λ hence we must have l^T Q l ≤ 0
    # As Q is positive definite, it means l^T Q l = 0
    l = Polyhedra.coord(member)
    JuMP.build_constraint(_error, sublevel_eval(set, l), MOI.EqualTo(0.0))
end
function JuMP.build_constraint(_error::Function,
                               member::Polyhedra.Ray,
                               set::Union{Sets.EllipsoidAtOrigin,
                                          Sets.ConvexPolynomialSublevelSetAtOrigin})
    # We must have (λl)^T Q (λl) ≤ 1 for all λ > 0 hence we must have l^T Q l ≤ 0
    # As Q is positive definite, it means l^T Q l = 0
    r = Polyhedra.coord(member)
    JuMP.build_constraint(_error, sublevel_eval(set, r), MOI.EqualTo(0.0))
end

function JuMP.build_constraint(_error::Function,
                               member::Point{T},
                               set::Union{Sets.PerspectiveEllipsoid,
                                          Sets.PerspectiveConvexPolynomialSet}) where {T <: Number}
    p = member
    val = sublevel_eval(set, coord(p), scaling(p))
    JuMP.build_constraint(_error, val, MOI.LessThan(0.0))
end
function JuMP.build_constraint(_error::Function,
                               member::SymScaledPoint{T},
                               set::Union{Sets.PerspectiveEllipsoid,
                                          Sets.PerspectiveConvexPolynomialSet}) where {T <: Number}
    p = member
    val = sublevel_eval(set, coord(p), scaling(p))
    JuMP.build_constraint(_error, val, MOI.EqualThan(0.0))
end
