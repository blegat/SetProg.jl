"""
    LpSum(first, second, p)

Firey Lₚ sum of two convex sets containing the origin, for finite `p ≥ 1`.
Its support function satisfies
`h(LpSum(C, S, p), x)^p = h(C, x)^p + h(S, x)^p`.
In particular, `p = 1` is Minkowski addition.

The exponent is independent of the representation of either operand.
See also [`PolarLpSum`](@ref).
"""
struct LpSum{T,A,B,P<:Real} <: AbstractSet{T}
    first::A
    second::B
    p::P
end

"""
    PolarLpSum(first, second, p)

Polar Firey Lₚ sum of two convex gauge sets, for finite `p ≥ 1`:
`polar(LpSum(polar(first), polar(second), p))`.
Its gauge satisfies
`g(PolarLpSum(C, S, p), x)^p = g(C, x)^p + g(S, x)^p`.

For example, `PolarLpSum(Ellipsoid(Q), Ellipsoid(R), 2)` represents the
ellipsoid with matrix `Q + R`. The operands here are concrete `Sets.Ellipsoid`
sets. An exponent of two also makes sense for polyhedral operands: their
squared gauges are added, rather than their gauges.

The operands must have convex, nonnegative gauges. This operation uses a
positive exponent; it corresponds to parameter `-p` in the convention of
Gardner, Hug and Weil, *Operations between sets in geometry*, Section 5.5.
"""
struct PolarLpSum{T,A,B,P<:Real} <: AbstractSet{T}
    first::A
    second::B
    p::P
end

const AbstractLpSum = Union{LpSum,PolarLpSum}
_lp_scalar_type(::AbstractSet{T}) where {T} = T
_lp_dimension(set::AbstractSet) = dimension(set)
function _lp_dimension(first, second)
    a, b = _lp_dimension(first), _lp_dimension(second)
    a === nothing && return b
    b === nothing && return a
    a == b || throw(DimensionMismatch("Lₚ sum operands have different dimensions"))
    return a
end
function _lp_variables(first, second)
    a, b = space_variables(first), space_variables(second)
    a === nothing && return b
    b === nothing && return a
    a == b || throw(ArgumentError("Lₚ sum operands must share polynomial variables"))
    return a
end
function _check_lp_sum(first, second, p)
    isfinite(p) && p >= 1 || throw(ArgumentError("the exponent must be finite and ≥ 1"))
    _lp_dimension(first, second)
    _lp_variables(first, second)
    return promote_type(_lp_scalar_type(first), _lp_scalar_type(second))
end
function LpSum(first, second, p::Real)
    T = _check_lp_sum(first, second, p)
    return LpSum{T,typeof(first),typeof(second),typeof(p)}(first, second, p)
end
function PolarLpSum(first, second, p::Real)
    T = _check_lp_sum(first, second, p)
    return PolarLpSum{T,typeof(first),typeof(second),typeof(p)}(first, second, p)
end

dimension(set::AbstractLpSum) = _lp_dimension(set.first, set.second)
space_variables(set::AbstractLpSum) = _lp_variables(set.first, set.second)
function perspective_variable(set::AbstractLpSum)
    a, b = perspective_variable(set.first), perspective_variable(set.second)
    a === nothing && return b
    b === nothing && return a
    a == b || throw(ArgumentError("Lₚ sum operands must share their perspective variable"))
    return a
end

# For concrete operands, distribute polarity through the expression. Expressions
# involving set variables retain the generic Polar wrapper until loading.
function Polyhedra.polar(set::LpSum{T,A,B}) where {T,A<:AbstractSet,B<:AbstractSet}
    return PolarLpSum(polar(set.first), polar(set.second), set.p)
end
function Polyhedra.polar(set::PolarLpSum{T,A,B}) where {T,A<:AbstractSet,B<:AbstractSet}
    return LpSum(polar(set.first), polar(set.second), set.p)
end
