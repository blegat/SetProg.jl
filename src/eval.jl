function sublevel_eval(set::Sets.ShiftedEllipsoid, a::AbstractVector)
    return quad_form(set.Q, a) + 2dot(set.b, a) + set.β + 1
end

function sublevel_eval(ell::Sets.Ellipsoid, a::AbstractVector)
    return quad_form(ell.Q, a)
end
function sublevel_eval(set::Sets.Householder, a::AbstractVector, β)
    x = Sets.space_variables(set)
    z = Sets.perspective_variable(set)
    # Avoid large values, with high degree polynomials, it might cause issues
    if false
        scale_a = norm(a)
        scale_β = abs(β)
        scaling = max(scale_a, scale_β) / sqrt(scale_a * scale_β)
    else
        scaling = 1.0
    end
    return Sets.perspective_gauge0(set)(z => β / scaling, x => a / scaling)
end

function sublevel_eval(
    set::Sets.PolarOf{<:Sets.Piecewise{<:Any,<:Sets.PolarPoint}},
    x::AbstractVector,
)
    inner = polar(set)
    return maximum(
        dot(x, Polyhedra.coord(ray)) / dot(piece.a, Polyhedra.coord(ray)) for
        (piece, cone) in zip(inner.sets, inner.pieces) for ray in rays(cone)
    )
end

function sublevel_eval(set::Sets.PolarPoint, x::AbstractVector)
    return dot(set.a, x)
end
function sublevel_eval(set::Sets.Piecewise{<:Any,<:Sets.PolarPoint}, x::AbstractVector)
    return maximum(sublevel_eval(piece, x) for piece in set.sets)
end

# Unlike sublevel_eval, this always evaluates a power of the actual gauge.
# The exponent of the operation need not equal the defining polynomial degree.
function _gauge_power(set::Sets.Ellipsoid, x, p)
    return sublevel_eval(set, x)^(p / 2)
end
function _gauge_power(set::Sets.PolarOf{<:Sets.Ellipsoid}, x, p)
    return dot(x, polar(set).Q \ x)^(p / 2)
end
function _gauge_power(set::Sets.PolarPoint, x, p)
    return max(zero(eltype(x)), sublevel_eval(set, x))^p
end
function _gauge_power(set::Sets.Piecewise, x, p)
    i = findfirst(piece -> x in piece, set.pieces)
    return _gauge_power(set.sets[i], x, p)
end
function _gauge_power(set::Sets.PolarOf{<:Sets.Piecewise{<:Any,<:Sets.PolarPoint}}, x, p)
    return sublevel_eval(set, x)^p
end
function sublevel_eval(set::Sets.PolarLpSum, x::AbstractVector)
    return _gauge_power(set.first, x, set.p) + _gauge_power(set.second, x, set.p)
end
function _gauge_power(set::Sets.PolarLpSum, x, p)
    return sublevel_eval(set, x)^(p / set.p)
end
function Polyhedra.support_function(x::AbstractVector, set::Sets.LpSum)
    return (
        _gauge_power(polar(set.first), x, set.p) +
        _gauge_power(polar(set.second), x, set.p)
    )^(1 / set.p)
end
function _gauge_power(set::Sets.Intersection, x, p)
    return x in set.domain ? _gauge_power(set.set, x, p) : Inf
end
