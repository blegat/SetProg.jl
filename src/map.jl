need_variablify(p::Sets.Projection) = need_variablify(p.set)
Sets.perspective_variable(p::Sets.Projection) = Sets.perspective_variable(p.set)
clear_spaces(p::Sets.Projection) = clear_spaces(p.set)

need_variablify(lm::Sets.LinearImage) = need_variablify(lm.set)
function variablify(lm::Sets.LinearImage)
    return Sets.LinearImage(variablify(lm.set), lm.A)
end

"""
    apply_map(li::Sets.LinearImage{<:Sets.Polar}, new_vars)

The set ``(AS)^\\circ``, the polar of the set ``AS``, is ``A^{-\\top}S^\\circ``.
"""
function apply_map(li::Sets.LinearImage{<:Sets.Polar}, new_vars)
    return Sets.polar(apply_map(Sets.LinearPreImage(Sets.polar(li.set), li.A'), new_vars))
end

function apply_map(li::Sets.LinearPreImage{<:Sets.PolarPoint}, new_vars)
    return Sets.PolarPoint(li.A' * li.set.a)
end

"""
    apply_map(li::Sets.LinearPreImage{<:Sets.Ellipsoid}, new_vars)

Given ``S = \\{\\, x \\mid x^\\top Q x \\le 1\\,\\}``, we have
``A^{-1}S = \\{\\, x \\mid x^\\top A^\\top Q A x \\le 1\\,\\}``.
"""
function apply_map(li::Sets.LinearPreImage{<:Sets.Ellipsoid}, new_vars)
    return Sets.Ellipsoid(Symmetric(li.A' * li.set.Q * li.A))
end

function apply_map(li::Sets.LinearPreImage{<:Sets.Piecewise}, new_vars)
    return Sets.Piecewise(
        [apply_map(Sets.LinearPreImage(set, li.A), new_vars) for set in li.set.sets],
        li.A \ li.set.polytope,
        [li.A \ piece for piece in li.set.pieces],
        li.set.graph # /!\ FIXME The nij are now incorrect
    )
    return Sets.Ellipsoid(Symmetric(li.A' * li.set.Q * li.A))
end

"""
    apply_map(li::Sets.LinearPreImage{<:Sets.PolySet}, new_vars)

Given ``S = \\{\\, x \\mid p(x) \\le 1\\,\\}``, we have
``A^{-1}S = \\{\\, x \\mid p(Ax) \\le 1\\,\\}``.
"""
function apply_map(li::Sets.LinearPreImage{<:Sets.PolySet}, new_vars)
    deg = li.set.degree
    @assert iseven(deg)
    q = apply_matrix(li.set.p, li.A, new_vars, div(deg, 2))
    return Sets.PolySet(deg, q)
end
function apply_map(li::Sets.LinearPreImage{<:Sets.ConvexPolySet}, new_vars)
    deg = li.set.degree
    @assert iseven(deg)
    q = apply_matrix(li.set.p, li.A, new_vars, div(deg, 2))
    return Sets.ConvexPolySet(deg, q, nothing)
end

abstract type SymbolicVariable end

# FIXME, for Sets.AbstractSet, we should apply it directly
function Base.:(*)(A::AbstractMatrix, set::Union{SetVariableRef, Sets.AbstractSet, SymbolicVariable})
    return Sets.LinearImage(set, A)
end

needs_polyvars(::Sets.AbstractSet) = true
needs_polyvars(::Union{Sets.Ellipsoid,Sets.PolarPoint}) = false
needs_polyvars(set::Union{Sets.Polar,Sets.PerspectiveDual}) = needs_polyvars(set.set)
needs_polyvars(set::Sets.Piecewise) = any(needs_polyvars, set.sets)
