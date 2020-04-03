need_variablify(lm::Sets.LinearImage) = need_variablify(lm.set)
function variablify(lm::Sets.LinearImage)
    return Sets.LinearImage(variablify(lm.set), lm.A)
end
need_variablify(lm::Sets.AffineImage) = need_variablify(lm.set)
function variablify(lm::Sets.AffineImage)
    return Sets.AffineImage(variablify(lm.set), lm.A, lm.b)
end


"""
    apply_map(li::Sets.LinearImage{<:Sets.PolarOf{<:Sets.EllipsoidAtOrigin}}, new_vars)

The set ``(AS)^\\circ``, the polar of the set ``AS``, is ``A^{-\\top}S^\\circ``
and given ``S^\\circ = \\{\\, x \\mid x^\\top Q x \\le 1\\,\\}``, we have
``A^{-\\top}S^\\circ = \\{\\, x \\mid x^\\top AQA^\\top x \\le 1\\,\\}``.
"""
function apply_map(li::Sets.LinearImage{<:Sets.PolarOf{<:Sets.EllipsoidAtOrigin}}, new_vars)
    return Sets.polar(Sets.EllipsoidAtOrigin(Symmetric(li.A * Sets.polar(li.set).Q * li.A')))
end

function apply_map(li::Sets.LinearPreImage{<:Sets.EllipsoidAtOrigin}, new_vars)
    return Sets.EllipsoidAtOrigin(Symmetric(li.A' * li.set.Q * li.A))
end

"""
    apply_map(li::Sets.LinearImage{<:Sets.PolarOf{<:Sets.ConvexPolynomialSublevelSetAtOrigin}}, new_vars)

The set ``(AS)^\\circ``, the polar of the set ``AS``, is ``A^{-\\top}S^\\circ``
and given ``S^\\circ = \\{\\, x \\mid p(x) \\le 1\\,\\}``, we have
``A^{-\\top}S^\\circ = \\{\\, x \\mid x^\\top p(A^\\top x) \\le 1\\,\\}``.
"""
function apply_map(li::Sets.LinearImage{<:Sets.PolarOf{<:Sets.ConvexPolynomialSublevelSetAtOrigin}}, new_vars)
    deg = Sets.polar(li.set).degree
    @assert iseven(deg)
    q = apply_matrix(Sets.polar(li.set).p, li.A', new_vars,
                     div(deg, 2))
    return Sets.polar(Sets.ConvexPolynomialSublevelSetAtOrigin(deg, q, nothing))
end

function apply_map(li::Sets.LinearPreImage{<:Sets.PolynomialSublevelSetAtOrigin},
                   new_vars)
    deg = li.set.degree
    @assert iseven(deg)
    q = apply_matrix(li.set.p, li.A, new_vars, div(deg, 2))
    return Sets.PolynomialSublevelSetAtOrigin(deg, q)
end
function apply_map(li::Sets.LinearPreImage{<:Sets.ConvexPolynomialSublevelSetAtOrigin},
                   new_vars)
    deg = li.set.degree
    @assert iseven(deg)
    q = apply_matrix(li.set.p, li.A, new_vars, div(deg, 2))
    return Sets.ConvexPolynomialSublevelSetAtOrigin(deg, q, nothing)
end

"""
    apply_map(li::Sets.LinearImage{<:Sets.PerspectiveDualOf{<:Union{Sets.PerspectiveEllipsoid,
                                                                    Sets.PerspectiveConvexPolynomialSet}}}, new_vars)

The set ``(AS)^\\circ``, the polar of the set ``AS``, is ``A^{-\\top}S^\\circ``
and given ..., we have
...
"""
function apply_map(li::Sets.LinearImage{<:Sets.PerspectiveDualOf{<:Sets.Householder{T}}}, new_vars) where T
    old_vars = Sets.space_variables(li.set)
    p = subs(Sets.perspective_gauge0(li.set.set), old_vars => li.A' * new_vars)
    dual = Sets.Householder(Sets.UnknownSet{T}(), p, li.A * li.set.set.h,
                            Sets.perspective_variable(li.set), new_vars)
    return Sets.perspective_dual(dual)
end

# [x'] = [A b] [x]
# [z'] = [0 1] [z]
# The transpose is
# [A' 0]
# [b' 1]
# So we substitute x for A' x' and then z for b' x' + z'
function apply_map(li::Sets.AffineImage{<:Sets.PerspectiveDualOf{<:Sets.Householder{T}}},
                   new_vars) where T
    z = Sets.perspective_variable(li.set)
    old_vars = Sets.space_variables(li.set)
    q = subs(Sets.perspective_gauge0(li.set.set), old_vars => li.A' * new_vars)
    # We do this substitution after in case `new_vars` and `old_vars` have an intersection
    p = subs(q, z => z + dot(li.b, new_vars))
    dual = Sets.Householder(Sets.UnknownSet{T}(), p, li.A * li.set.set.h + li.b,
                            z, new_vars)
    return Sets.perspective_dual(dual)
end

# FIXME, for Sets.AbstractSet, we should apply it directly
function Base.:(*)(A::AbstractMatrix, set::Union{SetVariableRef, Sets.AbstractSet})
    return Sets.LinearImage(set, A)
end
function Base.:(+)(li::Sets.LinearImage, b::AbstractVector)
    return Sets.AffineImage(li.set, li.A, b)
end
