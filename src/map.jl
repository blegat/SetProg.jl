mutable struct LinearImage{S, T, MT <: AbstractMatrix{T}}
    set::S
    A::MT
    space_index::Union{Nothing, SpaceIndex}
end


need_variablify(lm::LinearImage) = need_variablify(lm.set)
function variablify(lm::LinearImage)
    return LinearImage(variablify(lm.set), lm.A, lm.space_index)
end

function clear_spaces(li::LinearImage)
    li.space_index = nothing
end
Sets.perspective_variable(li::LinearImage) = Sets.perspective_variable(li.set)
function create_spaces(li::LinearImage, spaces::Spaces)
    if li.space_index === nothing
        li.space_index = new_space(spaces, size(li.A, 1))
    end
    return li.space_index
end
space_index(li::LinearImage) = li.space_index

"""
    apply_map(li::LinearImage{Sets.PolarEllipsoidAtOrigin})

The set ``(AS)^\\circ``, the polar of the set ``AS``, is ``A^{-\\top}S^\\circ``
and given ``S^\\circ = \\{\\, x \\mid x^\\top Q x \\le 1\\,\\}``, we have
``A^{-\\top}S^\\circ = \\{\\, x \\mid x^\\top AQA^\\top x \\le 1\\,\\}``.
"""
function apply_map(model, li::LinearImage{<:Sets.PolarEllipsoidAtOrigin})
    return Sets.PolarEllipsoidAtOrigin(Symmetric(li.A * li.set.Q * li.A'))
end

"""
    apply_map(li::LinearImage{Sets.PolarConvexPolynomialSublevelSetAtOrigin})

The set ``(AS)^\\circ``, the polar of the set ``AS``, is ``A^{-\\top}S^\\circ``
and given ``S^\\circ = \\{\\, x \\mid p(x) \\le 1\\,\\}``, we have
``A^{-\\top}S^\\circ = \\{\\, x \\mid x^\\top p(A^\\top x) \\le 1\\,\\}``.
"""
function apply_map(model,
                   li::LinearImage{<:Sets.PolarConvexPolynomialSublevelSetAtOrigin})
    d = data(model)
    new_vars = space_polyvars(d.spaces, li.space_index)
    @assert iseven(li.set.degree)
    q = apply_matrix(li.set.p, li.A', new_vars, div(li.set.degree, 2))
    return Sets.PolarConvexPolynomialSublevelSetAtOrigin(li.set.degree, q, nothing)
end

"""
    apply_map(li::LinearImage{Sets.InteriorDualQuadCone})

The set ``(AS)^\\circ``, the polar of the set ``AS``, is ``A^{-\\top}S^\\circ``
and given ..., we have
...
"""
function apply_map(model, li::LinearImage{<:Union{Sets.InteriorDualQuadCone,
                                                  Sets.DualConvexPolynomialCone}})
    d = data(model)
    old_vars = Sets.space_variables(li.set)
    new_vars = space_polyvars(d.spaces, li.space_index)
    q = subs(li.set.p, old_vars => li.A' * new_vars)
    return Sets.DualPolynomialSet(2, q, li.A * li.set.h, d.perspective_polyvar,
                                  new_vars)
end

# FIXME, for Sets.AbstractSet, we should apply it directly
function Base.:(*)(A::AbstractMatrix, set::Union{VariableRef, Sets.AbstractSet})
    return LinearImage(set, A, nothing)
end
