struct LinearImage{S, T, MT <: AbstractMatrix{T}}
    set::S
    A::MT
    space::Union{Nothing, SpaceIndex}
    variables::Union{Nothing, Vector{SpaceVariable}}
    space_index::Union{Nothing, SpaceIndex}
end

function clear_space(li::LinearImage)
    li.space_index = nothing
end
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
    new_vars = li.variables
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
    new_vars = li.variables
    q = subs(li.set.p, old_vars => li.A' * new_vars)
    return Sets.DualPolynomialSet(2, q, li.A * li.set.h, d.perspective_polyvar,
                                  new_vars)
end

# FIXME, for Sets.AbstractSet, we should apply it directly
function Base.:(*)(A::AbstractMatrix, set::Union{VariableRef, Sets.AbstractSet})
    return LinearImage(set, A, nothing, nothing)
end
