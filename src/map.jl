struct LinearImage{S, T, MT <: AbstractMatrix{T}}
    set::S
    A::MT
end

"""
    apply_map(lm::LinearImage{Sets.PolarEllipsoidAtOrigin})

The set ``(AS)^\\circ``, the polar of the set ``AS``, is ``A^{-\\top}S^\\circ``
and given ``S^\\circ = \\{\\, x \\mid x^\\top Q x \\le 1\\,\\}``, we have
``A^{-\\top}S^\\circ = \\{\\, x \\mid x^\\top AQA^\\top x \\le 1\\,\\}``.
"""
function apply_map(model, lm::LinearImage{<:Sets.PolarEllipsoidAtOrigin})
    return Sets.PolarEllipsoidAtOrigin(Symmetric(lm.A * lm.set.Q * lm.A'))
end

"""
    apply_map(lm::LinearImage{Sets.PolarConvexPolynomialSublevelSetAtOrigin})

The set ``(AS)^\\circ``, the polar of the set ``AS``, is ``A^{-\\top}S^\\circ``
and given ``S^\\circ = \\{\\, x \\mid p(x) \\le 1\\,\\}``, we have
``A^{-\\top}S^\\circ = \\{\\, x \\mid x^\\top p(A^\\top x) \\le 1\\,\\}``.
"""
function apply_map(model,
                   lm::LinearImage{<:Sets.PolarConvexPolynomialSublevelSetAtOrigin})
    new_vars = data(model).polyvars[1:size(lm.A, 1)]
    @assert iseven(lm.set.degree)
    q = apply_matrix(lm.set.p, lm.A', new_vars, div(lm.set.degree, 2))
    return Sets.PolarConvexPolynomialSublevelSetAtOrigin(lm.set.degree, q, nothing)
end

Base.:(*)(A::AbstractMatrix, set::VariableRef) = LinearImage(set, A)
