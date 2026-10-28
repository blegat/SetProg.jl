
"""
    apply_map(li::Sets.LinearImage{<:Sets.PerspectiveDualOf{<:Union{Sets.PerspectiveEllipsoid,
                                                                    Sets.PerspectiveConvexPolynomialSet}}}, new_vars)

The set ``(AS)^\\circ``, the polar of the set ``AS``, is ``A^{-\\top}S^\\circ``
and given ..., we have
...
"""
function apply_map(
    li::Sets.LinearImage{<:Sets.PerspectiveDualOf{<:Sets.Householder{T}}},
    new_vars,
) where {T}
    old_vars = Sets.space_variables(li.set)
    p = MP.subs(Sets.perspective_gauge0(li.set.set), old_vars => li.A' * new_vars)
    dual = Sets.Householder(
        Sets.UnknownSet{T}(),
        p,
        li.A * li.set.set.h,
        Sets.perspective_variable(li.set),
        new_vars,
    )
    return Sets.perspective_dual(dual)
end

function apply_matrix(p::MP.AbstractPolynomial, A::AbstractMatrix, new_vars, degree)
    return MP.subs(p, MP.variables(p) => A * new_vars)
end
