function poly_eval(p::AbstractPolynomial{JuMP.AffExpr},
                   a::AbstractVector{Float64})
    vars = variables(p)
    aff = zero(JuMP.AffExpr)
    for term in terms(p)
        mono = monomial(term)
        aff = JuMP.destructive_add!(aff, mono(vars => a), coefficient(term))
    end
    return aff
end
function sublevel_eval(ell::Sets.EllipsoidAtOrigin,
                       a::AbstractVector)
    return quad_form(ell.Q, a)
end
function sublevel_eval(set::Union{Sets.PolynomialSublevelSetAtOrigin,
                                  Sets.ConvexPolynomialSublevelSetAtOrigin},
                       a::AbstractVector)
    return poly_eval(polynomial(set.p), a)
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
