function poly_eval(p::MP.AbstractPolynomial{JuMP.AffExpr}, a::AbstractVector{Float64})
    vars = MP.variables(p)
    aff = zero(JuMP.AffExpr)
    for term in MP.terms(p)
        mono = MP.monomial(term)
        aff = MA.add_mul!!(aff, mono(vars => a), MP.coefficient(term))
    end
    return aff
end
function sublevel_eval(set::Union{Sets.PolySet,Sets.ConvexPolySet}, a::AbstractVector)
    return poly_eval(MP.polynomial(set.p), a)
end

function poly_eval(p::MP.AbstractPolynomial{<:Real}, a::AbstractVector{<:Real})
    return p(MP.variables(p) => a)
end
