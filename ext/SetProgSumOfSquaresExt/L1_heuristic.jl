function rectangle_integrate(m::MP.AbstractMonomialLike, vertex)
    exp = MP.exponents(m)
    if any(isodd, exp)
        return zero(power_integrate(2, vertex[1]))
    else
        @assert length(exp) == length(vertex)
        return prod(power_integrate.(exp, vertex))
    end
end
function rectangle_integrate(t::MP.AbstractTermLike, vertex)
    return MP.coefficient(t) * rectangle_integrate(MP.monomial(t), vertex)
end
function rectangle_integrate(p::MP.AbstractPolynomialLike, vertex)
    return sum(rectangle_integrate(t, vertex) for t in MP.terms(p))
end
function integrate_monomials(
    monos::AbstractVector{<:MP.AbstractMonomial},
    polytope::Polyhedron,
)
    return integrate_decompositions(decompose.(monos), polytope)
end
function integrate(p::MP.AbstractPolynomial, polytope::Polyhedron)
    return MP.coefficients(p)'integrate_monomials(MP.monomials(p), polytope)
end
decompose(mono::MP.AbstractMonomial) = decompose(MP.exponents(mono))

function all_exponents(set::Union{Sets.PolySet,Sets.ConvexPolySet})
    return MP.exponents.(MP.monomials(Sets.space_variables(set), set.degree))
end
function evaluate_monomials(
    monomial_value::Function,
    set::Union{Sets.PolySet{T},Sets.ConvexPolySet{T}},
) where {T}
    U = MA.promote_operation(*, Float64, T)
    total = zero(MA.promote_operation(+, U, U))
    for t in MP.terms(set.p)
        total = MA.add_mul!!(
            total,
            monomial_value(MP.exponents(MP.monomial(t))),
            MP.coefficient(t),
        )
    end
    return total
end
function l1_integral(set::Sets.HouseDualOf{<:Sets.ConvexPolynomialSet}, vertex)
    return rectangle_integrate(MP.subs(set.set.set.q, set.set.set.z => 0), 1 ./ vertex)
end
