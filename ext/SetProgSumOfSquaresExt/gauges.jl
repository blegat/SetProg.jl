function unit_set(family::SetProg.PolySet)
    d = div(family.degree, 2)
    monos = MP.monomials(family.variables, d)
    coefficients = [factorial(d) / prod(factorial, MP.exponents(mono)) for mono in monos]
    gram = GramMatrix(Matrix(Diagonal(coefficients)), monos)
    return Sets.PolySet(family.degree, gram)
end
function SetProg.gauge_value(set::Union{Sets.PolySet,Sets.ConvexPolySet}, x::AbstractVector)
    return MP.polynomial(set.p)(Sets.space_variables(set) => collect(x))
end
