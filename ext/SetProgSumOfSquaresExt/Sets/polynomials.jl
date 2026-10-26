function space_variables(set::PolySet)
    return MP.variables(set.p)
end
function ConvexPolySet(
    degree::Sets.PolynomialDegrees,
    p::GramMatrix{T,B,U},
    convexity_proof::SumOfSquares.SymMatrix{T},
) where {T,B,U}
    return ConvexPolySet{T,B,U,typeof(p),typeof(convexity_proof)}(
        degree,
        p,
        convexity_proof,
    )
end
function ConvexPolySet(
    degree::Sets.PolynomialDegrees,
    p::GramMatrix{S},
    convexity_proof::SumOfSquares.SymMatrix{T},
) where {S,T}
    U = promote_type(S, T)
    _convert(mat) = SumOfSquares.SymMatrix(convert(Vector{U}, mat.Q), mat.n)
    return ConvexPolySet(
        degree,
        GramMatrix(_convert(p.Q), p.basis),
        _convert(convexity_proof),
    )
end
function space_variables(set::ConvexPolySet)
    return MP.variables(set.p)
end
function dimension(set::ConvexPolySet)
    return MP.nvariables(set.p)
end
function zero_eliminate(set::ConvexPolySet, I)
    vars = space_variables(set)[I]
    basis_monos = MB.keys_as_monomials(set.p.basis)
    K = findall(mono -> all(var -> iszero(MP.degree(mono, var)), vars), basis_monos)
    Q = SumOfSquares.square_getindex(set.p.Q, K)
    monos = basis_monos[K]
    J = setdiff(1:dimension(set), I)
    monos =
        DynamicPolynomials.MonomialVector(monos.vars[J], Vector{Int}[z[J] for z in monos.Z])
    p = SumOfSquares.GramMatrix(Q, MB.SubBasis{MB.Monomial}(monos))
    return ConvexPolySet(set.degree, p, nothing)
end

function scaling_function(set::Union{PolySet,ConvexPolySet})
    @assert dimension(set) == 2
    # We convert the GramMatrix to a polynomial to avoid having to do the
    # conversion for every substitution.
    p = MP.polynomial(set.p)
    vars = MP.variables(p)
    @assert length(vars) == 2
    vx, vy = vars
    return (x, y) -> p(vx => x, vy => y)^(1 / set.degree)
end

function _print_gauge_function(set::ConvexPolySet; digits = 6)
    print(" ")
    p = MP.polynomial(set.p)
    if digits !== nothing
        p = round(p, digits = digits)
    end
    println(p)
end
function gauge1(set::ConvexPolynomialSet{T}) where {T}
    return MP.subs(set.q, perspective_variable(set) => one(T))
end
