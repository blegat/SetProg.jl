# Efficient implementation a' * Q * b that avoid unnecessary type promotion as
# well as unnecessary allocation
function quad_form(
    a::AbstractVector{<:Real},
    Q::SymMatrix{<:Real},
    b::AbstractVector{<:Real},
)
    n = length(a)
    @assert n == LinearAlgebra.checksquare(Q)
    @assert n == length(b)
    out = zero(typeof(zero(eltype(a)) * zero(eltype(Q)) * zero(eltype(b))))
    k = 0
    for j = 1:n
        α = a[j]
        β = b[j]
        for i = 1:(j-1)
            k += 1
            out += a[i] * Q.Q[k] * β
            out += α * Q.Q[k] * b[i]
        end
        k += 1
        out += a[j] * Q.Q[k] * b[j]
    end
    return out
end

# We have x' Q x and we want y' Q y where y is obtained by substituting
# vars for A * new_vars in x. We want to compute the matrix M such that
# y = M z where z is the vector of monomials of degree d of new_vars
# Then we will have y' Q y = z' M' Q M z
struct GramTransformation{T,MT<:MP.AbstractMonomial,MVT<:AbstractVector{MT}}
    M::Matrix{T}
    monos::MVT
end

function apply_transformation(p::SumOfSquares.GramMatrix, t::GramTransformation)
    new_n = length(t.monos)
    new_Q = [quad_form(t.M[:, i], p.Q, t.M[:, j]) for j = 1:new_n for i = 1:j]
    return GramMatrix(SymMatrix(new_Q, new_n), t.monos)
end

function transformation(old_monos, A::AbstractMatrix, new_vars, d)
    new_monos = MP.monomials(new_vars, d)
    new_n = length(new_monos)
    M = zeros(eltype(A), length(old_monos), new_n)
    mapped_vars = A * new_vars
    # Cache the result of mapped_vars[i]^n
    powers =
        [Union{Nothing,eltype(mapped_vars)}[mapped_vars[i]] for i in eachindex(mapped_vars)]
    # Compute mapped_vars[i]^n by "Power by Squaring" and cache it in `powers`
    function _power(i, n)
        @assert n > 0
        while n > length(powers[i])
            push!(powers[i], nothing)
        end
        if powers[i][n] === nothing
            if isodd(n)
                powers[i][n] = mapped_vars[i] * _power(i, n - 1)
            else
                p_2 = _power(i, div(n, 2))
                powers[i][n] = p_2 * p_2
            end
        end
        return powers[i][n]::eltype(mapped_vars)
    end
    function _map(mono::MP.AbstractMonomial)
        exps = MP.exponents(mono)
        length(exps) == length(mapped_vars) ||
            throw(ArgumentError("A monomial have less variables than `new_vars`"))
        cur = one(eltype(mapped_vars))
        for i in eachindex(exps)
            if exps[i] > 0
                cur *= _power(i, exps[i])
            end
        end
        return cur
    end
    for i in eachindex(old_monos)
        y = _map(old_monos[i])
        j = 1
        for term in MP.terms(y)
            mono = MP.monomial(term)
            while j <= length(new_monos) && mono > new_monos[j]
                j += 1
            end
            M[i, j] = MP.coefficient(term)
        end
    end
    return GramTransformation(M, new_monos)
end

# computes p ∘ A or more precisely p(variables(p) => A * new_vars)
function apply_matrix(
    p::SumOfSquares.GramMatrix{T,<:MultivariateBases.SubBasis{MultivariateBases.Monomial}},
    A::AbstractMatrix,
    new_vars,
    d,
) where {T}
    return apply_transformation(
        p,
        transformation(MultivariateBases.keys_as_monomials(p.basis), A, new_vars, d),
    )
end

# computes A # μ or more precisely p(variables(p) => A * new_vars)



function sos_psd_constraint(Q::Symmetric)
    n = LinearAlgebra.checksquare(Q)
    q = [Q[i, j] for j = 1:n for i = 1:j]
    # For n == 0, it will create no constraint, for n == 1, it will simply
    # be a Nonnegatives constraint and for n == 2 it will be a rotated SOC.
    set = SumOfSquares.matrix_cone(MOI.PositiveSemidefiniteConeTriangle, n)
    return PolyJuMP.bridgeable(
        JuMP.build_constraint(error, q, set),
        JuMP.moi_function_type(typeof(q)),
        typeof(set),
    )
end
