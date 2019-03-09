# Efficient implementation a' * Q * b that avoid unnecessary type promotion as
# well as unnecessary allocation
function quad_form(a::AbstractVector{<:Real},
                   Q::Union{Symmetric{JuMP.VariableRef},
                            SymMatrix{JuMP.VariableRef}},
                   b::AbstractVector{<:Real})
    n = length(a)
    @assert n == LinearAlgebra.checksquare(Q)
    @assert n == length(b)
    aff = zero(JuMP.GenericAffExpr{eltype(a), JuMP.VariableRef})
    for j in 1:n
        for i in 1:n
            JuMP.add_to_expression!(aff, a[i] * b[j], Q[i, j])
        end
    end
    return aff
end
function quad_form(a::AbstractVector{<:Real},
                   Q::Union{Symmetric{<:Real},
                            SymMatrix{<:Real}},
                   b::AbstractVector{<:Real})
    n = length(a)
    @assert n == LinearAlgebra.checksquare(Q)
    @assert n == length(b)
    return sum(a[i] * Q[i, j] * b[j] for j in 1:n for i in 1:n)
end

# Same as quad_form(a, Q, a)
function quad_form(Q::Symmetric{<:JuMP.AbstractJuMPScalar},
                   a::AbstractVector{<:AbstractMonomialLike})
    n = length(a)
    @assert n == LinearAlgebra.checksquare(Q)
    return sum((i == j ? 1 : 2) * a[i] * Q[i, j] * a[j] for j in 1:n for i in 1:j)
end
function quad_form(Q::Symmetric{JuMP.VariableRef}, a::AbstractVector{<:Real})
    return quad_form(a, Q, a)
end


# computes p ∘ A or more precisely p(variables(p) => A * new_vars)
function apply_matrix(p::SumOfSquares.GramMatrix,
                      A::AbstractMatrix, new_vars, d)
    vars = variables(p)
    # We have x' Q x and we want y' Q y where y is obtained by substituting
    # vars for A * new_vars in x. We want to compute the matrix M such that
    # y = M z where z is the vector of monomials of degree d of new_vars
    # Then we will have y' Q y = z' M' Q M z
    z = monomials(new_vars, d)
    new_n = length(z)
    M = zeros(eltype(A), length(p.x), new_n)
    for i in 1:length(p.x)
        y = p.x[i](vars => A * new_vars)
        j = 1
        for term in terms(y)
            mono = monomial(term)
            while j <= length(z) && mono < z[j]
                j += 1
            end
            M[i, j] = coefficient(term)
        end
    end
    new_Q = [quad_form(M[:, i], p.Q, M[:, j]) for j in 1:new_n for i in 1:j]
    return GramMatrix(SymMatrix(new_Q, new_n), z)
end
