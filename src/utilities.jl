function quad_form(a::AbstractVector,
                   Q::AbstractMatrix,
                   b::AbstractVector)
    (n, m) = size(Q)
    n == length(a) || throw(DimensionMismatch())
    m == length(b) || throw(DimensionMismatch())
    U = MA.promote_operation(*, eltype(a), eltype(Q), eltype(b))
    out = zero(MA.promote_operation(+, U, U))
    for j in 1:m
        for i in 1:n
            out = MA.add_mul!!(out, a[i], Q[i, j], b[j])
        end
    end
    return out
end
# Same as quad_form(a, Q, a)
function quad_form(Q::Symmetric,
                   a::AbstractVector)
    return quad_form(a, Q, a)
end
function quad_form(Q::Symmetric{JuMP.VariableRef}, a::AbstractVector{<:Real})
    return quad_form(a, Q, a)
end

function psd_constraint(model, Q::Symmetric)
    return JuMP.add_constraint(model, psd_constraint(Q))
end

function psd_constraint(Q::Symmetric)
    ext = Base.get_extension(@__MODULE__, :SetProgSumOfSquaresExt)
    if ext !== nothing
        return ext.sos_psd_constraint(Q)
    end
    return JuMP.build_constraint(error, Q, JuMP.PSDCone())
end
