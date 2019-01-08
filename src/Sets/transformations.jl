# A^{-1} * S
struct LinearPreImage{S, T, MT <: AbstractMatrix{T}}
    set::S
    A::MT
end

# S + c
struct Translation{S, T, VT <: AbstractVector{T}}
    set::S
    c::VT
end
