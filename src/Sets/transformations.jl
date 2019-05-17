# A^{-1} * S
struct LinearImage{S, T, MT <: AbstractMatrix{T}} <: AbstractSet{T}
    set::S
    A::MT
end
perspective_variable(li::LinearImage) = perspective_variable(li.set)
space_variables(li::LinearImage) = nothing
dimension(li::LinearImage) = size(li.A, 1)

# A^{-1} * S
struct LinearPreImage{S, T, MT <: AbstractMatrix{T}} <: AbstractSet{T}
    set::S
    A::MT
end

# S + c
struct Translation{S, T, VT <: AbstractVector{T}} <: AbstractSet{T}
    set::S
    c::VT
end
dimension(t::Translation) = length(t.c)
space_variables(t::Translation) = space_variables(t.set)
function Polyhedra.project(t::Translation, I)
    return Translation(Polyhedra.project(t.set, I), t.c[I])
end

"""
    householder(x)

Householder reflection
```math
I - 2 v v^T / (v^T v)
```
It is symmetric and orthogonal.
"""
function householder(x)
    y = copy(x)
    t = LinearAlgebra.reflector!(y)
    v = [1; y[2:end]]
    I - t * v * v'
end
_householder(h) = householder([1.0; h]) # We add 1, for perspective variable z

struct Householder{T, S <: AbstractSet{T}, U} <: AbstractSet{T}
    set::S
    p::DynamicPolynomials.Polynomial{true, U}
    h::Vector{Float64}
    z::SpaceVariable
    x::Vector{SpaceVariable}
end
perspective_gauge0(set::Householder) = set.p
perspective_variable(set::Householder) = set.z
space_variables(set::Householder) = set.x
convexity_proof(set::Householder) = convexity_proof(set.set)

const HouseDualOf{S, T, U} = PerspectiveDualOf{Householder{T, S, U}}

function Polyhedra.project(set::Polar{T}, I) where T
    return polar(zero_eliminate(polar(set), setdiff(1:dimension(set), I)))
end

function zero_eliminate(set::Householder{T}, I) where T
    J = setdiff(1:dimension(set), I)
    p = subs(set.p, set.x[I] => zeros(T, length(I)))
    return Householder(UnknownSet{T}(), p, set.h[J], set.z, set.x[J])
end
function Polyhedra.project(set::PerspectiveDual, I)
    return perspective_dual(zero_eliminate(perspective_dual(set),
                                           setdiff(1:dimension(set), I)))
end

function _HPH(set::Householder)
    H = _householder(set.h)
    return H * _HPH(set.set) * H
end
