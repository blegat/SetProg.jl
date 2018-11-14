module Sets
using LinearAlgebra

abstract type AbstractSet{T} end

include("ellipsoids.jl")
include("polynomials.jl")

end
