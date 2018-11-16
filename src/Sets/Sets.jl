module Sets
using LinearAlgebra
using RecipesBase

abstract type AbstractSet{T} end

include("ellipsoids.jl")
include("polynomials.jl")

end
