module Sets
using LinearAlgebra
using RecipesBase

abstract type AbstractSet{T} end

"""
    polar(set::AbstractSet)

Return the polar of `set`.
"""
function polar end

include("ellipsoids.jl")
include("polynomials.jl")

end
