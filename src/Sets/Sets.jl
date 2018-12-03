module Sets
using LinearAlgebra
using RecipesBase
using MultivariatePolynomials
using DynamicPolynomials
using Polyhedra

const SpaceVariable = DynamicPolynomials.PolyVar{true}

abstract type AbstractSet{T} end

# TODO rename space_dimension to avoid confusion with Polyhedra.dimension
"""
    dimension(set::AbstractSet)

Return the dimension of the space where the set `set` is defined.
"""
function dimension end

"""
    space_variables(set::AbstractSet)

Return the variables of the space where the set `set` is defined or `nothing`
if none are used.
"""
function space_variables end

"""
    polar(set::AbstractSet)

Return the polar of `set`.
"""
function polar end

include("ellipsoids.jl")
include("polynomials.jl")

end
