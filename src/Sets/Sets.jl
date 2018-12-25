module Sets
using LinearAlgebra
using RecipesBase
using MultivariatePolynomials
using DynamicPolynomials
using Polyhedra

const SpaceVariable = DynamicPolynomials.PolyVar{true}

abstract type AbstractSet{T} end

"""
    struct Polar{T, S<:AbstractSet{T}} <: AbstractSet{T}
        set::S
    end

The polar of the set `set`.
"""
struct Polar{T, S<:AbstractSet{T}} <: AbstractSet{T}
    set::S
end

const PolarOf{S} = Polar{<:Any, S}

const PolarOrNot{S} = Union{S, PolarOf{S}}

"""
    polar(set::AbstractSet)

Return the polar of `set`.
"""
polar(set::AbstractSet) = Polar(set)
polar(set::Polar) = set.set

"""
    polar_representation(set::AbstractSet)

Return a representation of the same set but in the polar space of the current
representation of `set`.
"""
function polar_representation end

"""
    gauge1(set::AbstractSet)

Function `f(x)` such that `set` is the 1-sublevel set of `f`, i.e.
`{ x | f(x) ≤ 1 }`.
"""
function gauge1 end

"""
    struct PerspectiveDual{S}
        set::S
    end

Set determined by the dual of the perspective cone of `set`.
"""
struct PerspectiveDual{T, S <: AbstractSet{T}} <: AbstractSet{T}
    set::S
end

"""
    perspective_dual(set::AbstractSet)

Return the set determined by the dual of the perspective cone of `set`.
"""
perspective_dual(set::AbstractSet) = PerspectiveDual(set)
perspective_dual(set::PerspectiveDual) = set.set

# TODO rename space_dimension to avoid confusion with Polyhedra.dimension
"""
    dimension(set::AbstractSet)

Return the dimension of the space where the set `set` is defined.
"""
function dimension end

function dimension(set::Union{Polar, PerspectiveDual})
    return dimension(set.set)
end

"""
    space_variables(set::AbstractSet)

Return the variables of the space where the set `set` is defined or `nothing`
if none are used.
"""
function space_variables end

function space_variables(set::Union{Polar, PerspectiveDual})
    return space_variables(set.set)
end

function perspective_variable(set::Union{Polar, PerspectiveDual})
    return perspective_variable(set.set)
end

convexity_proof(set::Union{Polar, PerspectiveDual}) = convexity_proof(set.set)

include("ellipsoids.jl")
include("polynomials.jl")

end
