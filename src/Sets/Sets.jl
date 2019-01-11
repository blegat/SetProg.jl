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
    struct PerspectiveDual{T, S <: AbstractSet{T}} <: AbstractSet{T}
        set::S
    end

Set determined by the dual of the perspective cone of `set`.
"""
struct PerspectiveDual{T, S <: AbstractSet{T}} <: AbstractSet{T}
    set::S
end

const PerspectiveDualOf{S} = PerspectiveDual{<:Any, S}

const PerspectiveDualOrPolarOrNot{S} = Union{PerspectiveDualOf{S}, PolarOrNot{S}}

"""
    perspective_dual(set::AbstractSet)

Return the set determined by the dual of the perspective cone of `set`.
"""
perspective_dual(set::AbstractSet) = PerspectiveDual(set)
perspective_dual(set::PerspectiveDual) = set.set

function scaling_function(set::PerspectiveDual)
    @assert length(space_variables(set)) == 2
    vars = [perspective_variable(set); space_variables(set)]
    # z is a halfspace of the primal so a ray of the dual
    z = [1.0, 0.0, 0.0]
    in_set(Δ::Vector) = perspective_gauge0(set.set)(vars => z + Δ) < 0
    @assert in_set(zeros(3))
    return (Δz, Δx, Δy) -> begin
        Δ = [Δz, Δx, Δy]
        _in_set(λ::Real) = in_set(Δ * λ)
        λ = 1.0
        while _in_set(λ)
            if λ > 1e10
                error("Error in plotting : the `InteriorPoint` seems to be on the boundary")
            end
            λ *= 2
        end
        λmin = 0.0
        λmax = λ
        # Binary search. Invariant: in_set(λmin) and !in_set(λmax)
        while abs(λmin - λmax) > 1e-8
            λ = (λmin + λmax) / 2
            if _in_set(λ)
                λmin = λ
            else
                λmax = λ
            end
        end
        λ = (λmin + λmax) / 2
        return 1 / λ
    end
end

# TODO rename space_dimension to avoid confusion with Polyhedra.dimension
"""
    dimension(set::AbstractSet)

Return the dimension of the space where the set `set` is defined.
"""
dimension(set::AbstractSet) = length(space_variables(set))

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
function perspective_gauge0 end
function perspective_gauge1 end

convexity_proof(set::Union{Polar, PerspectiveDual}) = convexity_proof(set.set)

struct UnknownSet{T} <: AbstractSet{T} end
include("transformations.jl")

include("ellipsoids.jl")
include("polynomials.jl")
include("recipe.jl")

end
