module SetProg

using LinearAlgebra

include("Sets/Sets.jl")

import MutableArithmetics
const MA = MutableArithmetics

import Reexport
Reexport.@reexport using MultivariateBases
Reexport.@reexport using JuMP
const MOI = JuMP.MOI
Reexport.@reexport using Polyhedra

using SumOfSquares
using MultivariateMoments
import DynamicPolynomials
import MultivariatePolynomials
const MP = MultivariatePolynomials
const SpaceVariable = DynamicPolynomials.PolyVar{true}

export Ellipsoid, PolySet
export nth_root, L1_heuristic

@enum Space Undecided PrimalSpace DualSpace

include("utilities.jl")
include("spaces.jl")
include("macros.jl") # need to be before `variables.jl` and `constraints.jl` as they use `⊆` in `@constraint`
include("variables.jl")
include("map.jl")
include("constraints.jl")
include("objective.jl")
include("data.jl")
include("optimize.jl")

end # module
