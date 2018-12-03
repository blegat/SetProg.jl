module SetProg

using LinearAlgebra

include("Sets/Sets.jl")

import Reexport
Reexport.@reexport using JuMP
const MOI = JuMP.MOI
Reexport.@reexport using Polyhedra

using SumOfSquares
using MultivariateMoments
using DynamicPolynomials
const SpaceVariable = DynamicPolynomials.PolyVar{true}

export Ellipsoid, PolySet
export nth_root, L1_heuristic

@enum Space Undecided PrimalSpace DualSpace

include("utilities.jl")
include("spaces.jl")
include("variables.jl")
include("macros.jl") # need to be before constraint as it uses ⊆ in @constraint
include("map.jl")
include("constraints.jl")
include("objective.jl")
include("data.jl")
include("optimize.jl")

end # module
