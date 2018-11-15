module SetProg

using LinearAlgebra

include("Sets/Sets.jl")

import Reexport
Reexport.@reexport using JuMP
const MOI = JuMP.MOI
Reexport.@reexport using Polyhedra

export Ellipsoid
export nth_root

@enum Space Undecided PrimalSpace DualSpace

include("variables.jl")
include("macros.jl") # need to be before constraint as it uses ⊆ in @constraint
include("constraints.jl")
include("objective.jl")
include("data.jl")
include("optimize.jl")

end # module
