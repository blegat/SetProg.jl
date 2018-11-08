module SetProg

include("Sets/Sets.jl")

import Reexport
Reexport.@reexport using JuMP
Reexport.@reexport using Polyhedra

export Ellipsoid
export nth_root

@enum Space Undecided PrimalSpace DualSpace

include("variables.jl")
include("constraints.jl")
include("objective.jl")
include("macros.jl")
include("data.jl")
include("optimize.jl")

end # module
