module SetProg

include("Sets/Sets.jl")

import Reexport
Reexport.@reexport using JuMP
Reexport.@reexport using Polyhedra

export Ellipsoid
export nth_root

include("variables.jl")
include("constraints.jl")
include("objective.jl")
include("macros.jl")
include("data.jl")

end # module
