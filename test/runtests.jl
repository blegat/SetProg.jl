include("extensions.jl")
include("sos_extension.jl")

import GLPK
const lp_solver = GLPK.Optimizer
import Polyhedra
const lib = Polyhedra.DefaultLibrary{Float64}(lp_solver)

include("Sets.jl")
include("intersection.jl")
include("apply.jl")
include("variables.jl")
include("L1_heuristic.jl")
include("recipes.jl")
include("model.jl")
include("separation.jl")
include("hypatia.jl")
include("mock_tests.jl")
