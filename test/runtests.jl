import GLPK
const lp_solver = GLPK.Optimizer
import Polyhedra
const lib = Polyhedra.DefaultLibrary{Float64}(lp_solver)

include("Sets.jl")
include("apply.jl")
include("variables.jl")
include("L1_heuristic.jl")
include("recipes.jl")
include("mock_tests.jl")
