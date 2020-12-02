module Tests

using JuMP
import GLPK
const lp_solver = optimizer_with_attributes(GLPK.Optimizer, "presolve" => GLPK.GLP_ON)
import Polyhedra
const lib = Polyhedra.DefaultLibrary{Float64}(lp_solver)

include("utilities.jl")

include("square.jl")
include("invariant.jl")
include("controlled_invariant.jl")
include("switched_system.jl")

end
