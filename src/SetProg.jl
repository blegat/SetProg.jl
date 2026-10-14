module SetProg

using LinearAlgebra

include("Sets/Sets.jl")

import MutableArithmetics
const MA = MutableArithmetics

import Reexport
Reexport.@reexport using JuMP
const MOI = JuMP.MOI
Reexport.@reexport using Polyhedra

function sos_extension()
    ext = Base.get_extension(@__MODULE__, :SetProgSumOfSquaresExt)
    if ext === nothing
        throw(ArgumentError(
            "Polynomial sets and operations require SumOfSquares. " *
            "Load it with `using SumOfSquares`.",
        ))
    end
    return ext
end

new_polyvars(n::Int) = sos_extension().new_polyvars(n)
new_perspective_variable() = sos_extension().new_perspective_variable()

export Polytope, Ellipsoid, PolySet
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

function _moment_matrix end
function apply_matrix end
function apply_transformation end
function constrain_convex end
function integrate_monomials end
function perspective_dual_polyset end
function poly_eval end
function rectangle_integrate end
function s_procedure end
function transformation end

end # module
