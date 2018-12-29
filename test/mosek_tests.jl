include("square.jl")
include("controlled_invariant.jl")

using MathOptInterfaceMosek
@testset "Mosek" begin
    optimizer = MOI.Bridges.fullbridgeoptimizer(MosekOptimizer(QUIET = true),
                                                Float64)
    config = MOIT.TestConfig(atol=1e-3, rtol=1e-3, query=false)
    # TODO Need to improve SOSConvex by splitting the Newton polytope for y and x
    # in order for quartic tests to pass
    @testset "Square" begin
        squaretest(optimizer, config, ["quartic_inner_homogeneous",
                                       "quartic_outer_homogeneous"])
    end
    @testset "Controlled Invariant in Square" begin
        citest(optimizer, config)
    end
end
