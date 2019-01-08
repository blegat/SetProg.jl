include("square.jl")
include("controlled_invariant.jl")

using MathOptInterfaceMosek
@testset "Mosek" begin
    optimizer = MOI.Bridges.full_bridge_optimizer(MosekOptimizer(QUIET = true),
                                                  Float64)
    config = MOIT.TestConfig(atol=1e-3, rtol=1e-3, query=false)
    @testset "Square" begin
        squaretest(optimizer, config)
    end
    @testset "Controlled Invariant in Square" begin
        citest(optimizer, config)
    end
end
