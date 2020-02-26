include("solver_preamble.jl")
using MosekTools
optimizer_constructor = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true)
config = MOI.Test.TestConfig(atol=1e-3, rtol=1e-3)
@testset "Mosek" begin
    @testset "Square" begin
        Tests.square_test(optimizer_constructor, config)
    end
    @testset "Invariant in Square" begin
        Tests.invariant_test(optimizer_constructor, config)
    end
    @testset "Controlled Invariant in Square" begin
        Tests.ci_test(optimizer_constructor, config)
    end
    @testset "Switched System" begin
        Tests.switched_system_test(optimizer_constructor, config)
    end
end
