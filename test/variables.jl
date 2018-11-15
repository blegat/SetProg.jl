using Test
using SetProg

@testset "Variables" begin
    @testset "Ellipsoid" begin
        err = ErrorException("Dimension of Ellipsoid not specified, use Ellipsoid(dimension=...)")
        @test_throws err Ellipsoid()
    end
end
