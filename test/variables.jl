using Test
using SetProg

@testset "Variables" begin
    @testset "Ellipsoid" begin
        err = ErrorException("Dimension of Ellipsoid not specified, use Ellipsoid(dimension=...)")
        @test_throws err Ellipsoid()
    end
    @testset "PolySet" begin
        err = ErrorException("Degree of PolySet not specified, use PolySet(degree=..., ...)")
        @test_throws err PolySet(dimension=1)
        @test_throws ArgumentError("Degree of PolySet not even") PolySet(degree=1)
        err = ErrorException("Dimension of Ellipsoid not specified, use Ellipsoid(dimension=..., ...)")
        @test_throws err PolySet(degree=2)
        @testset "Convex" begin
            model = Model()
            @variable(model, E, PolySet(degree=2, dimension=2))
            @objective(model, Max, nth_root(volume(E)))
            err = ErrorException("Cannot optimize volume of non-convex polynomial sublevel set. Use PolySet(convex=true, ...)")
            @test_throws err begin
                JuMP.optimize!(model)
            end
        end
    end
end
