using LinearAlgebra
using Test
using RecipesBase
using DynamicPolynomials
using SetProg

RecipesBase.is_key_supported(k::Symbol) = true # Plots normally defines this
function recipe(set)
    result = RecipesBase.apply_recipe(Dict{Symbol, Any}(:npoints => 4), set)
    return result[1].args
end
function recipe_test(set, exp_x, exp_y)
    x, y = recipe(set)
    @test x ≈ exp_x
    @test y ≈ exp_y
end

@testset "Recipe" begin
    Q = [1.0 0.0; 0.0 1.0]
    @testset "Ellipsoid" begin
        @testset "Circle" begin
            circle = SetProg.Sets.EllipsoidAtOrigin(Symmetric(Q))
            recipe_test(circle,
                        [1.0, 0.0, -1.0, 0.0], [0.0, 1.0, 0.0, -1.0])
            recipe_test(SetProg.Sets.polar(circle),
                        [1.0, 0.0, -1.0, 0.0], [0.0, 1.0, 0.0, -1.0])
        end
        @testset "Scaled circle" begin
            scaled_circle = SetProg.Sets.EllipsoidAtOrigin(Symmetric(2Q))
            recipe_test(scaled_circle,
                        [1/√2, 0.0, -1/√2, 0.0], [0.0, 1/√2, 0.0, -1/√2])
            recipe_test(SetProg.Sets.polar(scaled_circle),
                        [√2, 0.0, -√2, 0.0], [0.0, √2, 0.0, -√2])
        end
    end
    @testset "Polynomial" begin
        @polyvar x y
        @testset "Circle" begin
            p = SetProg.MatPolynomial{Float64}((i, j) -> convert(Float64, i == j),
                                               monovec([x, y]))
            circle = SetProg.Sets.PolynomialSublevelSetAtOrigin(2, p, nothing)
            recipe_test(circle,
                        [1.0, 0.0, -1.0, 0.0], [0.0, 1.0, 0.0, -1.0])
            hr = recipe(SetProg.Sets.polar(circle))[1]
            @test !hashyperplanes(hr)
            hss = collect(halfspaces(hr))
            @test hss[1].a ≈ [1.0, 0.0]
            @test hss[1].β ≈ 1.0
            @test hss[2].a ≈ [0.0, 1.0]
            @test hss[2].β ≈ 1.0
            @test hss[3].a ≈ [-1.0, 0.0]
            @test hss[3].β ≈ 1.0
            @test hss[4].a ≈ [0.0, -1.0]
            @test hss[4].β ≈ 1.0
        end
        @testset "Scaled circle" begin
            p = SetProg.MatPolynomial{Float64}((i, j) -> 2convert(Float64, i == j),
                                               monovec([x, y]))
            circle = SetProg.Sets.PolynomialSublevelSetAtOrigin(2, p, nothing)
            recipe_test(circle,
                        [1/√2, 0.0, -1/√2, 0.0], [0.0, 1/√2, 0.0, -1/√2])
            hr = recipe(SetProg.Sets.polar(circle))[1]
            @test !hashyperplanes(hr)
            hss = collect(halfspaces(hr))
            @test hss[1].a ≈ [1.0, 0.0]
            @test hss[1].β ≈ √2
            @test hss[2].a ≈ [0.0, 1.0]
            @test hss[2].β ≈ √2
            @test hss[3].a ≈ [-1.0, 0.0]
            @test hss[3].β ≈ √2
            @test hss[4].a ≈ [0.0, -1.0]
            @test hss[4].β ≈ √2
        end
    end
end
