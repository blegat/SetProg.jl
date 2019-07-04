using LinearAlgebra
using Test
using RecipesBase
using DynamicPolynomials
using SetProg, SetProg.Sets

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
            for circle in [Sets.HyperSphere(2),
                           Sets.EllipsoidAtOrigin(Symmetric(Q))]
                recipe_test(circle,
                            [1.0, 0.0, -1.0, 0.0], [0.0, 1.0, 0.0, -1.0])
                recipe_test(Sets.polar(circle),
                            [1.0, 0.0, -1.0, 0.0], [0.0, 1.0, 0.0, -1.0])
            end
        end
        @testset "Shifted Circle" begin
            for circle in [Sets.HyperSphere(2),
                           Sets.EllipsoidAtOrigin(Symmetric(Q))]
                shifted = Sets.Translation(circle, [1.0, 2.0])
                recipe_test(shifted, [2.0, 1.0, 0.0, 1.0], [2.0, 3.0, 2.0, 1.0])
            end
        end
        @testset "Scaled circle" begin
            scaled_circle = Sets.EllipsoidAtOrigin(Symmetric(2Q))
            recipe_test(scaled_circle,
                        [1/√2, 0.0, -1/√2, 0.0], [0.0, 1/√2, 0.0, -1/√2])
            recipe_test(Sets.polar(scaled_circle),
                        [√2, 0.0, -√2, 0.0], [0.0, √2, 0.0, -√2])
        end
    end
    @testset "Polynomial" begin
        @polyvar z x y
        @testset "Circle" begin
            p = SetProg.GramMatrix{Float64}((i, j) -> convert(Float64, i == j),
                                               monovec([x, y]))
            circle = Sets.PolynomialSublevelSetAtOrigin(2, p)
            recipe_test(circle,
                        [1.0, 0.0, -1.0, 0.0], [0.0, 1.0, 0.0, -1.0])
        end
        @testset "Convex Circle" begin
            p = SetProg.GramMatrix{Float64}((i, j) -> convert(Float64, i == j),
                                               monovec([x, y]))
            circle = Sets.ConvexPolynomialSublevelSetAtOrigin(2, p, nothing)
            recipe_test(circle,
                        [1.0, 0.0, -1.0, 0.0], [0.0, 1.0, 0.0, -1.0])
            recipe_test(Sets.polar(circle),
                        [-1.0, 1.0, 1.0, -1.0, -1.0], [1.0, 1.0, -1.0, -1.0, 1.0])
        end
        @testset "Scaled circle" begin
            p = SetProg.GramMatrix{Float64}((i, j) -> 2convert(Float64, i == j),
                                               monovec([x, y]))
            circle = Sets.ConvexPolynomialSublevelSetAtOrigin(2, p, nothing)
            recipe_test(circle,
                        [1/√2, 0.0, -1/√2, 0.0], [0.0, 1/√2, 0.0, -1/√2])
            recipe_test(Sets.polar(circle),
                        [-√2, √2,  √2, -√2, -√2],
                        [ √2, √2, -√2, -√2,  √2])
            hr = recipe(Sets.polar(circle))[1]
        end
        @testset "Non-homogeneous Circle" begin
            @testset "Basic" begin
                q = SetProg.GramMatrix(Float64[0 0 0
                                               0 1 0
                                               0 0 1], monovec([z, x, y]))
                shifted_circle = SetProg.perspective_dual_polyset(2, q, SetProg.InteriorPoint(zeros(2)), z, [x, y])
                recipe_test(shifted_circle,
                            [-1.0, 1.0, 1.0, -1.0, -1.0],
                            [1.0, 1.0, -1.0, -1.0, 1.0])
            end
            @testset "Scaled" begin
                q = SetProg.GramMatrix(Float64[0 0 0
                                               0 2 0
                                               0 0 2], monovec([z, x, y]))
                shifted_circle = SetProg.perspective_dual_polyset(2, q, SetProg.InteriorPoint(zeros(2)), z, [x, y])
                recipe_test(shifted_circle,
                           [-√2, √2,  √2, -√2, -√2],
                           [ √2, √2, -√2, -√2,  √2])
            end
            @testset "z-Scaled" begin
                # z: -1/2 + 1 = 1/2
                q = SetProg.GramMatrix([1/2 0 0
                                        0   1 0
                                        0   0 1], monovec([z, x, y]))
                shifted_circle = SetProg.perspective_dual_polyset(2, q, SetProg.InteriorPoint(zeros(2)), z, [x, y])
				recipe_test(shifted_circle,
							[-√2, √2,  √2, -√2, -√2],
							[ √2, √2, -√2, -√2,  √2])
            end
        end
    end
end
