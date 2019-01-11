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
            circle = Sets.EllipsoidAtOrigin(Symmetric(Q))
            recipe_test(circle,
                        [1.0, 0.0, -1.0, 0.0], [0.0, 1.0, 0.0, -1.0])
            recipe_test(Sets.polar(circle),
                        [1.0, 0.0, -1.0, 0.0], [0.0, 1.0, 0.0, -1.0])
        end
        @testset "Shifted Circle" begin
            circle = Sets.Translation(Sets.EllipsoidAtOrigin(Symmetric(Q)),
                                      [1.0, 2.0])
            recipe_test(circle,
                        [2.0, 1.0, 0.0, 1.0], [2.0, 3.0, 2.0, 1.0])
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
            p = SetProg.MatPolynomial{Float64}((i, j) -> convert(Float64, i == j),
                                               monovec([x, y]))
            circle = Sets.ConvexPolynomialSublevelSetAtOrigin(2, p, nothing)
            recipe_test(circle,
                        [1.0, 0.0, -1.0, 0.0], [0.0, 1.0, 0.0, -1.0])
            hr = recipe(Sets.polar(circle))[1]
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
            circle = Sets.ConvexPolynomialSublevelSetAtOrigin(2, p, nothing)
            recipe_test(circle,
                        [1/√2, 0.0, -1/√2, 0.0], [0.0, 1/√2, 0.0, -1/√2])
            hr = recipe(Sets.polar(circle))[1]
            @test !hashyperplanes(hr)
            hss = collect(halfspaces(hr))
            @test hss[1].a ≈ [1/√2, 0.0]
            @test hss[1].β ≈ 1.0
            @test hss[2].a ≈ [0.0, 1/√2]
            @test hss[2].β ≈ 1.0
            @test hss[3].a ≈ [-1/√2, 0.0]
            @test hss[3].β ≈ 1.0
            @test hss[4].a ≈ [0.0, -1/√2]
            @test hss[4].β ≈ 1.0
        end
        @testset "Non-homogeneous Circle" begin
            @testset "Basic" begin
                q = SetProg.MatPolynomial(Float64[0 0 0
                                                  0 1 0
                                                  0 0 1], monovec([z, x, y]))
                shifted_circle = SetProg.perspective_dual_polyset(2, q, SetProg.InteriorPoint(zeros(2)), z, [x, y])
                hr = recipe(shifted_circle)[1]
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
            @testset "Scaled" begin
                q = SetProg.MatPolynomial(Float64[0 0 0
                                                  0 2 0
                                                  0 0 2], monovec([z, x, y]))
                shifted_circle = SetProg.perspective_dual_polyset(2, q, SetProg.InteriorPoint(zeros(2)), z, [x, y])
                hr = recipe(shifted_circle)[1]
                @test !hashyperplanes(hr)
                hss = collect(halfspaces(hr))
                @test hss[1].a ≈ [1/√2, 0.0]
                @test hss[1].β ≈ 1.0
                @test hss[2].a ≈ [0.0, 1/√2]
                @test hss[2].β ≈ 1.0
                @test hss[3].a ≈ [-1/√2, 0.0]
                @test hss[3].β ≈ 1.0
                @test hss[4].a ≈ [0.0, -1/√2]
                @test hss[4].β ≈ 1.0
            end
            @testset "z-Scaled" begin
                # z: -1/2 + 1 = 1/2
                q = SetProg.MatPolynomial([1/2 0 0
                                           0   1 0
                                           0   0 1], monovec([z, x, y]))
                shifted_circle = SetProg.perspective_dual_polyset(2, q, SetProg.InteriorPoint(zeros(2)), z, [x, y])
                hr = recipe(shifted_circle)[1]
                @test !hashyperplanes(hr)
                hss = collect(halfspaces(hr))
                @test hss[1].a ≈ [1/√2, 0.0]
                @test hss[1].β ≈ 1.0
                @test hss[2].a ≈ [0.0, 1/√2]
                @test hss[2].β ≈ 1.0
                @test hss[3].a ≈ [-1/√2, 0.0]
                @test hss[3].β ≈ 1.0
                @test hss[4].a ≈ [0.0, -1/√2]
                @test hss[4].β ≈ 1.0
            end
        end
    end
end
