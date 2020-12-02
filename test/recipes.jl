using LinearAlgebra
using Test
using RecipesBase
using DynamicPolynomials
using SetProg, SetProg.Sets

RecipesBase.is_key_supported(k::Symbol) = true # Plots normally defines this
function recipe(set, npoints)
    result = RecipesBase.apply_recipe(Dict{Symbol, Any}(:npoints => npoints), set)
    return result[1].args
end
function recipe_test(set, exp_x, exp_y, npoints=4)
    x, y = recipe(set, npoints)
    @test x ≈ exp_x
    @test y ≈ exp_y
end

@testset "Recipe" begin
    Q = [1.0 0.0; 0.0 1.0]
    @testset "Ellipsoid" begin
        @testset "Circle" begin
            for circle in [Sets.HyperSphere(2),
                           Sets.Ellipsoid(Symmetric(Q))]
                recipe_test(circle,
                            [1.0, 0.0, -1.0, 0.0], [0.0, 1.0, 0.0, -1.0])
                recipe_test(Sets.polar(circle),
                            [1.0, 0.0, -1.0, 0.0], [0.0, 1.0, 0.0, -1.0])
            end
        end
        @testset "Shifted Circle" begin
            for circle in [Sets.HyperSphere(2),
                           Sets.Ellipsoid(Symmetric(Q))]
                shifted = Sets.Translation(circle, [1.0, 2.0])
                recipe_test(shifted, [2.0, 1.0, 0.0, 1.0], [2.0, 3.0, 2.0, 1.0])
            end
        end
        @testset "Scaled circle" begin
            scaled_circle = Sets.Ellipsoid(Symmetric(2Q))
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
            circle = Sets.PolySet(2, p)
            recipe_test(circle,
                        [1.0, 0.0, -1.0, 0.0], [0.0, 1.0, 0.0, -1.0])
        end
        @testset "Convex Circle" begin
            p = SetProg.GramMatrix{Float64}((i, j) -> convert(Float64, i == j),
                                               monovec([x, y]))
            circle = Sets.ConvexPolySet(2, p, nothing)
            recipe_test(circle,
                        [1.0, 0.0, -1.0, 0.0], [0.0, 1.0, 0.0, -1.0])
            recipe_test(Sets.polar(circle),
                        [-1.0, -1.0, 1.0,  1.0, -1.0],
                        [-1.0,  1.0, 1.0, -1.0, -1.0])
        end
        @testset "Scaled circle" begin
            p = SetProg.GramMatrix{Float64}((i, j) -> 2convert(Float64, i == j),
                                               monovec([x, y]))
            circle = Sets.ConvexPolySet(2, p, nothing)
            recipe_test(circle,
                        [1/√2, 0.0, -1/√2, 0.0], [0.0, 1/√2, 0.0, -1/√2])
            recipe_test(Sets.polar(circle),
                        [-√2, -√2, √2,  √2, -√2],
                        [-√2,  √2, √2, -√2, -√2])
        end
        @testset "Non-homogeneous Circle" begin
            @testset "Basic" begin
                q = SetProg.GramMatrix(Float64[0 0 0
                                               0 1 0
                                               0 0 1], monovec([z, x, y]))
                shifted_circle = SetProg.perspective_dual_polyset(2, q, SetProg.InteriorPoint(zeros(2)), z, [x, y])
                recipe_test(shifted_circle,
                            [-1.0, -1.0, 1.0,  1.0, -1.0],
                            [-1.0,  1.0, 1.0, -1.0, -1.0])
            end
            @testset "Scaled" begin
                q = SetProg.GramMatrix(Float64[0 0 0
                                               0 2 0
                                               0 0 2], monovec([z, x, y]))
                shifted_circle = SetProg.perspective_dual_polyset(2, q, SetProg.InteriorPoint(zeros(2)), z, [x, y])
                recipe_test(shifted_circle,
							[-√2, -√2, √2,  √2, -√2],
							[-√2,  √2, √2, -√2, -√2])
            end
            @testset "z-Scaled" begin
                # z: -1/2 + 1 = 1/2
                q = SetProg.GramMatrix([1/2 0 0
                                        0   1 0
                                        0   0 1], monovec([z, x, y]))
                shifted_circle = SetProg.perspective_dual_polyset(2, q, SetProg.InteriorPoint(zeros(2)), z, [x, y])
				recipe_test(shifted_circle,
							[-√2, -√2, √2,  √2, -√2],
							[-√2,  √2, √2, -√2, -√2])
            end
        end
    end
    @testset "Piecewise" begin
        polytope = polyhedron(HalfSpace([1, 0], 1) ∩ HalfSpace([0, 1], 1) ∩ HalfSpace([-1, 1], 1) ∩ HalfSpace([-1, -1], 1) ∩ HalfSpace([1, -1], 1), lib)
        Q_1 = [1.0 0.0
               0.0 0.0]
        Q_2 = [0.0 0.0
               0.0 1.0]
        Q_3 = [1.0 0.0
               0.0 1.0]
        Q_4 = [1.0 1.0
               1.0 1.0]
        Q_5 = [ 1.0 -0.5
               -0.5  1.0]
        set = Sets.Piecewise(Sets.Ellipsoid.(Symmetric.([Q_1, Q_2, Q_3, Q_4, Q_5])), polytope)
        recipe_test(
            set,
            [1.0, 0.0, -1.0, 0.0],
            [0.0, 1.0, 0.0, -1.0]
        )
        recipe_test(
            set,
            [1.0, 1.0, 0.0, -1/√2, -1.0, -0.5, -0.0,  1/√3],
            [0.0, 1.0, 1.0,  1/√2,  0.0, -0.5, -1.0, -1/√3],
            8
        )
        α = 0.4142135623730951
        β = 0.732050807568877
        recipe_test(
            Sets.polar(set),
            [-1.0, -1.0, -α, 0.0, 1.0, 1.0, β, -1.0],
            [-1.0, α, 1.0, 1.0, 0.0, -β, -1.0, -1.0],
            8)
    end
end
