using LinearAlgebra
using Test

using SetProg
using Polyhedra
using MultivariatePolynomials

using JuMP
const MOI = JuMP.MOI

@testset "Square" begin
    □ = polyhedron(HalfSpace([1, 0], 1.0) ∩ HalfSpace([-1, 0], 1) ∩ HalfSpace([0, 1], 1) ∩ HalfSpace([0, -1], 1))
    function square_test(inner::Bool, variable::SetProg.AbstractVariable,
                         metric::Function, mock_optimize!::Function,
                         objective_value, set_test)
        mock = MOI.Utilities.MockOptimizer(JuMP.JuMPMOIModel{Float64}())
        model = JuMP.direct_model(mock);
        # Q = [1 0
        #      0 1]
        # t = √det(Q) = 1                                                                  Q11  Q12  Q22  t
        MOI.Utilities.set_mock_optimize!(mock, mock_optimize!)

        @variable(model, ◯, variable)
        if inner
            cref = @constraint(model, ◯ ⊆ □)
        else
            cref = @constraint(model, □ ⊆ ◯)
        end
        @objective(model, inner ? MOI.MaxSense : MOI.MinSense,
                   metric(volume(◯)))

        SetProg.optimize!(model)
        @test JuMP.termination_status(model) == MOI.Success
        @test JuMP.objective_sense(model) == MOI.MaxSense
        @test JuMP.objective_value(model) == objective_value
        set_test(JuMP.value(◯))
    end
    @testset "Ellipsoid" begin
        @testset "John" begin
            # Q = [1 0
            #      0 1]
            # t = √det(Q) = 1
            Q = [1.0, 0.0, 1.0]
            t = 1.0
            @testset "Homogeneous" begin
                square_test(true, Ellipsoid(symmetric=true, dimension=2), nth_root,
                            mock -> MOI.Utilities.mock_optimize!(mock, [Q; t]),
                            1.0,
                            ◯ -> begin
                                @test ◯ isa SetProg.Sets.PolarEllipsoidAtOrigin
                                @test ◯.Q == Symmetric([1.0 0.0; 0.0 1.0])
                            end)
            end
            @testset "Non-homogeneous" begin
                @testset "Ellipsoid" begin
                    β = 1.0
                    b = [0.0, 0.0]
                    square_test(true,
                                Ellipsoid(point=SetProg.InteriorPoint([0.0, 0.0])),
                                nth_root,
                                mock -> MOI.Utilities.mock_optimize!(mock, [Q; β; b; t]),
                                1.0,
                                ◯ -> begin
                                    @test ◯ isa SetProg.Sets.InteriorDualQuadCone{Float64,Float64}
                                    z, x, y = variables(◯.p)
                                    @test ◯.p == z^2 + x^2 + y^2
                                    @test ◯.Q == Symmetric([1.0 0.0; 0.0 1.0])
                                    @test ◯.b == [0.0, 0.0]
                                    @test ◯.β == 1.0
                                    @test ◯.H ≈ [-1.0 0.0 0.0
                                                  0.0 1.0 0.0
                                                  0.0 0.0 1.0]
                                end)
                end
                @testset "PolySet" begin
                    β = 1.0
                    b = [0.0, 0.0]
                    square_test(true,
                                PolySet(degree=2, convex=true, point=SetProg.InteriorPoint([0.0, 0.0])),
                                set -> L1_heuristic(set, [1.0, 1.0]),
                                mock -> begin
                                    # β-1 b[1] b[2]
                                    #  .  Q[1] Q[2]
                                    #  .   .   Q[3]
                                    MOI.Utilities.mock_optimize!(mock, [β-1; b[1]; Q[1]; b[2]; Q[2]; Q[3]; 2Q])
                                end,
                                8/3,
                                ◯ -> begin
                                    @test ◯ isa SetProg.Sets.DualConvexPolynomialCone{Float64,Float64}
                                    z, x, y = variables(◯.p)
                                    @test ◯.p == -z^2 + x^2 + y^2
                                end)
                end
            end
        end
        @testset "Löwner" begin
            mock = MOI.Utilities.MockOptimizer(JuMP.JuMPMOIModel{Float64}())
            model = JuMP.direct_model(mock);
            # Q = [√2  0
            #       0 √2]
            # t = √det(Q) = 2                                                                  Q11  Q12  Q22  t
            square_test(false, Ellipsoid(symmetric=true, dimension=2), nth_root,
                        mock -> MOI.Utilities.mock_optimize!(mock, [0.5, 0.0, 0.5, 0.5]),
                        0.5,
                        ◯ -> begin
                            @test ◯ isa SetProg.Sets.EllipsoidAtOrigin
                            @test ◯.Q == Symmetric([0.5 0.0; 0.0 0.5])
                        end)
        end
    end
    @testset "Quartic" begin
        @testset "Inner" begin
            # The PSD matrix for the variable  is 3 x 3 so 3 * (3+1) / 2 = 6
            # The PSD matrix for the convexity is 6 x 6 so 6 * (6+1) / 2 = 21
            # entries
            # 1 variable for t
            # hence 28 variables
            square_test(true,
                        PolySet(symmetric=true, degree=4, dimension=2, convex=true),
                        nth_root,
                        mock -> MOI.Utilities.mock_optimize!(mock, ones(28)),
                        1.0,
                        ◯ -> begin
                            @test ◯ isa SetProg.Sets.PolarConvexPolynomialSublevelSetAtOrigin{Float64}
                            @test ◯.degree == 4
                            x, y = variables(◯.p)
                            @test polynomial(◯.p) == x^4 + 2x^3*y + 3x^2*y^2 + 2x*y^3 + y^4
                            @test ◯.convexity_proof.n == 6
                            @test ◯.convexity_proof.Q == ones(21)
                        end)
        end
        @testset "Outer" begin
            square_test(false,
                        PolySet(symmetric=true, degree=4, dimension=2, convex=true),
                        nth_root,
                        mock -> MOI.Utilities.mock_optimize!(mock, ones(28)),
                        1.0,
                        ◯ -> begin
                            @test ◯ isa SetProg.Sets.ConvexPolynomialSublevelSetAtOrigin{Float64}
                            @test ◯.degree == 4
                            x, y = variables(◯.p)
                            @test polynomial(◯.p) == x^4 + 2x^3*y + 3x^2*y^2 + 2x*y^3 + y^4
                            @test ◯.convexity_proof.n == 6
                            @test ◯.convexity_proof.Q == ones(21)
                        end)
        end
    end
end
