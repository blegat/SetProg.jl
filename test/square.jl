using LinearAlgebra
using Test

using SetProg
using Polyhedra

using JuMP
const MOI = JuMP.MOI

@testset "Square" begin
    □ = polyhedron(HalfSpace([1, 0], 1.0) ∩ HalfSpace([-1, 0], 1) ∩ HalfSpace([0, 1], 1) ∩ HalfSpace([0, -1], 1))
    @testset "Ellipsoid" begin
        @testset "John" begin
            mock = MOI.Utilities.MockOptimizer(JuMP.JuMPMOIModel{Float64}())
            model = JuMP.direct_model(mock);
            # Q = [1 0
            #      0 1]
            # t = √det(Q) = 1                                                                  Q11  Q12  Q22  t
            MOI.Utilities.set_mock_optimize!(mock, mock -> MOI.Utilities.mock_optimize!(mock, [1.0, 0.0, 1.0, 1.0]));

            @variable(model, ◯, Ellipsoid(dimension=2))
            cref = @constraint(model, ◯ ⊆ □)
            @objective(model, Max, nth_root(volume(◯)))

            SetProg.optimize!(model)
            @test JuMP.termination_status(model) == MOI.Success
            @test JuMP.objective_sense(model) == MOI.MaxSense
            @test JuMP.objective_value(model) == 1.0
            @test JuMP.value(◯) isa SetProg.Sets.PolarEllipsoidAtOrigin
            @test JuMP.value(◯).Q == Symmetric([1.0 0.0; 0.0 1.0])
        end
        @testset "Löwner" begin
            mock = MOI.Utilities.MockOptimizer(JuMP.JuMPMOIModel{Float64}())
            model = JuMP.direct_model(mock);
            # Q = [√2  0
            #       0 √2]
            # t = √det(Q) = 2                                                                  Q11  Q12  Q22  t
            MOI.Utilities.set_mock_optimize!(mock, mock -> MOI.Utilities.mock_optimize!(mock, [0.5, 0.0, 0.5, 0.5]));

            @variable(model, ◯, Ellipsoid(dimension=2))
            cref = @constraint(model, □ ⊆ ◯)
            @objective(model, Min, nth_root(volume(◯)))

            SetProg.optimize!(model)
            @test JuMP.termination_status(model) == MOI.Success
            @test JuMP.objective_sense(model) == MOI.MaxSense
            @test JuMP.objective_value(model) == 0.5
            @test JuMP.value(◯) isa SetProg.Sets.EllipsoidAtOrigin
            @test JuMP.value(◯).Q == Symmetric([0.5 0.0; 0.0 0.5])
        end
    end
    @testset "Quartic" begin
        @testset "Inner" begin
            mock = MOI.Utilities.MockOptimizer(JuMP.JuMPMOIModel{Float64}())
            model = JuMP.direct_model(mock);
            # 5 variables for the polynomials
            # The PSD matrix for the convexity is 6 x 6 so 6 * (6+1) / 2 = 21
            # entries
            # 1 variable for t
            # hence 27 variables
            MOI.Utilities.set_mock_optimize!(mock, mock -> MOI.Utilities.mock_optimize!(mock, ones(27)))

            @variable(model, ◯, PolySet(degree=4, dimension=2, convex=true))
            cref = @constraint(model, ◯ ⊆ □)
            @objective(model, Max, nth_root(volume(◯)))

            SetProg.optimize!(model)
            @test JuMP.termination_status(model) == MOI.Success
            @test JuMP.objective_sense(model) == MOI.MaxSense
            @test JuMP.objective_value(model) == 1.0
            @test JuMP.value(◯) isa SetProg.Sets.PolarPolynomialSublevelSetAtOrigin{Float64}
            @test JuMP.value(◯).degree == 4
            x, y = SetProg.data(model).polyvars
            @test JuMP.value(◯).p == x^4 + x^3*y + x^2*y^2 + x*y^3 + y^4
            @test JuMP.value(◯).convexity_proof.n == 6
            @test JuMP.value(◯).convexity_proof.Q == ones(21)
            @test JuMP.objective_value(model) == 1.0
        end
        @testset "Outer" begin
            mock = MOI.Utilities.MockOptimizer(JuMP.JuMPMOIModel{Float64}())
            model = JuMP.direct_model(mock);
            MOI.Utilities.set_mock_optimize!(mock, mock -> MOI.Utilities.mock_optimize!(mock, ones(27)))

            @variable(model, ◯, PolySet(degree=4, dimension=2, convex=true))
            cref = @constraint(model, □ ⊆ ◯)
            @objective(model, Min, nth_root(volume(◯)))

            SetProg.optimize!(model)
            @test JuMP.termination_status(model) == MOI.Success
            @test JuMP.objective_sense(model) == MOI.MaxSense
            @test JuMP.objective_value(model) == 1.0
            @test JuMP.value(◯) isa SetProg.Sets.PolynomialSublevelSetAtOrigin{Float64}
            @test JuMP.value(◯).degree == 4
            x, y = SetProg.data(model).polyvars
            @test JuMP.value(◯).p == x^4 + x^3*y + x^2*y^2 + x*y^3 + y^4
            @test JuMP.value(◯).convexity_proof.n == 6
            @test JuMP.value(◯).convexity_proof.Q == ones(21)
        end
    end
end
