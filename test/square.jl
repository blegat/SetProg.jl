using LinearAlgebra
using Test

using Polyhedra
using SetProg

using JuMP
const MOI = JuMP.MOI

@testset "Square" begin
    □ = polyhedron(HalfSpace([1, 0], 1.0) ∩ HalfSpace([-1, 0], 1) ∩ HalfSpace([0, 1], 1) ∩ HalfSpace([0, -1], 1))
    @testset "John ellispoid" begin
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
    @testset "Löwner ellispoid" begin
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
