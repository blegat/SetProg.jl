using Polyhedra
using SetProg

using JuMP
const MOI = JuMP.MOI

@testset "Square" begin
    h = HalfSpace([1, 0], 1.0) ∩ HalfSpace([-1, 0], 1) ∩ HalfSpace([0, 1], 1) ∩ HalfSpace([0, -1], 1)
    p = polyhedron(h)
    mock = MOI.Utilities.MockOptimizer(JuMP.JuMPMOIModel{Float64}())
    model = JuMP.direct_model(mock);
    # Q = [1 0
    #      0 1]
    # t = det(Q)^(1/n) = 1
    #                                                                                  Q11  Q12  Q22  t
    MOI.Utilities.set_mock_optimize!(mock, mock -> MOI.Utilities.mock_optimize!(mock, [1.0, 0.0, 1.0, 1.0]));

    @variable(model, E, Ellipsoid(2))
    cref = @constraint(model, E ⊆ p)
    @objective(model, Max, nth_root(volume(E)))

    SetProg.optimize!(model)
    @test JuMP.termination_status(model) == MOI.Success
    @test JuMP.objective_value(model) == 1.0
    @test JuMP.value(E).Q == Symmetric([1.0 0.0; 0.0 1.0])
end
