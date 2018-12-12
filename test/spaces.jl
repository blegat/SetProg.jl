using LinearAlgebra, Test
using DynamicPolynomials
using SetProg
const MOI = SetProg.JuMP.MOI

@testset "Spaces" begin
    B = Symmetric([1.0 0.0; 0.0 1.0])
    b = [0.0, 0.0]
    β = 1.0
    h = [0.0, 0.0]
    @polyvar x y z
    ◯ = SetProg.Sets.InteriorDualQuadCone(B, b, β, [z, x, y], h)
    mock = MOI.Utilities.MockOptimizer(JuMP.JuMPMOIModel{Float64}())
    model = JuMP.direct_model(mock);
    Q = [1.0, 0.0, 1.0]
    t = 1.0
    # The difference will be zero hence we put `zeros(6)`
    mock_optimize!(mock) = MOI.Utilities.mock_optimize!(mock, [Q; β; b; zeros(6); t])
    MOI.Utilities.set_mock_optimize!(mock, mock_optimize!)
    @variable(model, ◯◯, Ellipsoid(point=SetProg.InteriorPoint(h)))
    cref = @constraint(model, ◯◯ ⊆ ◯)
    @objective(model, Max, nth_root(volume(◯◯)))

    SetProg.optimize!(model)
    @test JuMP.termination_status(model) == MOI.Success
    @test JuMP.objective_sense(model) == MOI.MaxSense
    @test JuMP.objective_value(model) == t
    @test JuMP.value(◯◯).p == ◯.p
end
