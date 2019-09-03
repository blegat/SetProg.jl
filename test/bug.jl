using SetProg
using Polyhedra

using JuMP
const MOI = JuMP.MOI

□ = polyhedron(HalfSpace([1, 0], 1.0) ∩ HalfSpace([-1, 0], 1) ∩ HalfSpace([0, 1], 1) ∩ HalfSpace([0, -1], 1))

mock = MOI.Utilities.MockOptimizer(MOI.Utilities.Model{Float64}())
model = JuMP.direct_model(mock);
# Q = [1 0
#      0 1]
# t = √det(Q) = 1                                                                  Q11  Q12  Q22  t
MOI.Utilities.set_mock_optimize!(mock, mock -> MOI.Utilities.mock_optimize!(mock, [1.0, 0.0, 1.0, 1.0]));

@variable(model, ◯, PolySet(degree=4, dimension=2, convex=true))
cref = @constraint(model, ◯ ⊆ □)
@objective(model, Max, nth_root(volume(◯)))

SetProg.optimize!(model)
