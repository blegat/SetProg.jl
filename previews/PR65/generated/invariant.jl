using Polyhedra
h = HalfSpace([1, 0], 1.0) ∩ HalfSpace([-1, 0], 1) ∩ HalfSpace([0, 1], 1) ∩ HalfSpace([0, -1], 1)
□ = polyhedron(h)

using SetProg
import CSDP
sdp_solver = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true)

A = [0.0 -1.0
     1.0  0.0]

model = Model(sdp_solver)
@variable(model, S, Ellipsoid(symmetric=true))
@constraint(model, S ⊆ □)
@constraint(model, A * S ⊆ S)
@objective(model, Max, nth_root(volume(S)))
optimize!(model)
@show solve_time(model)
@show termination_status(model)
@show objective_value(model)
maximal = value(S)

model = Model(sdp_solver)
@variable(model, S, Ellipsoid(symmetric=true))
@constraint(model, □ ⊆ S)
@constraint(model, A * S ⊆ S)
@objective(model, Min, nth_root(volume(S)))
optimize!(model)
@show solve_time(model)
@show termination_status(model)
@show objective_value(model)
minimal = value(S)

using Plots
plot(ratio=:equal)
plot!(minimal)
plot!(□)
plot!(maximal)

model = Model(sdp_solver)
@variable(model, S, PolySet(degree=4, convex=true, symmetric=true))
@constraint(model, S ⊆ □)
@constraint(model, A * S ⊆ S)
@objective(model, Max, nth_root(volume(S)))
optimize!(model)
@show solve_time(model)
@show termination_status(model)
@show objective_value(model)
maximal_convex = value(S)

model = Model(sdp_solver)
@variable(model, S, PolySet(degree=4, convex=true, symmetric=true))
@constraint(model, □ ⊆ S)
@constraint(model, A * S ⊆ S)
@objective(model, Min, nth_root(volume(S)))
optimize!(model)
@show solve_time(model)
@show termination_status(model)
@show objective_value(model)
minimal_convex = value(S)

plot(ratio=:equal)
plot!(minimal_convex)
plot!(□)
plot!(maximal_convex)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

