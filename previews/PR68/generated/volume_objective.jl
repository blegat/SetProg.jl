using Polyhedra
h = HalfSpace([1, 0], 1.0) ∩ HalfSpace([-1, 0], 1) ∩ HalfSpace([0, 1], 1) ∩ HalfSpace([0, -1], 1)
p = polyhedron(h)

using SetProg
import CSDP
sdp_solver = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true)

model = Model(sdp_solver)
@variable(model, john, Ellipsoid(symmetric=true, dimension=2))
@constraint(model, john ⊆ p)
@objective(model, Max, nth_root(volume(john)))
optimize!(model)
@show solve_time(model)
@show termination_status(model)
@show objective_value(model)
SetProg.Sets.print_support_function(value(john))

model = Model(sdp_solver)
@variable(model, löwner, Ellipsoid(symmetric=true, dimension=2))
@constraint(model, p ⊆ löwner)
@objective(model, Min, nth_root(volume(löwner)))
optimize!(model)
@show solve_time(model)
@show termination_status(model)
@show objective_value(model)
löwner_value = value(löwner)

using Plots
plot(ratio=:equal)
plot!(löwner_value)
plot!(p)
plot!(value(john))

model = Model(sdp_solver)
@variable(model, quartic_inner, PolySet(degree=4, symmetric=true, convex=true))
@constraint(model, quartic_inner ⊆ p)
@objective(model, Max, nth_root(volume(quartic_inner)))
optimize!(model)
@show solve_time(model)
@show termination_status(model)
@show objective_value(model)
quartic_inner_value = value(quartic_inner)

model = Model(sdp_solver)
@variable(model, quartic_outer, PolySet(symmetric=true, degree=4, convex=true))
@constraint(model, p ⊆ quartic_outer)
@objective(model, Min, nth_root(volume(quartic_outer)))
optimize!(model)
@show solve_time(model)
@show termination_status(model)
@show objective_value(model)
quartic_outer_value = value(quartic_outer)

plot(ratio=:equal)
plot!(quartic_outer_value)
plot!(p)
plot!(quartic_inner_value)
plot!(value(john))

function inner_L1(d)
    model = Model(sdp_solver)
    @variable(model, S, PolySet(symmetric=true, degree=d, convex=true))
    @constraint(model, S ⊆ p)
    @objective(model, Max, L1_heuristic(volume(S), [1.0, 1.0]))
    optimize!(model)
    @show solve_time(model)
    @show termination_status(model)
    @show objective_value(model)
    return value(S)
end

S2 = inner_L1(2)
S4 = inner_L1(4)
S6 = inner_L1(6)
S8 = inner_L1(8)
S10 = inner_L1(10)

plot(ratio=:equal)
plot!(quartic_outer_value)
plot!(p)
plot!(S10)
plot!(S8)
plot!(S6)
plot!(S4)
plot!(S2)

shift = 1.1
h_shift = HalfSpace([1, 0], 1.0 + shift) ∩ HalfSpace([-1, 0], 1.0 - shift) ∩ HalfSpace([0, 1], 1) ∩ HalfSpace([0, -1], 1)
p_shift = polyhedron(h_shift)

model = Model(sdp_solver)
@variable(model, john_shift, Ellipsoid(point=SetProg.InteriorPoint([shift, 0.0])))
@constraint(model, john_shift ⊆ p_shift)
@objective(model, Max, nth_root(volume(john_shift)))
optimize!(model)
@show solve_time(model)
@show termination_status(model)
@show objective_value(model)

plot(ratio=:equal)
plot!(p_shift)
plot!(value(john_shift))

model = Model(sdp_solver)
@variable(model, quartic_shift, PolySet(convex=true, degree=2, point=SetProg.InteriorPoint([shift, 0.0])))
@constraint(model, quartic_shift ⊆ p_shift)
@objective(model, Max, L1_heuristic(volume(quartic_shift), [1.0, 1.0]))
optimize!(model)
@show solve_time(model)
@show termination_status(model)
@show objective_value(model)

plot(ratio=:equal)
plot!(p_shift)
plot!(value(quartic_shift))

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

