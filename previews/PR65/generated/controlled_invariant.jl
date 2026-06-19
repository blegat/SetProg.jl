using Polyhedra
h = HalfSpace([1, 0], 1.0) ∩ HalfSpace([-1, 0], 1) ∩ HalfSpace([0, 1], 1) ∩ HalfSpace([0, -1], 1)
□ = polyhedron(h)

using SetProg
import CSDP
sdp_solver = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true)

Δt = 0.5
A = [1.0 Δt]
E = [1.0 0.0]

model = Model(sdp_solver)
@variable(model, S, Ellipsoid(symmetric=true))
@constraint(model, S ⊆ □)
@constraint(model, A * S ⊆ E * S)
@objective(model, Max, nth_root(volume(S)))
optimize!(model)
@show solve_time(model)
@show termination_status(model)
@show objective_value(model)
ell = value(S)

function Sd(d)
    model = Model(sdp_solver)
    @variable(model, S, PolySet(degree=d, convex=true, symmetric=true))
    @constraint(model, S ⊆ □)
    @constraint(model, A * S ⊆ E * S)
    @objective(model, Max, L1_heuristic(volume(S), [1.0, 1.0]))
    optimize!(model)
    @show solve_time(model)
    @show termination_status(model)
    @show objective_value(model)
    return value(S)
end

S4 = Sd(4)
S8 = Sd(8)
S12 = Sd(12)
S16 = Sd(16)
S20 = Sd(20)
S22 = Sd(22)

mci = □ ∩ HalfSpace([1.0, 0.5], 1.0) ∩ HalfSpace([-1.0, -0.5], 1.0)
polar_□ = polar(□)
polar_mci = polar(mci)

using Plots
plot(ratio=:equal)
plot!(□)
plot!(mci)
plot!(S22)
plot!(S16)
plot!(S8)
plot!(S4)
plot!(ell)

plot(ratio=:equal)
plot!(SetProg.Sets.polar(S4); npoints=512)
plot!(SetProg.Sets.polar(S8); npoints=512)
plot!(SetProg.Sets.polar(S16); npoints=512)
plot!(SetProg.Sets.polar(S22); npoints=512)
plot!(polar_mci)
plot!(polar_□)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

