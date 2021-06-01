using SetProg
import Clp
lp_solver = optimizer_with_attributes(Clp.Optimizer, MOI.Silent() => true)
import CSDP
sdp_solver = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true)
using Polyhedra
lib = DefaultLibrary{Float64}(lp_solver)
h = HalfSpace([1, 0], 1.0) ∩ HalfSpace([-1, 0], 1) ∩ HalfSpace([0, 1], 1) ∩ HalfSpace([0, -1], 1)
□ = polyhedron(h, lib)
◇ = polar(□)
A1 = [-1 -1
      -4  0] / 4
A2 = [ 3  3
      -2  1] / 4
A = [0.0 -1.0
     1.0  0.0]

function maximal_invariant(template, heuristic::Function)
    solver = if template isa Polytope
        lp_solver
    else
        sdp_solver
    end
    model = Model(solver)
    @variable(model, S, template)
    @constraint(model, S ⊆ □)
    @constraint(model, A1 * S ⊆ S)
    @constraint(model, A2 * S ⊆ S)
    @objective(model, Max, heuristic(volume(S)))
    optimize!(model)
    @show solve_time(model)
    @show termination_status(model)
    @show objective_value(model)
    return value(S)
end
function minimal_invariant(template, heuristic::Function)
    solver = if template isa Polytope
        lp_solver
    else
        sdp_solver
    end
    model = Model(solver)
    @variable(model, S, template)
    @constraint(model, □ ⊆ S)
    @constraint(model, A1 * S ⊆ S)
    @constraint(model, A2 * S ⊆ S)
    @objective(model, Min, heuristic(volume(S)))
    optimize!(model)
    @show solve_time(model)
    @show termination_status(model)
    @show objective_value(model)
    @show JuMP.objective_sense(model)
    return value(S)
end

using Plots
function hexcolor(rgb::UInt32)
    r = ((0xff0000 & rgb) >> 16) / 255
    g = ((0x00ff00 & rgb) >>  8) / 255
    b = ((0x0000ff & rgb)      ) / 255
    Plots.RGBA(r, g, b)
end # Values taken from http://www.toutes-les-couleurs.com/code-couleur-rvb.php
lichen = hexcolor(0x85c17e)
canard = hexcolor(0x048b9a)
aurore = hexcolor(0xffcb60)
frambo = hexcolor(0xc72c48)

function primal_plot(min_sol, max_sol; npoints=256, args...)
    plot(ratio=:equal, tickfont=Plots.font(12); args...)
    plot!(min_sol, color=canard, npoints=npoints)
    plot!(□, color=lichen)
    plot!(max_sol, color=aurore, npoints=npoints)
end
function polar_plot(min_sol, max_sol; npoints=256, args...)
    plot(ratio=:equal, tickfont=Plots.font(12); args...)
    plot!(polar(max_sol), color=aurore, npoints=npoints)
    plot!(◇, color=lichen)
    plot!(polar(min_sol), color=canard, npoints=npoints)
end

function backward_fixed_point_iteration(set::Polyhedron)
    new_set = set ∩ (A1 \ set) ∩ (A2 \ set)
    removehredundancy!(new_set)
    return new_set
end

max_polytope_1 = backward_fixed_point_iteration(□)

backward_fixed_point_iteration(max_polytope_1)

function forward_fixed_point_iteration(set::Polyhedron)
    new_set = convexhull(set, A1 * set, A2 * set)
    removevredundancy!(new_set)
    return new_set
end

min_polytope_1 = forward_fixed_point_iteration(□)

min_polytope_2 = forward_fixed_point_iteration(min_polytope_1)

min_polytope_3 = forward_fixed_point_iteration(min_polytope_2)

forward_fixed_point_iteration(min_polytope_3)

plot(ratio=:equal, tickfont=Plots.font(12))
plot!(min_polytope_3, color=frambo)
plot!(min_polytope_2, color=canard)
plot!(min_polytope_1, color=aurore)
plot!(□, color=lichen)
plot!(max_polytope_1, color=frambo)

pieces8 = convexhull(□, 1.25 * ◇)
max_polytope = maximal_invariant(Polytope(symmetric=true, piecewise=pieces8), L1_heuristic)
min_polytope = minimal_invariant(Polytope(symmetric=true, piecewise=polar(pieces8)), L1_heuristic)
primal_plot(min_polytope, max_polytope)

max_ell_vol = maximal_invariant(Ellipsoid(symmetric=true), nth_root)
min_ell_vol = minimal_invariant(Ellipsoid(symmetric=true), nth_root)
primal_plot(min_ell_vol, max_ell_vol)
polar_plot(min_ell_vol, max_ell_vol)

max_ell_L1 = maximal_invariant(Ellipsoid(symmetric=true), vol -> L1_heuristic(vol, ones(2)))
min_ell_L1 = minimal_invariant(Ellipsoid(symmetric=true), vol -> L1_heuristic(vol, ones(2)))
primal_plot(min_ell_L1, max_ell_L1)
savefig("AJPR14e54ell.png")

max_4 = maximal_invariant(PolySet(symmetric=true, convex=true, degree=4), vol -> L1_heuristic(vol, ones(2)))
min_4 = minimal_invariant(PolySet(symmetric=true, degree=4), vol -> L1_heuristic(vol, ones(2)))
primal_plot(min_4, max_4)

max_6 = maximal_invariant(PolySet(symmetric=true, convex=true, degree=6), vol -> L1_heuristic(vol, ones(2)))
min_6 = minimal_invariant(PolySet(symmetric=true, degree=6), vol -> L1_heuristic(vol, ones(2)))
primal_plot(min_6, max_6)

max_8 = maximal_invariant(PolySet(symmetric=true, convex=true, degree=8), vol -> L1_heuristic(vol, □))
min_8 = minimal_invariant(PolySet(symmetric=true, degree=8), vol -> L1_heuristic(vol, □))
primal_plot(min_8, max_8, npoints=1024)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

