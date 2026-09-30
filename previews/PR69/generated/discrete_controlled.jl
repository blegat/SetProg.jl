using SetProg
import GLPK
lp_solver = optimizer_with_attributes(GLPK.Optimizer, MOI.Silent() => true)
import CSDP
sdp_solver = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true)
A = [1.0 0.5]
E = [1.0 0.0]
using Polyhedra
lib = DefaultLibrary{Float64}(lp_solver)
h = HalfSpace([1, 0], 1.0) ∩ HalfSpace([-1, 0], 1) ∩ HalfSpace([0, 1], 1) ∩ HalfSpace([0, -1], 1)
□ = polyhedron(h, lib) # [0, 1]^2
v = convexhull([1.0, 0], [0, 1], [-1, 0], [0, -1])
◇ = polyhedron(v, lib); # polar of [0, 1]^2

function fixed_point_iteration(set::Polyhedron)
    new_set = set ∩ (A \ (E * set))
    removehredundancy!(new_set)
    return new_set
end

mci = fixed_point_iteration(□)

fixed_point_iteration(mci)

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

plot(ratio=:equal, tickfont=Plots.font(12))
plot!(□, color=lichen)
plot!(mci, color=aurore)

polar_mci = polar(mci)
plot(ratio=:equal, tickfont=Plots.font(12))
plot!(polar_mci, color=aurore)
plot!(◇, color=lichen)

function maximal_invariant(template, heuristic::Function)
    solver = if template isa Polytope
        lp_solver
    else
        sdp_solver
    end
    model = Model(solver)
    @variable(model, S, template)
    @constraint(model, S ⊆ □)
    @constraint(model, A * S ⊆ E * S)
    @objective(model, Max, heuristic(volume(S)))
    optimize!(model)
    @show solve_time(model)
    @show termination_status(model)
    @show objective_value(model)
    return value(S)
end
function primal_plot(set; npoints=256, args...)
    plot(ratio=:equal, tickfont=Plots.font(12); args...)
    plot!(□, color=lichen)
    plot!(mci, color=aurore)
    plot!(set, color=canard, npoints=npoints)
end
function polar_plot(set; npoints=256, args...)
    plot(ratio=:equal, tickfont=Plots.font(12); args...)
    plot!(SetProg.Sets.polar(set), color=canard, npoints=npoints)
    plot!(polar_mci, color=aurore)
    plot!(◇, color=lichen)
end

sol_polytope_◇ = maximal_invariant(Polytope(symmetric=true, piecewise=◇), L1_heuristic)
SetProg.Sets.print_support_function(sol_polytope_◇)

sol_polytope_□ = maximal_invariant(Polytope(symmetric=true, piecewise=□), L1_heuristic)
SetProg.Sets.print_support_function(sol_polytope_□)

primal_plot(sol_polytope_□, xlim=(-1.05, 1.05), ylim=(-1.05, 1.05))

polar_plot(sol_polytope_□, xlim=(-1.05, 1.05), ylim=(-1.05, 1.05))

pieces8 = convexhull(□, 1.25 * ◇)
sol_polytope_8 = maximal_invariant(Polytope(symmetric=true, piecewise=pieces8), L1_heuristic)
SetProg.Sets.print_support_function(sol_polytope_8)

primal_plot(sol_polytope_8, xlim=(-1.05, 1.05), ylim=(-1.05, 1.05))

polar_plot(sol_polytope_8, xlim=(-1.05, 1.05), ylim=(-1.05, 1.05))

sol_polytope_mci = maximal_invariant(Polytope(symmetric=true, piecewise=polar_mci), L1_heuristic)
SetProg.Sets.print_support_function(sol_polytope_mci)

primal_plot(sol_polytope_mci)

polar_plot(sol_polytope_mci)

sol_ell_vol = maximal_invariant(Ellipsoid(symmetric=true, dimension=2), nth_root)
SetProg.Sets.print_support_function(sol_ell_vol)

primal_plot(sol_ell_vol)

polar_plot(sol_ell_vol)

sol_ell_L1 = maximal_invariant(Ellipsoid(symmetric=true, dimension=2), vol -> L1_heuristic(vol, ones(2)))
SetProg.Sets.print_support_function(sol_ell_L1)

primal_plot(sol_ell_L1)

polar_plot(sol_ell_L1)

sol_piece_◇ = maximal_invariant(Ellipsoid(symmetric=true, piecewise=◇), L1_heuristic)
SetProg.Sets.print_support_function(sol_piece_◇)

primal_plot(sol_piece_◇, xlim=(-1.05, 1.05), ylim=(-1.05, 1.05))

polar_plot(sol_piece_◇, xlim=(-1.05, 1.05), ylim=(-1.05, 1.05))

sol_piece_□ = maximal_invariant(Ellipsoid(symmetric=true, piecewise=□), L1_heuristic)
SetProg.Sets.print_support_function(sol_piece_□)

primal_plot(sol_piece_□, xlim=(-1.05, 1.05), ylim=(-1.05, 1.05))

polar_plot(sol_piece_□, xlim=(-1.6, 1.6), ylim=(-1.6, 1.6))

sol_piece_8 = maximal_invariant(Ellipsoid(symmetric=true, piecewise=pieces8), L1_heuristic)
SetProg.Sets.print_support_function(sol_piece_8)

primal_plot(sol_piece_8, xlim=(-1.05, 1.05), ylim=(-1.05, 1.05))

polar_plot(sol_piece_8, xlim=(-1.1, 1.1), ylim=(-1.05, 1.05))

sol_piece_mci = maximal_invariant(Ellipsoid(symmetric=true, piecewise=polar_mci), L1_heuristic)
SetProg.Sets.print_support_function(sol_piece_mci)

primal_plot(sol_piece_mci)

polar_plot(sol_piece_mci)

sol4 = maximal_invariant(PolySet(symmetric=true, degree=4, convex=true), vol -> L1_heuristic(vol, ones(2)))
SetProg.Sets.print_support_function(sol4)

primal_plot(sol4)

polar_plot(sol4)

sol6 = maximal_invariant(PolySet(symmetric=true, degree=6, convex=true), vol -> L1_heuristic(vol, ones(2)))
SetProg.Sets.print_support_function(sol6)

primal_plot(sol6)

polar_plot(sol6)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

