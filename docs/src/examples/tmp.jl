function maximal_invariant(template, heuristic = v -> L1_heuristic(v, [1.0, 1.0]))
    model = Model(sdp_solver)
    @variable(model, S, template)
    @constraint(model, S ⊆ □)
    @constraint(model, A * S ⊆ E * S)
    @objective(model, Max, heuristic(volume))
    JuMP.optimize!(model)
    @show solve_time(model)
    @show JuMP.termination_status(model)
    @show JuMP.objective_value(model)
    return JuMP.value(S)
end
Sd(d) = maximal_invariant(PolySet(degree=d, convex=true, symmetric=true))

# ## Ellipsoidal template
#
# We start with the ellipsoidal template. We consider two different objectives (see Section 4.2 of [L20]):
# * the volume of the set (which corresponds to $\log(\det(Q))$ or $\sqrt[n]{\det(Q)}$ in the objective function) and
# * the sum of the squares of the length of its semi-axes of the polar (which corresponds to the trace of $Q$ in the objective function).

sol_ell = maximal_invariant(Ellipsoid(symmetric=true), nth_root)

using Plots
function hexcolor(rgb::UInt32)
    r = ((0xff0000 & rgb) >> 16) / 255
    g = ((0x00ff00 & rgb) >>  8) / 255
    b = ((0x0000ff & rgb)      ) / 255
    Plots.RGBA(r, g, b)
end
# Values taken from http://www.toutes-les-couleurs.com/code-couleur-rvb.php
lichen = hexcolor(0x85c17e)
canard = hexcolor(0x048b9a)
aurore = hexcolor(0xffcb60)
frambo = hexcolor(0xc72c48)
cols = [canard, frambo]

x2 = range(0, stop=1, length=20)
x1 = 1 .- x2.^2 / 2
upper = [[[-1, 1]]; [[x1[i], x2[i]] for i in eachindex(x2)]]
mci = polyhedron(vrep([upper; (-).(upper)]), lib)
polar_mci = polar(mci)

function _print_gauge_function(ell::SetProg.Sets.Ellipsoid, x)
    print(" ")
    println(x' * round.(ell.Q, digits=3) * x)
end
function print_support_function(set::SetProg.Sets.Polar)
    SetProg.@polyvar x[1:SetProg.Sets.dimension(set)]
    print("h(S, x) =")
    _print_gauge_function(polar(set), x)
end
print_support_function(project(sol_ell, 1:2))

# We can plot the primal solution as follows:

function primal_plot(set; npoints=256, xlim=(-1.05, 1.05), ylim=(-1.05, 1.05), args...)
    plot(ratio=:equal, tickfont=Plots.font(12); xlim=xlim, ylim=ylim, args...)
    plot!(□, color=lichen)
    plot!(mci, color=aurore)
    plot!(set, color=canard, npoints=npoints)
end
primal_plot(project(sol_ell, 1:2))

# and the dual plot as follows:

function polar_plot(set; npoints=256, xlim=(-1.5, 1.5), ylim=(-1.5, 1.5), args...)
    plot(ratio=:equal, tickfont=Plots.font(12); xlim=xlim, ylim=ylim, args...)
    plot!(polar(set), color=canard, npoints=npoints)
    plot!(polar_mci, color=aurore)
    plot!(polar(□), color=lichen)
end
polar_plot(project(sol_ell, 1:2))

# ## Polyset template

p4, γ4 = maximal_invariant(PolySet(symmetric=true, degree=4, convex=true), 0.91)
γ4

# Below is the primal plot:

primal_plot(project(p4, 1:2), γ4)

# and here is the polar plot:

polar_plot(project(p4, 1:2), γ4)

# ## Piecewise semi-ellipsoidal template

sol_piece_◇, γ_piece_◇ = maximal_invariant(Ellipsoid(symmetric=true, piecewise=polar(□_3)), dirs=all_dirs)
γ_piece_◇

# Below is the primal plot:

primal_plot(project(sol_piece_◇, 1:2), γ_piece_◇)

# and here is the polar plot:

polar_plot(project(sol_piece_◇, 1:2), γ_piece_◇)
