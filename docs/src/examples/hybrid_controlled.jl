using Test     #src

A = [0.0 1.0 0.0
     0.0 0.0 1.0
     0.0 0.0 0.0]
B = reshape([0.0, 0.0, 1.0], 3, 1)
E = [1.0 0.0 0.0
     0.0 1.0 0.0]
C = A[1:2, :]
U = [-1.0  0.0  1/4
      0.0  1.0 -1/4]

using SetProg
function maximal_invariant(family, γ = nothing; dirs=dirs)
    model = Model(sdp_solver)
    @variable(model, S, family)
    @constraint(model, S ⊆ □_3)
    @variable(model, T, family)
    @constraint(model, T ⊆ □_3)
    x = boundary_point(S, :x)
    @constraint(model, C * x in E * tangent_cone(S, x))
    @constraint(model, E * S ⊆ E * T)
    @constraint(model, U * T ⊆ E * S)
    S_2 = project(S, 1:2)
    if γ === nothing
        @variable(model, γ)
    end
    @constraint(model, [point in dirs], γ * point in S_2)
    @objective(model, Max, γ)
    JuMP.optimize!(model)
    @show solve_time(model)
    @show JuMP.termination_status(model)
    @show JuMP.objective_value(model)
    if JuMP.termination_status(model) == MOI.OPTIMAL
        return JuMP.value(S), JuMP.objective_value(model)
    else
        return
    end
end

import GLPK
lp_solver = optimizer_with_attributes(GLPK.Optimizer, MOI.Silent() => true, "presolve" => GLPK.GLP_ON)
import CSDP
sdp_solver = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true)
using Polyhedra
interval = HalfSpace([1.0], 1.0) ∩ HalfSpace([-1.0], 1.0)
lib = Polyhedra.DefaultLibrary{Float64}(lp_solver)
□_2 = polyhedron(interval * interval, lib)
□_3 = □_2 * interval
dirs = [[-1 + √3, -1 + √3], [-1 + √3, 1 - √3]]
all_dirs = [dirs; (-).(dirs)]
inner = polyhedron(vrep(all_dirs), lib)
outer = polar(inner)

# ## Ellipsoidal template
#
# We start with the ellipsoidal template. We consider two different objectives (see Section 4.2 of [L20]):
# * the volume of the set (which corresponds to $\log(\det(Q))$ or $\sqrt[n]{\det(Q)}$ in the objective function) and
# * the sum of the squares of the length of its semi-axes of the polar (which corresponds to the trace of $Q$ in the objective function).

sol_ell, γ_ell = maximal_invariant(Ellipsoid(symmetric=true))

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
cols = [canard, frambo]

x2 = range(0, stop=1, length=20)
x1 = 1 .- x2.^2 / 2
upper = [[[-1, 1]]; [[x1[i], x2[i]] for i in eachindex(x2)]]
mci = polyhedron(vrep([upper; (-).(upper)]), lib)
polar_mci = polar(mci)

SetProg.Sets.print_support_function(project(sol_ell, 1:2))

# We can plot the primal solution as follows:

function primal_plot(set, γ=nothing; npoints=256, xlim=(-1.05, 1.05), ylim=(-1.05, 1.05), args...)
    plot(ratio=:equal, tickfont=Plots.font(12); xlim=xlim, ylim=ylim, args...)
    plot!(□_2, color=lichen)
    plot!(mci, color=aurore)
    plot!(set, color=canard, npoints=npoints)
    γ === nothing || plot!(γ * inner, color=frambo)
    plot!()
end
primal_plot(project(sol_ell, 1:2), γ_ell)

# and the dual plot as follows:

function polar_plot(set, γ; npoints=256, xlim=(-1.5, 1.5), ylim=(-1.5, 1.5), args...)
    plot(ratio=:equal, tickfont=Plots.font(12); xlim=xlim, ylim=ylim, args...)
    γ === nothing || plot!(inv(γ) * outer, color=frambo)
    plot!(polar(set), color=canard, npoints=npoints)
    plot!(polar_mci, color=aurore)
    plot!(polar(□_2), color=lichen)
end
polar_plot(project(sol_ell, 1:2), γ_ell)

# ## Polyset template

p4, γ4 = maximal_invariant(PolySet(symmetric=true, degree=4, convex=true), 0.71)
γ4

p6, γ6 = maximal_invariant(PolySet(symmetric=true, degree=6, convex=true), 0.73)
γ6

p8, γ8 = maximal_invariant(PolySet(symmetric=true, degree=6, convex=true), 0.73)
γ8
