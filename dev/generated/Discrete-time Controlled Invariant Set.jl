using SetProg
import CSDP
solver = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true)
A = [1.0 0.5]
E = [1.0 0.0]
using Polyhedra
lib = DefaultLibrary{Float64}(solver)
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
end

lichen = hexcolor(0x85c17e)
canard = hexcolor(0x048b9a)
aurore = hexcolor(0xffcb60)

plot(ratio=:equal, tickfont=Plots.font(12))
plot!(□, color=lichen)
plot!(mci, color=aurore)

polar_mci = polyhedron(convexhull(
    [1.0, 0.0], [-1.0, 0.0],
        [0.0, 1.0], [0.0, -1.0],
        [1.0, 0.5], [-1.0, -0.5]
), lib)
plot(ratio=:equal, tickfont=Plots.font(12))
plot!(polar_mci, color=aurore)
plot!(◇, color=lichen)

function maximal_invariant(template, heuristic::Function)
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
function _print_gauge_function(ell::SetProg.Sets.Ellipsoid, x)
    print(" ")
    println(x' * round.(ell.Q, digits=6) * x)
end
function print_support_function(set::SetProg.Sets.Polar)
    SetProg.@polyvar x[1:2]
    print("h(S, x) =")
    _print_gauge_function(SetProg.Sets.polar(set), x)
end

sol_ell_vol = maximal_invariant(Ellipsoid(symmetric=true, dimension=2), nth_root)
print_support_function(sol_ell_vol)

primal_plot(sol_ell_vol)

polar_plot(sol_ell_vol)

sol_ell_L1 = maximal_invariant(Ellipsoid(symmetric=true, dimension=2), vol -> L1_heuristic(vol, ones(2)))
print_support_function(sol_ell_L1)

primal_plot(sol_ell_L1)

polar_plot(sol_ell_L1)

# Piecewise semi-ellipsoidal template

#We now study the maximal piecewise semi-ellipsoidal control invariant sets of a given conic partition.
#The volume is not directly maximized. Instead, for each cone, we compute the sum $s$ of the normalized rays and consider the polytope obtained by intersecting the cone with the halfspace $s^\top x \le \|s\|_2^2$. We integrate the quadratic form corresponding to this cone, i.e. $h^2(S, x)$ over the polytope. The sum of the integrals over each polytope is the objective function we use. This can be seen as the generalization of the sum of the squares of the semi-axes of the polar of the ellipsoid.

#Note that the constraint (29) of Program 1 of [LRJ20] is implemented with Proposition 2 of [LRJ20] for all results of this capsule.

function _print_gauge_function(set, x)
    println()
    for (set, piece) in zip(set.sets, set.pieces)
        print("         ")
        _print_gauge_function(set, x)
        print("              if ")
        for (i, h) in enumerate(halfspaces(piece))
            if i > 1
                print(", ")
            end
            a = -h.a
            if count(!iszero, a) == 1
                a /= abs(sum(a)) # Simplify printing
            end
            print(a'x)
            print(" ≥ 0")
        end
        println()
    end
end

sol_piece_◇ = maximal_invariant(Ellipsoid(symmetric=true, piecewise=◇), L1_heuristic)
print_support_function(sol_piece_◇)

primal_plot(sol_piece_◇, xlim=(-1.05, 1.05), ylim=(-1.05, 1.05))

polar_plot(sol_piece_◇, xlim=(-1.05, 1.05), ylim=(-1.05, 1.05))

sol_piece_mci = maximal_invariant(Ellipsoid(symmetric=true, piecewise=polar_mci), L1_heuristic)
print_support_function(sol_piece_mci)

primal_plot(sol_piece_mci)

polar_plot(sol_piece_mci)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

