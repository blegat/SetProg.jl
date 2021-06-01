using Test     #src
# # Discrete-time Controlled Invariant Set
#
#md # [![Binder](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Discrete-time Controlled Invariant Set.ipynb)
#md # [![nbviewer](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Discrete-time Controlled Invariant Set.ipynb)
#
# An similar example is available in a [codeocean capsule](https://doi.org/10.24433/CO.6396918.v1).
#
# ## Introduction
#
# This example considers the linear control system already introduced in the Example 2 of [LTJ20]; see also Example 2 of [LRJ20]:
# ```math
# \begin{aligned}
# x_{k+1} & = x_k + u_k / 2\\
# u_{k+1} & = u_k'
# \end{aligned}
# ```
# with $(x_k, u_k) \in [-1, 1]^2$.
#
# The system is $x_{k+1} = Ax_k + Bu_k$ where $B = (0, 1)$ and
# ```math
# A = \begin{bmatrix}
# 1 & 1/2\\
# 0 & 0
# \end{bmatrix}.
# ```
# As shown in Example 4 of [LTJ18], a set $S \subseteq [-1, 1]^2$ is controlled invariant for this system if
# ```math
# \begin{bmatrix}
# 1 & 1/2
# \end{bmatrix}
# S \subseteq
# \begin{bmatrix}
# 1 & 0
# \end{bmatrix}
# S
# ```
#
# [LRJ20] Legat, Benoît, Saša V. Raković, and Raphaël M. Jungers.
# *Piecewise semi-ellipsoidal control invariant sets.*
# IEEE Control Systems Letters 5.3 (2020): 755-760.
#
# [LTJ18] B. Legat, P. Tabuada and R. M. Jungers.
# *Computing controlled invariant sets for hybrid systems with applications to model-predictive control*.
# 6th IFAC Conference on Analysis and Design of Hybrid Systems ADHS (2018).

# We need to pick an LP and an SDP solver, see [here](https://jump.dev/JuMP.jl/stable/installation/#Supported-solvers) for a list of available ones.

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

# ## Polyhedral template
#
# ### Fixed point approach
#
# This section shows that the maximal control invariant set of this simple control system is polyhedral. Moreover, this polyhedron is obtained after only one fixed point iteration of the standard viability
# kernel algorithm. We implement the fixed point iteration with the following function:

function fixed_point_iteration(set::Polyhedron)
    new_set = set ∩ (A \ (E * set))
    removehredundancy!(new_set)
    return new_set
end

# We start with $[-1, 1]^2$ and obtain the polytope `mci` after one iteration.

mci = fixed_point_iteration(□)

# One additional iteration gives the same polytope, showing that this polytope is control invariant.

fixed_point_iteration(mci)

# We plot this polytope `mci` in yellow below along with the safe set $[-1, 1]^2$ in green.

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

# We also plot their respective polars. Note that the inclusion is reversed in the polar space.

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

# ### Linear programming approach
#
# As introduced in [R21], given a fixed polyhedral conic partition of the state-space,
# we can search over all polyhedra for which the gauge or support function
# is linear over each piece of the partition.
# The partition can be defined by a polytope as described in [Eq. (4.6), R21].
# In this example, because the sets are represented with the support function,
# the pieces will be used for the support function/polar set and
# the pieces of the gauge function/primal will depend on the value of
# the support function on each piece.
#
# [R21] Raković, S. V.
# *Control Minkowski–Lyapunov functions*
# Automatica, Elsevier BV, 2021, 128, 109598

# Let's start with setting the partition using the cross polytope `◇`.

sol_polytope_◇ = maximal_invariant(Polytope(symmetric=true, piecewise=◇), L1_heuristic)
SetProg.Sets.print_support_function(sol_polytope_◇)

# We can see that the polar set is unbounded so the primal set has zero volume.
# More precisely, the polar is $[-1, 1] \times [-\infty, \infty]$ and the primal
# set is $[-1, 1] \times \{0\}$.

# Let's now try setting the partition using the square `□`.

sol_polytope_□ = maximal_invariant(Polytope(symmetric=true, piecewise=□), L1_heuristic)
SetProg.Sets.print_support_function(sol_polytope_□)

# The primal plot is below:

primal_plot(sol_polytope_□, xlim=(-1.05, 1.05), ylim=(-1.05, 1.05))

# And the polar plot is below:

polar_plot(sol_polytope_□, xlim=(-1.05, 1.05), ylim=(-1.05, 1.05))

# Let's now use a partition with 8 pieces.

pieces8 = convexhull(□, 1.25 * ◇)
sol_polytope_8 = maximal_invariant(Polytope(symmetric=true, piecewise=pieces8), L1_heuristic)
SetProg.Sets.print_support_function(sol_polytope_8)

# The primal plot is below:

primal_plot(sol_polytope_8, xlim=(-1.05, 1.05), ylim=(-1.05, 1.05))

# And the polar plot is below:

polar_plot(sol_polytope_8, xlim=(-1.05, 1.05), ylim=(-1.05, 1.05))

# We can use the polar of the polytope resulting from the first fixed point iteration to generate a refined conic partition.
# With this partition, the maximal piecewise semi-ellipsoidal control invariant set matches the maximal control invariant set.

sol_polytope_mci = maximal_invariant(Polytope(symmetric=true, piecewise=polar_mci), L1_heuristic)
SetProg.Sets.print_support_function(sol_polytope_mci)

# The primal plot is below:

primal_plot(sol_polytope_mci)

# And the polar plot is below:

polar_plot(sol_polytope_mci)

# ## Ellipsoidal template
#
# We now compute the maximal ellipsoidal control invariant set. Here we can either maximize its volume (which corresponds to $\log(\det(Q))$ or $\sqrt[n]{\det(Q)}$ in the objective function) or minimize the sum of the squares of the length of its semi-axes of the polar (which corresponds to the trace of $Q$ in the objective function).

# We can see below that the control invariant set of maximal volume has support function
# $$h^2(S, x) = x_1^2 - \frac{1}{2} x_1x_2 + x_2^2.$$

sol_ell_vol = maximal_invariant(Ellipsoid(symmetric=true, dimension=2), nth_root)
SetProg.Sets.print_support_function(sol_ell_vol)

# We can see below the ellipsoid in blue along with the maximal control invariant set in yellow and the safe set in green.

primal_plot(sol_ell_vol)

# Below is are the corresponding sets in the polar space.

polar_plot(sol_ell_vol)

# We can see below that the control invariant set of minimal sum of squares of the length of the semi-axes of the polar has support function
# $$h^2(S, x) = x_1^2 - \alpha x_1x_2 + x_2^2.$$
# with $\alpha \approx 1.377456$.
# Note that the sum of the squares of the length of the semi-axes of the polar is equal to the integral of $h^2(S, x)$ over the hypercube $[-1, 1]^2$ hence the use of the `L1_heuristic(vol, ones(2))` as objective which means taking the integral over the hyperrectangle with vertex `ones(2) = [1, 1]`.

sol_ell_L1 = maximal_invariant(Ellipsoid(symmetric=true, dimension=2), vol -> L1_heuristic(vol, ones(2)))
SetProg.Sets.print_support_function(sol_ell_L1)

# The primal plot is below:

primal_plot(sol_ell_L1)

# And the polar plot is below:

polar_plot(sol_ell_L1)

# ## Piecewise semi-ellipsoidal template
#
# We now study the maximal piecewise semi-ellipsoidal control invariant sets of a given conic partition.
# The volume is not directly maximized. Instead, for each cone, we compute the sum $s$ of the normalized rays and consider the polytope obtained by intersecting the cone with the halfspace $s^\top x \le \|s\|_2^2$. We integrate the quadratic form corresponding to this cone, i.e. $h^2(S, x)$ over the polytope. The sum of the integrals over each polytope is the objective function we use. This can be seen as the generalization of the sum of the squares of the semi-axes of the polar of the ellipsoid.
#
# Note that the constraint (29) of Program 1 of [LRJ20] is implemented with Proposition 2 of [LRJ20] for all results of this example.

# The conic partition are obtained by considering the conic hull of each facets of a given polytope.
# We first consider the conic partition corresponding to the polar of the safe set $[-1, 1]^2$. This gives the four quadrants as cones of the conic partition.
# The maximal piecewise semi-ellipsoidal control invariant set with this partition has the following support function:
# $$
# h^2(S, x) = \begin{cases}
# (x_1 - x_2)^2 & \text{ if }x_1x_2 \le 0,\\
# x_1^2 - x_1x_2/2 + x_2^2 & \text{ if }x_1x_2 \ge 0.
# \end{cases}
# $$
# For the cones $x_1x_2 \le 0$, the semi-ellipsoid matches the safe set and the maximal control invariant set.
# For the cones $x_1x_2 \ge 0$, the semi-ellipsoid matches the maximal volume control invariant ellipsoid.
# This illustrates one key feature of piecewise semi-ellipsoidal sets, they can combines advantages of both polyhedra and ellipsoids. It can be polyhedral on the directions where the maximal control invariant set is polyhedral and be ellipsoidal on the directions where the maximal control invariant set is smooth or requires many halfspaces in its representation.

sol_piece_◇ = maximal_invariant(Ellipsoid(symmetric=true, piecewise=◇), L1_heuristic)
SetProg.Sets.print_support_function(sol_piece_◇)

# The primal plot is below:

primal_plot(sol_piece_◇, xlim=(-1.05, 1.05), ylim=(-1.05, 1.05))

# And the polar plot is below:

polar_plot(sol_piece_◇, xlim=(-1.05, 1.05), ylim=(-1.05, 1.05))

# Let's now try setting the partition with the square `□`.

sol_piece_□ = maximal_invariant(Ellipsoid(symmetric=true, piecewise=□), L1_heuristic)
SetProg.Sets.print_support_function(sol_piece_□)

# The primal plot is below:

primal_plot(sol_piece_□, xlim=(-1.05, 1.05), ylim=(-1.05, 1.05))

# And the polar plot is below:

polar_plot(sol_piece_□, xlim=(-1.6, 1.6), ylim=(-1.6, 1.6))

# Let's now try with the partition of 8 pieces.

sol_piece_8 = maximal_invariant(Ellipsoid(symmetric=true, piecewise=pieces8), L1_heuristic)
SetProg.Sets.print_support_function(sol_piece_8)

# The primal plot is below:

primal_plot(sol_piece_8, xlim=(-1.05, 1.05), ylim=(-1.05, 1.05))

# And the polar plot is below:

polar_plot(sol_piece_8, xlim=(-1.1, 1.1), ylim=(-1.05, 1.05))

# We can use the polar of the polytope resulting from the first fixed point iteration to generate a refined conic partition.
# With this partition, the maximal piecewise semi-ellipsoidal control invariant set matches the maximal control invariant set. The support function is
# $$
# h^2(S, x) = \begin{cases}
# (x_1 - x_2)^2 & \text{ if }x_1x_2 \le 0,\\
# x_1^2 & \text{ if }x_2(x_1-2x_2) \ge 0,\\
# (x_1/2 + x_2)^2 & \text{ if }x_1(2x_2-x_1) \ge 0.\\
# \end{cases}
# $$

sol_piece_mci = maximal_invariant(Ellipsoid(symmetric=true, piecewise=polar_mci), L1_heuristic)
SetProg.Sets.print_support_function(sol_piece_mci)

# The primal plot is below:

primal_plot(sol_piece_mci)

# And the polar plot is below:

polar_plot(sol_piece_mci)

# ## Polyset template
#
# We now use the polyset templates. As details in [LRJ20], for inclusion
# constraints of the form `A * S ⊆ E * S`, the polar representation is used
# so we need the set to be convex. For this, we set the argument `convex=true`.
# We directly try `degree=4` because `degree=2` would simply give the same as
# `sol_ell_L1`.

sol4 = maximal_invariant(PolySet(symmetric=true, degree=4, convex=true), vol -> L1_heuristic(vol, ones(2)))
SetProg.Sets.print_support_function(sol4)

# The primal plot is below:

primal_plot(sol4)

# And the polar plot is below:

polar_plot(sol4)

# Let's now try degree 6.

sol6 = maximal_invariant(PolySet(symmetric=true, degree=6, convex=true), vol -> L1_heuristic(vol, ones(2)))
SetProg.Sets.print_support_function(sol6)

# The primal plot is below:

primal_plot(sol6)

# And the polar plot is below:

polar_plot(sol6)
