using Test     #src
# # Controlled Invariant Set
#
#md # [![Binder](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/controlled_invariant.ipynb)
#md # [![nbviewer](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/controlled_invariant.ipynb)
#
# In this example, we compute the maximal (resp. minimal) controlled invariant
# set contained in the square with vertices $(\pm 1, \pm 1)$ for the system
# ```math
# \begin{aligned}
# v_{k+1} & = v_k + a_k \Delta t\\
# a_{k+1} & = u_k
# \end{aligned}
# ```
# using the technique developed in [LTJ18].
#
# The system is $x_{k+1} = A x_k + B u_k$ where $B = (0, 1)$ and
# ```math
# A = \begin{bmatrix}
# 1 & \Delta t\\
# 0 & 0
# \end{bmatrix}.
# ```
# As shown in [LTJ18], a set $S$ is controlled invariant if
# ```math
# \begin{bmatrix}
# 1 & \Delta t
# \end{bmatrix} S
# \subseteq
# \begin{bmatrix}
# 1 & 0
# \end{bmatrix} S.
# ```
#
# [LTJ18] B. Legat, P. Tabuada and R. M. Jungers.
# *Computing controlled invariant sets for hybrid systems with applications to model-predictive control*.
# 6th IFAC Conference on Analysis and Design of Hybrid Systems ADHS 2018, **2018**.

using Polyhedra
h = HalfSpace([1, 0], 1.0) ∩ HalfSpace([-1, 0], 1) ∩ HalfSpace([0, 1], 1) ∩ HalfSpace([0, -1], 1)
□ = polyhedron(h)

# We need to pick an SDP solver, see
# [here](https://jump.dev/JuMP.jl/stable/installation/#Supported-solvers)
# for a list of available ones.

using SetProg
import CSDP
sdp_solver = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true)

Δt = 0.5
A = [1.0 Δt]
E = [1.0 0.0]

# ## Ellipsoidal template
#
# We first look for the largest controlled invariant ellipsoid contained in the
# square.

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

# ## Polynomial sublevel sets
#
# To allow for tighter approximations, we now use sublevel sets of polynomials
# of degree $d$. The `L1_heuristic` objective is a tractable proxy for the
# volume; see Section 4.2 of [L20] for details.
#
# [L20] Legat, B. (2020). *Set programming : theory and computation*. Ph.D. thesis, UCLouvain.

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

# We compute the maximal controlled invariant sublevel set for several degrees.

S4 = Sd(4)
S8 = Sd(8)
S12 = Sd(12)
S16 = Sd(16)
S20 = Sd(20)
S22 = Sd(22)

# ## Plotting
#
# Below are reference polytopes used in the plots: the maximal controlled
# invariant polytope `mci` and its polar.

mci = □ ∩ HalfSpace([1.0, 0.5], 1.0) ∩ HalfSpace([-1.0, -0.5], 1.0)
polar_□ = polar(□)
polar_mci = polar(mci)

# We plot the polynomial sublevel sets together with the square and `mci`.

using Plots
plot(ratio=:equal)
plot!(□)
plot!(mci)
plot!(S22)
plot!(S16)
plot!(S8)
plot!(S4)
plot!(ell)

# And the corresponding polar plot:

plot(ratio=:equal)
plot!(SetProg.Sets.polar(S4); npoints=512)
plot!(SetProg.Sets.polar(S8); npoints=512)
plot!(SetProg.Sets.polar(S16); npoints=512)
plot!(SetProg.Sets.polar(S22); npoints=512)
plot!(polar_mci)
plot!(polar_□)
