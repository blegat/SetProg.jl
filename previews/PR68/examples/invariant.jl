using Test     #src
# # Discrete-time Invariant Set
#
#md # [![Binder](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/invariant.ipynb)
#md # [![nbviewer](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/invariant.ipynb)
#
# In this example we compute the maximal (resp. minimal) invariant set contained
# in (resp. containing) the square with vertices $(\pm 1, \pm 1)$ for the system
# ```math
# \begin{aligned}
# x_{k+1} & = -y_k\\
# y_{k+1} & =  x_k.
# \end{aligned}
# ```
#
# The system is $x_{k+1} = A x_k$ where
# ```math
# A = \begin{bmatrix}
# 0 & -1\\
# 1 & 0
# \end{bmatrix}.
# ```
# A set $S$ is invariant if $A S \subseteq S$.

using Polyhedra
h = HalfSpace([1, 0], 1.0) ∩ HalfSpace([-1, 0], 1) ∩ HalfSpace([0, 1], 1) ∩ HalfSpace([0, -1], 1)
□ = polyhedron(h)

# We need to pick an SDP solver, see
# [here](https://jump.dev/JuMP.jl/stable/installation/#Supported-solvers)
# for a list of available ones.

using SetProg
import CSDP
sdp_solver = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true)

A = [0.0 -1.0
     1.0  0.0]

# ## Ellipsoidal template
#
# We first look for the largest ellipsoid contained in the square that is
# invariant under `A`.

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

# We then look for the smallest invariant ellipsoid containing the square.

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

# We plot the two ellipsoids together with the square.

using Plots
plot(ratio=:equal)
plot!(minimal)
plot!(□)
plot!(maximal)

# ## Quartic template
#
# We can also look for convex sublevel sets of a quartic polynomial.

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

# And the minimal containing the square:

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

# We plot the two quartic sublevel sets together with the square.

plot(ratio=:equal)
plot!(minimal_convex)
plot!(□)
plot!(maximal_convex)
