using Test     #src
# # Volume Objective
#
#md # [![Binder](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/volume_objective.ipynb)
#md # [![nbviewer](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/volume_objective.ipynb)
#
# In this example, we compute the maximal (resp. minimal) volume ellipsoids
# and polynomial sublevel sets contained in (resp. containing) the square with
# vertices $(\pm 1, \pm 1)$. We start by defining the square with
# [Polyhedra](https://github.com/JuliaPolyhedra/Polyhedra.jl).

using Polyhedra
h = HalfSpace([1, 0], 1.0) ∩ HalfSpace([-1, 0], 1) ∩ HalfSpace([0, 1], 1) ∩ HalfSpace([0, -1], 1)
p = polyhedron(h)

# We need to pick an SDP solver, see
# [here](https://jump.dev/JuMP.jl/stable/installation/#Supported-solvers)
# for a list of available ones.

using SetProg
import Hypatia
sdp_solver = optimizer_with_attributes(Hypatia.Optimizer, MOI.Silent() => true)

# ## John ellipsoid
#
# The maximal volume ellipsoid contained in a convex body is called its John
# ellipsoid. The John ellipsoid for our square can be computed as follows.

model = Model(sdp_solver)
@variable(model, john, Ellipsoid(symmetric=true, dimension=2))
@constraint(model, john ⊆ p)
@objective(model, Max, nth_root(volume(john)))
optimize!(model)
@show solve_time(model)
@show termination_status(model)
@show objective_value(model)
SetProg.Sets.print_support_function(value(john))

# ### Which cone does the solver receive ?
#
# The objective `nth_root(volume(john))` is reformulated by SetProg into a
# new variable `t` constrained by `[t; vec(Q)] in MOI.RootDetConeTriangle(2)`
# where `Q` is the matrix of the ellipsoid. We can see these reformulated
# constraints with `list_of_constraint_types`:

list_of_constraint_types(model)

# The solver may not support every one of these constraints natively. In that
# case, JuMP transforms them using *bridges*. For instance, a solver only
# supporting positive semidefinite constraints would receive the
# `RootDetConeTriangle` reformulated with a `GeometricMeanCone` and a
# `PositiveSemidefiniteConeTriangle`, so with additional variables and
# constraints. Hypatia supports the root-determinant cone natively, and we can
# check that with `print_active_bridges`. It shows, for each constraint
# type of the model, which bridges are used and which constraints the solver
# receives in the end.

print_active_bridges(model)

# To only look at the root-determinant cone, we give its function and set type:

print_active_bridges(model, Vector{VariableRef}, MOI.RootDetConeTriangle)

# The only bridge is `SetDotScalingBridge`. It only rescales the off-diagonal
# entries into the `MOI.Scaled{MOI.RootDetConeTriangle}` set, which is the
# form in which Hypatia supports the cone natively. There is no
# `GeometricMeanCone` or `PositiveSemidefiniteConeTriangle`.

# ### Log-determinant objective
#
# Instead of the `n`th root of the determinant, we can also maximize its
# logarithm with `log(volume(john))`. The optimal ellipsoid is the same, but
# the objective value is now `log(det(Q)) = log(1) = 0` instead of
# `det(Q)^(1/2) = 1`.

model = Model(sdp_solver)
@variable(model, john_log, Ellipsoid(symmetric=true, dimension=2))
@constraint(model, john_log ⊆ p)
@objective(model, Max, log(volume(john_log)))
optimize!(model)
@show termination_status(model)
@show objective_value(model)
@test objective_value(model) ≈ 0 atol=1e-5 #src
SetProg.Sets.print_support_function(value(john_log))

# The objective is now reformulated into `[t; 1; vec(Q)] in MOI.LogDetConeTriangle(2)`
# which means `t ≤ log(det(Q))`. As the constant `1` is part of the function,
# it is a `Vector{AffExpr}` instead of a `Vector{VariableRef}`.
# Hypatia also supports this cone natively:

print_active_bridges(model, Vector{AffExpr}, MOI.LogDetConeTriangle)

# ## Löwner ellipsoid
#
# The minimal volume ellipsoid containing a convex body is called its Löwner
# ellipsoid. The Löwner ellipsoid for our square can be computed as follows.

model = Model(sdp_solver)
@variable(model, löwner, Ellipsoid(symmetric=true, dimension=2))
@constraint(model, p ⊆ löwner)
@objective(model, Min, nth_root(volume(löwner)))
optimize!(model)
@show solve_time(model)
@show termination_status(model)
@show objective_value(model)
löwner_value = value(löwner)

# We can visualize the Löwner and John ellipsoids as follows.

using Plots
plot(ratio=:equal)
plot!(löwner_value)
plot!(p)
plot!(value(john))

# ## Higher degree polynomials
#
# Ellipsoids are the sublevel sets of positive definite *quadratic* forms. To
# allow for more sophisticated shapes, we instead look for sublevel sets of
# *quartic* forms. For this, we simply replace `Ellipsoid(dimension=2)` by
# `PolySet(degree=4, dimension=2)`. Note that the quantities optimized are
# not exactly the volume anymore but provide a reasonable heuristic.
#
# ### Maximal volume quartic sublevel set contained in the square

model = Model(sdp_solver)
@variable(model, quartic_inner, PolySet(degree=4, symmetric=true, convex=true))
@constraint(model, quartic_inner ⊆ p)
@objective(model, Max, nth_root(volume(quartic_inner)))
optimize!(model)
@show solve_time(model)
@show termination_status(model)
@show objective_value(model)
quartic_inner_value = value(quartic_inner)

# ### Minimal volume quartic sublevel set containing the square

model = Model(sdp_solver)
@variable(model, quartic_outer, PolySet(symmetric=true, degree=4, convex=true))
@constraint(model, p ⊆ quartic_outer)
@objective(model, Min, nth_root(volume(quartic_outer)))
optimize!(model)
@show solve_time(model)
@show termination_status(model)
@show objective_value(model)
quartic_outer_value = value(quartic_outer)

# We can visualize the quartic sublevel sets as follows.

plot(ratio=:equal)
plot!(quartic_outer_value)
plot!(p)
plot!(quartic_inner_value)
plot!(value(john))

# ### Inner sublevel sets of increasing degree
#
# We can also explore how the inner sublevel set tightens as the degree grows,
# using the `L1_heuristic` as a tractable volume proxy.

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

# ## Non-homogeneous case
#
# For non-symmetric bodies, the John/Löwner ellipsoids are not centered at the
# origin. We need to provide a `point` as a hint of an interior point. We
# illustrate on a horizontally-shifted square.

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

# Likewise for the quartic sublevel set.

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
