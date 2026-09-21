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
import CSDP
sdp_solver = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true)

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
