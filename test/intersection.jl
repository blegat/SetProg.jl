module TestIntersection

using Test
using LinearAlgebra
using SetProg
import GLPK

const Sets = SetProg.Sets

@testset "Intersection metadata and values" begin
    domain = polyhedron(hrep([-1.0 0.0; 0.0 -1.0], zeros(2)))
    ellipsoid = Sets.Ellipsoid(Symmetric([2.0 0.0; 0.0 3.0]))
    set = Sets.Intersection(ellipsoid, domain)
    @test Sets.dimension(set) == 2
    @test Sets.space_variables(set) === nothing
    @test Sets.perspective_variable(set) === nothing
    @test !SetProg.needs_polyvars(set)
    fitted = value(set)
    @test fitted.set.Q == ellipsoid.Q
    @test fitted.domain === domain

    model = Model(GLPK.Optimizer)
    set_silent(model)
    @variable(model, a[i=1:2] == i)
    variable_set = Sets.Intersection(Sets.PolarPoint(a), domain)
    optimize!(model)
    fitted = value(variable_set)
    @test fitted.set.a == [1.0, 2.0]
    @test fitted.domain === domain
end

@testset "Membership respects both the set and its domain" begin
    # The offset checks that the domain bound scales with a ScaledPoint.
    domain = polyhedron(hrep([-1.0 0.0; 0.0 -1.0; 1.0 0.0], [0.0, 0.0, 0.5]))
    set = Sets.Intersection(Sets.PolarPoint([1.0, 1.0]), domain)
    for (point, feasible) in (
        ([0.25, 0.25], true),
        ([-0.25, 0.25], false),
        ([0.75, 0.0], false),
        ([0.25, 1.0], false),
        (SetProg.ScaledPoint([1.0, 0.5], 2.0), true),
        (SetProg.ScaledPoint([1.5, 0.0], 2.0), false),
    )
        model = Model(GLPK.Optimizer)
        set_silent(model)
        @constraint(model, point in set)
        optimize!(model)
        @test termination_status(model) == (feasible ? MOI.OPTIMAL : MOI.INFEASIBLE)
    end
end

@testset "A variable point in a lower-dimensional domain" begin
    domain = polyhedron(hrep([HyperPlane([1.0, -1.0], 0.5)], [HalfSpace([-1.0, 0.0], 0.0)]))
    set = Sets.Intersection(Sets.PolarPoint([1.0, 1.0]), domain)
    model = Model(GLPK.Optimizer)
    set_silent(model)
    @variable(model, x[1:2])
    @constraint(model, x in set)
    @objective(model, Max, x[2])
    optimize!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test value.(x) ≈ [0.75, 0.25]

    # Scaling must apply to the affine hyperplane offset too.
    for (coords, feasible) in (([1.5, 0.5], true), ([1.0, 1.0], false))
        model = Model(GLPK.Optimizer)
        set_silent(model)
        point = SetProg.ScaledPoint(coords, 2.0)
        @constraint(model, point in set)
        optimize!(model)
        @test termination_status(model) == (feasible ? MOI.OPTIMAL : MOI.INFEASIBLE)
    end
end

end
