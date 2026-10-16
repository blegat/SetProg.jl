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

@testset "Linear preimages of intersections" begin
    # The affine upper bound distinguishes the domain from its recession cone.
    domain = polyhedron(hrep([-1.0 0.0; 0.0 -1.0; 1.0 0.0], [0.0, 0.0, 0.5]))
    for inner in
        (Sets.PolarPoint([1.0, 1.0]), Sets.Ellipsoid(Symmetric(Matrix{Float64}(I, 2, 2))))
        set = Sets.Intersection(inner, domain)
        for (A, cases) in (
            (
                [2.0 0.0; 0.0 1.0],
                (
                    ([0.125, 0.25], true),
                    ([-0.125, 0.25], false),
                    ([0.375, 0.0], false),
                    ([0.25, 1.0], false),
                ),
            ),
            (
                [1.0 0.0 0.0; 0.0 1.0 0.0],
                (
                    ([0.25, 0.25, 10.0], true),
                    ([-0.25, 0.25, 0.0], false),
                    ([0.75, 0.0, 0.0], false),
                    ([0.5, 1.0, -5.0], false),
                ),
            ),
            (
                reshape([1.0, 1.0], 2, 1),
                (([0.25], true), ([-0.25], false), ([0.75], false)),
            ),
            (
                [1.0 0.0; 0.0 0.0],
                (
                    ([0.25, 100.0], true),
                    ([0.25, -100.0], true),
                    ([-0.25, 0.0], false),
                    ([0.75, 0.0], false),
                ),
            ),
        )
            preimage = SetProg.apply_map(Sets.LinearPreImage(set, A), nothing)
            @test Sets.dimension(preimage) == size(A, 2)
            for (point, feasible) in cases
                model = Model(GLPK.Optimizer)
                set_silent(model)
                @constraint(model, point in preimage)
                optimize!(model)
                @test termination_status(model) == (feasible ? MOI.OPTIMAL : MOI.INFEASIBLE)
            end
        end
    end
end

@testset "Linear preimages retain variable coefficients" begin
    domain = polyhedron(hrep([-1.0 0.0; 0.0 -1.0], zeros(2)))
    for (point, feasible) in
        (([0.25, 0.0], true), ([0.75, 0.0], false), ([-0.25, 0.0], false))
        model = Model(GLPK.Optimizer)
        set_silent(model)
        @variable(model, c[i=1:2] == i)
        set = Sets.Intersection(Sets.PolarPoint(c), domain)
        preimage = SetProg.apply_map(Sets.LinearPreImage(set, [2.0 0.0; 0.0 1.0]), nothing)
        @constraint(model, point in preimage)
        optimize!(model)
        @test termination_status(model) == (feasible ? MOI.OPTIMAL : MOI.INFEASIBLE)
        if feasible
            @test value(preimage).set.a ≈ [2.0, 2.0]
        end
    end
end

@testset "Inclusions of copositive linear gauges" begin
    lib = Polyhedra.DefaultLibrary{Float64}(
        optimizer_with_attributes(
            GLPK.Optimizer,
            "presolve" => GLPK.GLP_ON,
            MOI.Silent() => true,
        ),
    )
    orthant = polyhedron(hrep(-Matrix{Float64}(I, 2, 2), zeros(2)), lib)
    halfplane = polyhedron(hrep(reshape([-1.0, 0.0], 1, 2), [0.0]), lib)
    diagonal =
        polyhedron(hrep([HyperPlane([1.0, -1.0], 0.0)], [HalfSpace([-1.0, 0.0], 0.0)]), lib)
    source = Sets.Intersection(Sets.PolarPoint([2.0, 3.0]), orthant)
    target = Sets.Intersection(Sets.PolarPoint([1.0, 2.0]), orthant)
    strip = Sets.Intersection(Sets.PolarPoint([1.0, 0.0]), halfplane)
    ray = Sets.Intersection(Sets.PolarPoint([2.0, 0.0]), diagonal)
    origin = Sets.Intersection(
        Sets.PolarPoint([1.0, 1.0]),
        polyhedron(
            hrep(
                [HyperPlane([1.0, 0.0], 0.0), HyperPlane([0.0, 1.0], 0.0)],
                HalfSpace{Float64,Vector{Float64}}[],
            ),
            lib,
        ),
    )
    for (subset, supset, feasible) in (
        (source, target, true),
        (target, source, false),
        (source, source, true),
        ([0.5 0.1; 0.0 0.5] * source, target, true),
        ([3.0 0.0; 0.0 1.0] * source, target, false),
        ([0.5 0.0; 0.0 0.0] * source, target, true),
        (zeros(2, 2) * source, target, true),
        # Gauge dominance holds here, but the image leaves the orthant.
        ([-0.5 0.0; 0.0 0.5] * source, target, false),
        (source, strip, true),
        (strip, strip, true),
        (origin, target, true),
        # A line must belong to the target domain in both directions.
        (strip, Sets.Intersection(strip.set, orthant), false),
        (ray, Sets.Intersection(Sets.PolarPoint([1.0, 1.0]), orthant), true),
        (source, ray, false),
    )
        model = Model(GLPK.Optimizer)
        set_silent(model)
        @constraint(model, subset ⊆ supset)
        optimize!(model)
        @test termination_status(model) == (feasible ? MOI.OPTIMAL : MOI.INFEASIBLE)
    end

    for rate in (0.6, 0.4)
        model = Model(GLPK.Optimizer)
        set_silent(model)
        @variable(model, c[1:2] >= 1)
        S = Sets.Intersection(Sets.PolarPoint(c), orthant)
        @constraint(model, ([0.5 0.0; 0.0 0.25] / rate) * S ⊆ S)
        optimize!(model)
        @test termination_status(model) == (rate == 0.6 ? MOI.OPTIMAL : MOI.INFEASIBLE)
    end

    for side in (:source, :target, :both)
        model = Model(GLPK.Optimizer)
        set_silent(model)
        @variable(model, c[1:2] >= 1)
        variable = Sets.Intersection(Sets.PolarPoint(c), orthant)
        if side == :source
            @constraint(model, variable ⊆ target)
            @objective(model, Min, sum(c))
        elseif side == :target
            @constraint(model, source ⊆ variable)
            @objective(model, Max, sum(c))
        else
            @variable(model, d[i=1:2] == i)
            other = Sets.Intersection(Sets.PolarPoint(d), orthant)
            @constraint(model, variable ⊆ other)
            @objective(model, Min, sum(c))
        end
        optimize!(model)
        @test termination_status(model) == MOI.OPTIMAL
        @test value.(c) ≈ (side == :target ? [2.0, 3.0] : [1.0, 2.0])
    end

    shifted = polyhedron(hrep([1.0 0.0; 0.0 1.0], ones(2)), lib)
    plane = polyhedron(
        hrep([HyperPlane([1.0, 0.0], 1.0)], HalfSpace{Float64,Vector{Float64}}[]),
        lib,
    )
    for domain in (shifted, plane)
        affine = Sets.Intersection(source.set, domain)
        for (subset, supset) in ((source, affine), (affine, source))
            model = Model()
            @test_throws ArgumentError @constraint(model, subset ⊆ supset)
        end
    end
end

end
