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

@testset "Domain variable wrapper" begin
    lib = Polyhedra.DefaultLibrary{Float64}(
        optimizer_with_attributes(
            GLPK.Optimizer,
            "presolve" => GLPK.GLP_ON,
            MOI.Silent() => true,
        ),
    )
    orthant = polyhedron(hrep(-Matrix{Float64}(I, 2, 2), zeros(2)), lib)
    @test Domain(Polytope(symmetric = true), orthant).dimension == 2
    @test Domain(Polytope(symmetric = true), hrep(orthant)).dimension == 2
    @test Domain(Polytope(symmetric = true), MOI.Nonnegatives(3)).dimension == 3
    @test Domain(Polytope(symmetric = true), MOI.Nonnegatives).dimension === nothing
    @test Domain(Polytope(true, 2, nothing), MOI.Nonnegatives).dimension == 2
    @test_throws DimensionMismatch Domain(Polytope(dimension = 3), orthant)
    @test_throws DimensionMismatch Domain(Polytope(dimension = 2), MOI.Nonnegatives(3))

    unit = Sets.Intersection(Sets.PolarPoint(ones(2)), orthant)
    for domain in (orthant, hrep(orthant), MOI.Nonnegatives(2), MOI.Nonnegatives),
        rate in (0.6, 0.4)

        model = Model(GLPK.Optimizer)
        set_silent(model)
        # The unit simplex supplies the dimension for the deferred orthant.
        @variable(model, S, Domain(Polytope(symmetric = true), domain))
        @constraint(model, S ⊆ unit)
        @constraint(model, ([0.5 0.1; 0.0 0.25] / rate) * S ⊆ S)
        optimize!(model)
        @test termination_status(model) == (rate == 0.6 ? MOI.OPTIMAL : MOI.INFEASIBLE)
        if rate == 0.6
            fitted = value(S)
            @test fitted isa Sets.Intersection
            @test all(>=(1 - 1e-7), fitted.set.a)
            @test Sets.dimension(fitted) == 2
            @test [-1.0, 1.0] ∉ fitted.domain
            @test SetProg.data(model).space == SetProg.PrimalSpace
        end
    end

    # The domain may also be bounded and affine; membership needs no conic certificate.
    bounded = polyhedron(hrep([-1.0 0.0; 0.0 -1.0; 1.0 0.0], [0.0, 0.0, 0.5]), lib)
    for (point, feasible) in
        (([0.25, 0.25], true), ([-0.25, 0.25], false), ([0.75, 0.0], false))
        model = Model(GLPK.Optimizer)
        set_silent(model)
        @variable(model, S, Domain(Polytope(symmetric = true), bounded))
        @constraint(model, point in S)
        optimize!(model)
        @test termination_status(model) == (feasible ? MOI.OPTIMAL : MOI.INFEASIBLE)
        if feasible
            @test value(S).domain === bounded
        end
    end

    model = Model(GLPK.Optimizer)
    @variable(model, S, Domain(Polytope(symmetric = true, dimension = 2), MOI.Nonnegatives))
    @constraint(model, S ⊆ hrep(Matrix{Float64}(I, 2, 2), ones(2)))
    @test_throws ArgumentError("domain-restricted variables require primal modeling") optimize!(
        model,
    )
end

@testset "Domain forwards ellipsoid construction" begin
    model = Model()
    family = Ellipsoid(symmetric = true, dimension = 2)
    domain = polyhedron(hrep(-Matrix{Float64}(I, 2, 2), zeros(2)))
    wrapper = Domain(family, domain)
    @test wrapper.set === family
    @test wrapper.dimension == 2
    @test Sets.space_variables(wrapper) === nothing
    @test SetProg.variable_polyvars(wrapper, SetProg.Spaces(), SetProg.SpaceIndex(1)) ===
          nothing
    restricted = SetProg.variable_set(model, wrapper, SetProg.PrimalSpace, 2, nothing)
    @test restricted isa Sets.Intersection{<:Any,<:Sets.Ellipsoid}
    @test restricted.domain === domain
    @test size(restricted.set.Q) == (2, 2)
    @test num_variables(model) == 3
    # SumOfSquares specializes the 2×2 PSD cone when its extension is loaded.
    plain_model = Model()
    SetProg.variable_set(plain_model, family, SetProg.PrimalSpace, 2, nothing)
    @test list_of_constraint_types(model) == list_of_constraint_types(plain_model)
    @test num_constraints(model; count_variable_in_set_constraints = true) ==
          num_constraints(plain_model; count_variable_in_set_constraints = true) ==
          1
    @test_throws DimensionMismatch Domain(family, MOI.Nonnegatives(3))
end

@testset "Polar variable families" begin
    inner = Domain(Polytope(symmetric = true), MOI.Nonnegatives)
    family = polar(inner)
    @test polar(family) === inner
    @test family.dimension === nothing
    @test polar(Ellipsoid(dimension = 3)).dimension == 3
    @test Sets.space_variables(family) === nothing
    @test SetProg.variable_polyvars(family, SetProg.Spaces(), SetProg.SpaceIndex(1)) ===
          nothing
    orthant = polyhedron(hrep(-Matrix{Float64}(I, 2, 2), zeros(2)))
    unit = polar(Sets.Intersection(Sets.PolarPoint(ones(2)), orthant))
    for rate in (0.6, 0.4)
        model = Model(GLPK.Optimizer)
        set_silent(model)
        @variable(model, S, family)
        @constraint(model, unit ⊆ S)
        @constraint(model, ([0.5 0.1; 0.0 0.25] / rate) * S ⊆ S)
        optimize!(model)
        @test termination_status(model) == (rate == 0.6 ? MOI.OPTIMAL : MOI.INFEASIBLE)
        if rate == 0.6
            fitted = value(S)
            @test fitted isa Sets.Polar{<:Any,<:Sets.Intersection}
            @test all(>=(1 - 1e-7), polar(fitted).set.a)
            @test Sets.dimension(fitted) == 2
        end
    end
    for space in (SetProg.PrimalSpace, SetProg.DualSpace)
        model = Model()
        set = SetProg.variable_set(
            model,
            polar(Ellipsoid(symmetric = true)),
            space,
            2,
            nothing,
        )
        @test (space == SetProg.PrimalSpace ? set : polar(set)) isa
              Sets.Polar{<:Any,<:Sets.Ellipsoid}
        @test num_variables(model) == 3
    end
    # A cross-polytope partition supplies the orthant fan of a weighted 1-norm.
    lib = Polyhedra.DefaultLibrary{Float64}(
        optimizer_with_attributes(
            GLPK.Optimizer,
            "presolve" => GLPK.GLP_ON,
            MOI.Silent() => true,
        ),
    )
    partition =
        polyhedron(convexhull([1.0, 0.0], [0.0, 1.0], [-1.0, 0.0], [0.0, -1.0]), lib)
    inner = Polytope(symmetric = true, piecewise = partition)
    unit = polar(
        Sets.Piecewise(
            [Sets.PolarPoint(h.a / h.β) for h in halfspaces(partition)],
            partition,
        ),
    )
    for rate in (0.6, 0.4)
        model = Model(GLPK.Optimizer)
        set_silent(model)
        @variable(model, S, polar(inner))
        @constraint(model, unit ⊆ S)
        @constraint(model, (-Matrix{Float64}(I, 2, 2)) * S ⊆ S)
        @constraint(model, ([0.5 -0.1; 0.0 0.25] / rate) * S ⊆ S)
        optimize!(model)
        @test termination_status(model) == (rate == 0.6 ? MOI.OPTIMAL : MOI.INFEASIBLE)
        if rate == 0.6
            fitted = value(S)
            @test fitted isa Sets.Polar{<:Any,<:Sets.Piecewise}
            pieces = polar(fitted).sets
            @test all(s -> any(t -> isapprox(s.a, -t.a; atol = 1e-7), pieces), pieces)
        end
    end

end


end
