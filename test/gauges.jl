module TestGauges
using Test, LinearAlgebra, SetProg
using SumOfSquares, DynamicPolynomials
import Hypatia, GLPK

@testset "Homogeneous gauges" begin
    A = [0.5 0.0; 0.0 0.25]
    families = (
        Ellipsoid(symmetric = true),
        Polytope(piecewise = :orthant),
        Polytope(piecewise = :orthant, polar = true),
        Polytope(symmetric = true, piecewise = :orthants, polar = true),
        Polytope(symmetric = true, piecewise = [[1.0 0.0; 0.0 1.0], [-1.0 0.0; 0.0 1.0]]),
    )
    @polyvar x[1:2]
    families = (families..., PolySet(symmetric = true, degree = 4, variables = x))
    for family in families, rate in (0.6, 0.4)
        optimizer = family isa Union{Ellipsoid,PolySet} ? Hypatia.Optimizer : GLPK.Optimizer
        model = Model(optimizer_with_attributes(optimizer, MOI.Silent() => true))
        S = SetProg.add_function_variables!(model, family, 2, 1)
        SetProg.add_normalization!(model, S.set, S)
        SetProg.add_decrease!(
            model,
            family,
            S,
            S,
            A;
            scale = rate^SetProg.homogeneity_degree(family),
        )
        optimize!(model)
        if rate > 0.5
            @test termination_status(model) in (MOI.OPTIMAL, MOI.ALMOST_OPTIMAL)
            fitted = SetProg.solution_value(family, S)
            @test fitted([1.0, 1.0]) > 0
            @test fitted([2.0, 2.0]) ≈
                  2^SetProg.homogeneity_degree(family) * fitted([1.0, 1.0])
            if family isa PolySet
                @test_throws ArgumentError SetProg.decrease_slack(family, fitted, fitted, A)
            else
                @test SetProg.decrease_slack(
                    family,
                    fitted,
                    fitted,
                    A;
                    scale = rate^SetProg.homogeneity_degree(family),
                ) >= -1e-6
            end
        else
            @test termination_status(model) in (MOI.INFEASIBLE, MOI.ALMOST_INFEASIBLE)
        end
    end
end

@testset "Gauge data and duality" begin
    @test_throws ArgumentError Polytope(piecewise = :orthant, polar = true, min_scale = 0)
    @test_throws ArgumentError Polytope(piecewise = :orthant, min_scale = 0)
    @test_throws ArgumentError Polytope(piecewise = :orthants, transform = zeros(2, 2))
    @test_throws ArgumentError Polytope(piecewise = Matrix{Float64}[])
    @test_throws DimensionMismatch SetProg.with_dimension(
        Polytope(piecewise = :orthants, transform = [1.0 0.0; 0.0 1.0]),
        3,
    )
    @test SetProg.supports_dual(Ellipsoid(symmetric = true))
    @test SetProg.dual_family(Polytope(piecewise = :orthant)).polar
    @test !SetProg.dual_family(Polytope(piecewise = :orthant, polar = true)).polar
    @test SetProg.dual_family(
        SetProg.dual_family(Polytope(piecewise = :orthant, polar = true, min_scale = 2)),
    ) == Polytope(piecewise = :orthant, polar = true, min_scale = 2)
    @test SetProg.homogeneity_degree([
        Polytope(piecewise = :orthant, polar = true),
        Polytope(piecewise = :orthant),
    ]) == 1
    ell = SetProg.Sets.Ellipsoid(Symmetric([2.0 0.0; 0.0 3.0]))
    @test SetProg.dual_value(SetProg.dual_value(ell)).Q ≈ ell.Q
    simplex = SetProg._linear_piecewise([[2.0, 3.0]], [Matrix{Float64}(I, 2, 2)])
    @test SetProg.dual_value(SetProg.dual_value(simplex)).sets[1].a == simplex.sets[1].a
end

@testset "Gauge normalization and edge directions" begin
    @polyvar x[1:2]
    family = SetProg.with_dimension(PolySet(symmetric = true, degree = 6), 2; variables = x)
    @test family.variables == x
    @test SetProg.unit_set(family)([2.0, 3.0]) ≈ 13^3
    @test_throws DimensionMismatch SetProg.with_dimension(
        PolySet(symmetric = true, degree = 4, variables = x),
        3,
    )
    @test_throws ArgumentError SetProg.with_dimension(Ellipsoid(), 2)
    @test_throws ArgumentError SetProg.homogeneity_degree(SetProg.AbstractVariable[])
    @test_throws ArgumentError SetProg.homogeneity_degree([
        Ellipsoid(symmetric = true),
        Polytope(piecewise = :orthant, polar = true),
    ])
    @test_throws ArgumentError Polytope(piecewise = [zeros(2, 0)])
    @test_throws ArgumentError Polytope(piecewise = :orthants, min_scale = 0)
    model = Model(optimizer_with_attributes(Hypatia.Optimizer, MOI.Silent() => true))
    S = SetProg.add_function_variables!(model, Ellipsoid(symmetric = true), 2, 1)
    SetProg.add_normalization!(model, S.set, S; conditioning_bound = 10)
    SetProg.add_decrease!(model, S.set, S, S, zeros(2, 2); scale = 0)
    optimize!(model)
    @test termination_status(model) in (MOI.OPTIMAL, MOI.ALMOST_OPTIMAL)
    Q = value(S).Q
    @test minimum(eigvals(Q)) >= 1 - 1e-6
    @test maximum(eigvals(Q)) <= 10 + 1e-6
    @test_throws ArgumentError SetProg.add_normalization!(
        model,
        S.set,
        S;
        conditioning_bound = 0.5,
    )
    @test_throws ArgumentError SetProg.add_decrease!(
        model,
        S.set,
        S,
        S,
        zeros(2, 2);
        scale = -1,
    )
    model = Model(optimizer_with_attributes(GLPK.Optimizer, MOI.Silent() => true))
    B = SetProg.add_function_variables!(
        model,
        Polytope(piecewise = :orthant, polar = true),
        2,
        1,
    )
    @test_throws ArgumentError SetProg.add_normalization!(
        model,
        B.set,
        B;
        conditioning_bound = 10,
    )
    SetProg.add_normalization!(model, B.set, B)
    SetProg.add_decrease!(model, B.set, B, B, Matrix{Float64}(I, 2, 2); scale = 0)
    optimize!(model)
    @test termination_status(model) == MOI.INFEASIBLE
    simplex = SetProg._linear_piecewise([[2.0, 3.0]], [Matrix{Float64}(I, 2, 2)])
    box = SetProg.polar(SetProg._linear_piecewise([[4.0, 5.0]], [Matrix{Float64}(I, 2, 2)]))
    @test SetProg.decrease_slack(
        Polytope(piecewise = :orthant),
        simplex,
        SetProg._linear_piecewise([[1.0, 2.0]], [Matrix{Float64}(I, 2, 2)]),
        [0.5 0.2; 0.0 0.5],
    ) ≈ 1.5
    @test SetProg.decrease_slack(
        Polytope(piecewise = :orthant, polar = true),
        box,
        SetProg.polar(SetProg._linear_piecewise([[3.0, 4.0]], [Matrix{Float64}(I, 2, 2)])),
        [0.5 0.2; 0.0 0.5],
    ) ≈ 0
    G = [2.0 0.0; 0.0 1.0]
    signs = ([1.0, 1.0], [-1.0, 1.0], [1.0, -1.0], [-1.0, -1.0])
    parallelotope = SetProg.polar(
        SetProg._linear_piecewise(
            [inv(G) * (sign .* [1.0, 4.0]) for sign in signs],
            [G' * Diagonal(sign) for sign in signs],
        ),
    )
    @test parallelotope([3.0, 2.0]) ≈ 6
    @test !SetProg.supports_dual(
        Polytope(piecewise = :orthants, transform = G', polar = true),
    )
end


@testset "Polytope configurations allow automatic space selection" begin
    for partition in (:orthant, :orthants),
        polar in (false, true),
        space in (SetProg.PrimalSpace, SetProg.DualSpace)

        family = Polytope(symmetric = true, piecewise = partition, polar = polar)
        model = Model(optimizer_with_attributes(GLPK.Optimizer, MOI.Silent() => true))
        @variable(model, S, SetProg.with_dimension(family, 2))
        if space == SetProg.PrimalSpace
            @constraint(model, vrep([[0.1, 0.0], [0.0, 0.1]]) ⊆ S)
        else
            @constraint(
                model,
                S ⊆ hrep([1.0 0.0; 0.0 1.0; -1.0 0.0; 0.0 -1.0], fill(3.0, 4))
            )
        end
        @constraint(model, S ⊆ S)
        SetProg.add_normalization!(model, S.set, S)
        SetProg.add_decrease!(model, S.set, S, S, [0.5 0.0; 0.0 0.0]; scale = 0.6)
        optimize!(model)
        @test termination_status(model) == MOI.OPTIMAL
        @test SetProg.data(model).space == space
        @test SetProg.homogeneity_degree(value(S)) == 1
        @test value(S)([0.1, 0.0]) <= 1 + 1e-6
        if partition == :orthant
            @test value(S)([-0.1, 0.0]) == Inf
        end
    end
end


@testset "Orthant domain is preserved" begin
    model = Model(GLPK.Optimizer)
    @variable(model, S, Polytope(piecewise = :orthant, dimension = 2, polar = true))
    @constraint(model, [0.0 -0.5; 0.5 0.0] * S ⊆ S)
    @test_throws ArgumentError optimize!(model)
end


@testset "Polar inclusions retain line directions of singular images" begin
    model = Model(optimizer_with_attributes(GLPK.Optimizer, MOI.Silent() => true))
    @variable(
        model,
        S,
        Polytope(symmetric = true, piecewise = :orthants, dimension = 2, polar = true)
    )
    @constraint(model, [0.5 0.0; 0.0 0.5] * S ⊆ [1.0 0.0; 0.0 0.0] * S)
    optimize!(model)
    @test termination_status(model) == MOI.INFEASIBLE
end

end
