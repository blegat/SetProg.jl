module TestGauges
using Test, LinearAlgebra, SetProg
import GLPK

@testset "Polytope configurations allow automatic space selection" begin
    for partition in (:orthant, :orthants),
        polar in (false, true),
        space in (SetProg.PrimalSpace, SetProg.DualSpace)

        family = Polytope(symmetric = true, piecewise = partition, polar = polar)
        model = Model(optimizer_with_attributes(GLPK.Optimizer, MOI.Silent() => true))
        @variable(
            model,
            S,
            Polytope(symmetric = true, piecewise = partition, polar = polar, dimension = 2)
        )
        if space == SetProg.PrimalSpace
            @constraint(model, vrep([[0.1, 0.0], [0.0, 0.1]]) ⊆ S)
        else
            @constraint(
                model,
                S ⊆ hrep([1.0 0.0; 0.0 1.0; -1.0 0.0; 0.0 -1.0], fill(3.0, 4))
            )
        end
        @constraint(model, S ⊆ S)
        @constraint(model, ([0.5 0.0; 0.0 0.0] / 0.6) * S ⊆ S)
        optimize!(model)
        @test termination_status(model) == MOI.OPTIMAL
        @test SetProg.data(model).space == space
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
