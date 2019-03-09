using Test
using SetProg
const MOIT = MOI.Test

function _model(optimizer::MOI.AbstractOptimizer)
    MOI.empty!(optimizer)
    return direct_model(optimizer)
end

function _model(factory::OptimizerFactory)
    return Model(factory)
end

function bridged_mock(mock_optimize!::Function...;
                      model = JuMP._MOIModel{Float64}())
    mock = MOI.Utilities.MockOptimizer(model)
    bridged = MOI.Bridges.full_bridge_optimizer(mock, Float64)
    MOI.Utilities.set_mock_optimize!(mock, mock_optimize!...)
    return bridged
end

@testset "apply_matrix" begin
    SetProg.@polyvar x[1:2]
    SetProg.@polyvar y[1:2]
    SetProg.@polyvar z[1:1]
    Q = [1 2 4
         2 3 5
         4 5 6]
    q = SetProg.GramMatrix(Q, SetProg.monomials(x, 2))
    p = SetProg.polynomial(q)
    @testset "2x2 Float64" begin
        B = [2.0 3.0
             4.0 5.0]
        qB = SetProg.apply_matrix(q, B, y, 2)
        @test qB isa SetProg.GramMatrix{Float64}
        @test qB == p(x => B * y)
    end
    @testset "2x1 Int" begin
        A = reshape([2, 3], 2, 1)
        qA = SetProg.apply_matrix(q, A, z, 2)
        @test qA isa SetProg.GramMatrix{Int}
        @test qA == p(x => A * z)
    end
end
