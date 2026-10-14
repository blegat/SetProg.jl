using LinearAlgebra
import MultivariatePolynomials
import DynamicPolynomials
using SumOfSquares
using Test
using SetProg

@testset "apply_matrix" begin
    DynamicPolynomials.@polyvar x[1:2]
    DynamicPolynomials.@polyvar y[1:2]
    DynamicPolynomials.@polyvar z[1:1]
    Q = [
        6 5 4
        5 3 2
        4 2 1
    ]
    q = SumOfSquares.GramMatrix(Q, MultivariatePolynomials.monomials(x, 2))
    p = MultivariatePolynomials.polynomial(q)
    @testset "2x2 Float64" begin
        B = [
            2.0 3.0
            4.0 5.0
        ]
        qB = SetProg.apply_matrix(q, B, y, 2)
        @test qB isa SumOfSquares.GramMatrix{Float64}
        @test qB == p(x => B * y)
    end
    @testset "2x1 Int" begin
        A = reshape([2, 3], 2, 1)
        qA = SetProg.apply_matrix(q, A, z, 2)
        @test qA isa SumOfSquares.GramMatrix{Int}
        @test qA == p(x => A * z)
    end
end

@testset "Rectangular polynomial images" begin
    DynamicPolynomials.@polyvar x[1:2]
    q = GramMatrix(Matrix{Float64}(I, 3, 3), MultivariatePolynomials.monomials(x, 2))
    set = SetProg.Sets.polar(SetProg.Sets.ConvexPolySet(4, q, nothing))
    model = Model()
    cref = @constraint(model, [1.0 2.0] * set ⊆ [3.0 4.0] * set)
    constraint = JuMP.constraint_object(cref)
    @test constraint.set isa SumOfSquares.SOSPolynomialSet
    @test MultivariatePolynomials.nvariables(constraint.set.basis) == 1
end
