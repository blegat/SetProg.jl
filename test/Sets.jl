using Test
using LinearAlgebra
using DynamicPolynomials
using SetProg, SetProg.Sets

@testset "Sets" begin
    @testset "ConvexPolySet" begin
        @polyvar x y
        P = SetProg.SumOfSquares.SymMatrix(Float64[1, 2, 3], 2)
        Q = SetProg.SumOfSquares.SymMatrix(BigInt[2, 3, 4], 2)
        basis = SetProg.Sets.MonoBasis(monomial_vector([x, y]))
        q = SetProg.Sets.ConvexPolySet(2, SetProg.SumOfSquares.GramMatrix(P, basis), Q)
        @test q isa SetProg.Sets.ConvexPolySet{BigFloat}
    end
    @testset "zero_eliminate" begin
        @polyvar x y z
        p = SetProg.GramMatrix{Float64}((i, j) -> convert(Float64, i + j),
                                           monomial_vector([x, y, z]))
        set = Sets.ConvexPolySet(2, p, nothing)
        el = Sets.zero_eliminate(set, 1:2)
        @test el.p.Q == 6ones(1, 1)
        el = Sets.zero_eliminate(set, 3:3)
        @test el.p.Q == [2 3; 3 4]
        el = Sets.zero_eliminate(set, 2:2)
        @test el.p.Q == [2 4; 4 6]

        @testset "Householder" begin
            p = SetProg.GramMatrix{Float64}((i, j) -> convert(Float64, i + j),
                                            monomial_vector([x, y]))
            set = SetProg.perspective_dual_polyset(2, p, SetProg.InteriorPoint(zeros(2)), z, [x, y])
            @test set.set.p == 2x^2 + 6x*y + 4y^2 - z^2
            @test set.set.h == zeros(2)
            @test set.set.x == [x, y]
            @test Sets.gauge1(set.set.set) == 2x^2 + 6x*y + 4y^2
            set2 = Sets.project(set, [2])
            @test set2.set.p == 4y^2 - z^2
            @test set2.set.h == zeros(1)
            @test set2.set.x == [y]
        end
    end
end
