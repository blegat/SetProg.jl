using Test
using DynamicPolynomials
using SetProg, SetProg.Sets

@testset "Sets" begin
    @testset "zero_eliminate" begin
        @polyvar x y z
        p = SetProg.GramMatrix{Float64}((i, j) -> convert(Float64, i + j),
                                           monovec([x, y, z]))
        set = Sets.ConvexPolySet(2, p, nothing)
        el = Sets.zero_eliminate(set, 1:2)
        @test el.p.Q == 6ones(1, 1)
        el = Sets.zero_eliminate(set, 3:3)
        @test el.p.Q == [2 3; 3 4]
        el = Sets.zero_eliminate(set, 2:2)
        @test el.p.Q == [2 4; 4 6]
    end
end
