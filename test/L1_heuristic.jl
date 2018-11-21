using Test
using DynamicPolynomials
using SetProg

@testset "L1 heuristic" begin
    @polyvar x y
    p = 2x^2*y + 3x + 4y^2 - 5x^4*y^2 + 6x^2
    @test SetProg.rectangle_integrate(p, [2, 3]) ≈ -672
end
