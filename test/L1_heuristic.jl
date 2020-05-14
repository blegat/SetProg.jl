using Test
using DynamicPolynomials
using SetProg
using Polyhedra

@testset "L1 heuristic" begin
    @polyvar x y
    p = 2x^2*y + 3x + 4y^2 - 5x^4*y^2 + 6x^2
    @test SetProg.rectangle_integrate(p, [2, 3]) ≈ -672

    Q = [1 2 4
         2 3 5
         4 5 6]
    q = SetProg.GramMatrix(Q, SetProg.monomials([x, y], 2))
    @test SetProg.rectangle_integrate(polynomial(q), [2, 3]) ≈ 3465.6

    subset = SetProg.Sets.PolySet(4, q)
    v_rectangle = polyhedron(convexhull([2.0, 3], [-2, -3], [-2, 3], [2, -3]))
    h_rectangle = polyhedron(HalfSpace([1.0, 0], 2) ∩ HalfSpace([-1, 0], 2) ∩
                             HalfSpace([0, 1], 3) ∩ HalfSpace([0, -1], 3))
    for rectangle in [v_rectangle, h_rectangle]
        set = SetProg.Sets.Piecewise([subset, subset, subset, subset], rectangle)
        @test SetProg.l1_integral(set, nothing) ≈ 3465.6
    end

    Q = [1 2
         2 3]
    @test SetProg.rectangle_integrate(polynomial(Q, [x, y]), [2, 3]) ≈ 248

    ell = SetProg.Sets.Ellipsoid(Symmetric(Q))
    for rectangle in [v_rectangle, h_rectangle]
        set = SetProg.Sets.Piecewise([ell, ell, ell, ell], rectangle)
        @test SetProg.l1_integral(set, nothing) ≈ 248
    end
end
