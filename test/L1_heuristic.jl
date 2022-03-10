using Test, LinearAlgebra
using DynamicPolynomials
using SetProg
using Polyhedra

@testset "L1 heuristic" begin
    @polyvar x y
    vars = [x, y]
    p = 2x^2*y + 3x + 4y^2 - 5x^4*y^2 + 6x^2
    @test SetProg.rectangle_integrate(p, vars, [2, 3]) ≈ -672

    Q = [1 2 4
         2 3 5
         4 5 6]
    q = SetProg.GramMatrix(Q, SetProg.monomials([x, y], 2))
    @test SetProg.rectangle_integrate(polynomial(q), vars, [2, 3]) ≈ 3465.6

    @test SetProg.rectangle_integrate(polynomial(q), vars, [1, 1]) ≈ 472/45

    subset = SetProg.Sets.PolySet(4, q)
    v_rep(a, b) = polyhedron(convexhull([a, b], [-a, -b], [-a, b], [a, -b]))
    h_rep(a, b) = polyhedron(HalfSpace([1, 0], a) ∩ HalfSpace([-1,  0], a) ∩
                             HalfSpace([0, 1], b) ∩ HalfSpace([ 0, -1], b))
    v_square = v_rep(1.0, 1.0)
    h_square = h_rep(1.0, 1.0)
    for square in [v_square, h_square]
        set = SetProg.Sets.Piecewise([subset, subset, subset, subset], square)
        @test SetProg.l1_integral(set, nothing) ≈ 472/45
    end

    Q = [1 2
         2 3]
    @test SetProg.rectangle_integrate(polynomial(Q, [x, y]), vars, [2, 3]) ≈ 248

    @test SetProg.rectangle_integrate(polynomial(Q, [x, y]), vars, [1, 1]) ≈ 16/3

    ell = SetProg.Sets.Ellipsoid(Symmetric(Q))
    for square in [v_square, h_square]
        set = SetProg.Sets.Piecewise([ell, ell, ell, ell], square)
        @test SetProg.l1_integral(set, nothing) ≈ 16/3
    end

    Δ = polyhedron(convexhull([0.0, 0.0], [1.0, 0.0], [1.0, 0.5]))
    Q = [0.0 0.5
         0.5 0.0]
    ell = SetProg.Sets.Ellipsoid(Symmetric(Q))
    @test SetProg.l1_integral(ell, Δ) ≈ 1/32
    Q = [1.0 0.0
         0.0 0.0]
    ell = SetProg.Sets.Ellipsoid(Symmetric(Q))
    @test SetProg.l1_integral(ell, Δ) ≈ 1/8
    Q = [0.0 0.0
         0.0 1.0]
    ell = SetProg.Sets.Ellipsoid(Symmetric(Q))
    @test SetProg.l1_integral(ell, Δ) ≈ 1/96
end
