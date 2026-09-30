using Test
using LinearAlgebra
using DynamicPolynomials
using SetProg, SetProg.Sets
using Polyhedra

function _stdout_string(f)
    return mktemp() do _, io
        redirect_stdout(f, io)
        seekstart(io)
        return read(io, String)
    end
end

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
        p = SetProg.GramMatrix{Float64}((i, j) -> convert(Float64, 8 - (i + j)),
                                           monomial_vector([x, y, z]))
        set = Sets.ConvexPolySet(2, p, nothing)
        el = Sets.zero_eliminate(set, 1:2)
        @test el.p.Q == 6ones(1, 1)
        el = Sets.zero_eliminate(set, 3:3)
        @test el.p.Q == [4 3; 3 2]
        el = Sets.zero_eliminate(set, 2:2)
        @test el.p.Q == [6 4; 4 2]

        @testset "Householder" begin
            p = SetProg.GramMatrix{Float64}((i, j) -> convert(Float64, 6 - (i + j)),
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
    @testset "Ellipsoid" begin
        Q = Symmetric([2.0 1.0; 1.0 3.0])
        ell = Sets.Ellipsoid(Q)
        @test Sets.space_variables(ell) === nothing
        @test Sets.polar_representation(ell) isa Sets.Polar
        @test Sets.polar(Sets.polar_representation(ell)).Q ≈ inv(Q)
        @test Sets.zero_eliminate(ell, [1]).Q == 3ones(1, 1)
        p = project(ell, [1])
        @test p isa Sets.Polar{Float64,Sets.Ellipsoid{Float64}}
        @test Sets.dimension(p) == 1
        # The projection of the ellipsoid is the zero-elimination of its polar
        @test Sets.polar(p).Q ≈ inv(Q)[1:1, 1:1]
        @test _stdout_string() do
            return Sets.print_support_function(Sets.polar(ell))
        end == "h(S, x) = 3.0*x[2]^2 + 2.0*x[1]*x[2] + 2.0*x[1]^2\n"
        @test _stdout_string() do
            return Sets.print_support_function(Sets.polar(ell), digits = nothing)
        end == "h(S, x) = 3.0*x[2]^2 + 2.0*x[1]*x[2] + 2.0*x[1]^2\n"
        @testset "Translation" begin
            c = [1.0, 2.0]
            t = Sets.Translation(ell, c)
            @test Sets.dimension(t) == 2
            @test Sets.space_variables(t) === nothing
            tp = project(t, [2])
            @test tp.c == [2.0]
            @test Sets.polar(tp.set).Q ≈ inv(Q)[2:2, 2:2]
            lifted = Sets.LiftedEllipsoid(t)
            @test Sets.dimension(lifted) == 2
            @test Sets.space_variables(lifted) === nothing
            @test Sets.perspective_variables(lifted) === nothing
            back = Sets.ellipsoid(lifted)
            @test back isa Sets.Translation{Sets.Ellipsoid{Float64}}
            @test back.set.Q ≈ Q
            @test back.c ≈ c
        end
    end
    @testset "PolarPoint" begin
        h = Sets.PolarPoint([1.0, 2.0])
        @test Sets.dimension(h) == 2
        @test Sets.space_variables(h) === nothing
        @test Sets.scaling_function(h)(1.0, 3.0) == 7.0
        @test _stdout_string() do
            return Sets.print_support_function(Sets.polar(h))
        end == "h(S, x) = 2.0*x[2] + x[1]\n"
    end
    @testset "Projection" begin
        @polyvar x y z
        p = SetProg.GramMatrix{Float64}((i, j) -> convert(Float64, i == j),
                                        monomial_vector([x, y, z]))
        set = Sets.ConvexPolySet(2, p, nothing)
        @test Sets.gauge1(set) === set.p
        proj = Sets.Projection(set, [1, 3])
        @test Sets.dimension(proj) == 2
        @test Sets.space_variables(proj) == [x, z]
        ell_proj = Sets.Projection(Sets.Ellipsoid(Symmetric(Matrix(1.0I, 3, 3))), [2])
        @test Sets.dimension(ell_proj) == 1
        @test Sets.space_variables(ell_proj) === nothing
    end
    @testset "Piecewise" begin
        □ = polyhedron(HalfSpace([1, 0], 1.0) ∩ HalfSpace([-1, 0], 1) ∩
                       HalfSpace([0, 1], 1) ∩ HalfSpace([0, -1], 1))
        ell = Sets.Ellipsoid(Symmetric([1.0 0.0; 0.0 2.0]))
        set = Sets.Piecewise([ell, ell, ell, ell], □)
        @test Sets.dimension(set) == 2
        @test Sets.space_variables(set) === nothing
        el = Sets.zero_eliminate(set, [2])
        @test Sets.dimension(el) == 1
        @test all(s -> s.Q == ones(1, 1), el.sets)
        @test all(adj -> all(iv -> length(iv[2]) == 1, adj), el.graph)
        str = _stdout_string() do
            return Sets.print_support_function(Sets.polar(set))
        end
        @test startswith(str, "h(S, x) =\n")
        @test count("2.0*x[2]^2 + x[1]^2", str) == 4
        @test occursin("if -x[2] + x[1] ≥ 0, x[2] + x[1] ≥ 0", str)
        points = Sets.Piecewise([Sets.PolarPoint(h.a) for h in halfspaces(□)], □)
        str = _stdout_string() do
            return Sets.print_support_function(Sets.polar(points), digits = nothing)
        end
        for line in [" x[1]", " -x[1]", " x[2]", " -x[2]"]
            @test occursin(line * "\n", str)
        end
        shifted = polyhedron(HalfSpace([1, 0], -1.0) ∩ HalfSpace([-1, 0], 2) ∩
                             HalfSpace([0, 1], 1) ∩ HalfSpace([0, -1], 1))
        err = ErrorException("The origin is not in the polytope")
        @test_throws err Sets.Piecewise([ell, ell, ell, ell], shifted)
    end
end
