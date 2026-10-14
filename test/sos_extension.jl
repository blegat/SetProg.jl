module TestSumOfSquaresExtension

using Test, LinearAlgebra, DynamicPolynomials, SumOfSquares, SetProg
import MultivariatePolynomials as MP
import Hypatia

const Sets = SetProg.Sets

function polynomial_constraint(cref)
    constraint = JuMP.constraint_object(cref)
    monos = SumOfSquares.MultivariateBases.keys_as_monomials(constraint.set.basis)
    return MP.polynomial(constraint.func, monos)
end

@testset "Polynomial extension coverage" begin
    @testset "Perspective synchronization" begin
        @polyvar z w
        @test SetProg.synchronize_perspective(nothing, z) === z
        @test SetProg.synchronize_perspective(z, nothing) === z
        @test SetProg.synchronize_perspective(z, z) === z
        @test_throws ArgumentError("Perspective variables do not match") SetProg.synchronize_perspective(
            z,
            w,
        )
    end

    @testset "Polynomial integration" begin
        @polyvar x y
        triangle = polyhedron(convexhull([0.0, 0.0], [1.0, 0.0], [0.0, 1.0]))
        @test SetProg.integrate(x^2 + 2x*y + 3y^2, triangle) ≈ 5 / 12
        @test SetProg.integrate_monomials(monomials([x, y], 2), triangle) ≈
              [1/12, 1/24, 1/12]
    end

    @testset "Polynomial gauge printing" begin
        @polyvar x y
        gram = GramMatrix([1.234 0.0; 0.0 2.345], monomial_vector([x, y]))
        set = Sets.ConvexPolySet(2, gram, nothing)
        for digits in (2, nothing)
            output = mktemp() do _, io
                redirect_stdout(io) do
                    Sets.print_support_function(Sets.polar(set); digits)
                end
                seekstart(io)
                read(io, String)
            end
            p = MP.polynomial(gram)
            expected = digits === nothing ? p : round(p; digits)
            @test output == "h(S, x) = " * string(expected) * "\n"
        end
    end

    @testset "Odd powers in Gram transformations" begin
        @polyvar x[1:2]
        @polyvar y[1:2]
        gram = GramMatrix(Matrix{Float64}(I, 4, 4), monomials(x, 3))
        A = [1.0 2.0; 3.0 4.0]
        mapped = SetProg.apply_matrix(gram, A, y, 3)
        @test mapped == MP.polynomial(gram)(x => A * y)
    end

    @testset "PolySet validation" begin
        @polyvar x[1:2]
        @polyvar y[1:2]
        superset = Sets.PolySet(2, GramMatrix(Matrix{Float64}(I, 2, 2), x))
        @test PolySet(degree = 2, dimension = 2, variables = x, superset = superset).variables ==
              x
        @test_throws DimensionMismatch PolySet(
            degree = 2,
            dimension = 3,
            superset = superset,
        )
        @test_throws ErrorException(
            "Space variables set does not correspond to superset space variables.",
        ) PolySet(degree = 2, variables = y, superset = superset)
        for (set, space, err) in [
            (
                PolySet(degree = 2, convex = true),
                SetProg.PrimalSpace,
                ErrorException("Non-symmetric PolySet in PrimalSpace not implemented yet"),
            ),
            (
                PolySet(degree = 2, convex = true),
                SetProg.DualSpace,
                ArgumentError(
                    "Specify a point for nonsymmetric polyset, e.g. `PolySet(point=InteriorPoint([1.0, 0.0]))",
                ),
            ),
            (
                PolySet(degree = 2, symmetric = true),
                SetProg.DualSpace,
                ErrorException("Non-convex PolySet not supported in DualSpace"),
            ),
            (
                PolySet(degree = 2),
                SetProg.PrimalSpace,
                ErrorException("Non-convex nonsymmetric PolySet not implemented yet."),
            ),
        ]
            model = Model()
            @polyvar z v[1:2]
            SetProg.data(model).perspective_polyvar = z
            @test_throws err SetProg.variable_set(model, set, space, 2, v)
        end
    end

    @testset "S-procedure scaling" begin
        @polyvar x y
        p = x^2 + y^2
        for constant_subset in (true, false)
            model = Model()
            @variable(model, t)
            subset = (p = constant_subset ? p : t * p,)
            superset = (p = constant_subset ? t * p : p,)
            constraint = SetProg.s_procedure(model, subset, superset)
            cref = add_constraint(model, constraint)
            λ = last(all_variables(model))
            expected = constant_subset ? (λ - t) * p : (t - λ) * p
            @test polynomial_constraint(cref) == expected
        end
    end

    @testset "Constant polynomial in halfspace" begin
        @polyvar x y
        gram = GramMatrix(Matrix{Float64}(I, 2, 2), monomial_vector([x, y]))
        set = Sets.ConvexPolySet(2, gram, nothing)
        model = Model()
        cref = @constraint(model, set ⊆ HalfSpace([1.0, 0.0], 1.0))
        λ = only(all_variables(model))
        @test lower_bound(λ) == 0.0
        @test polynomial_constraint(cref) == λ * (x^2 + y^2 - 1) - x + 1
    end

    @testset "Variable polynomial and halfspace offset" begin
        @polyvar x y
        model = Model()
        @variable(model, t)
        @variable(model, a[1:2])
        @variable(model, β)
        gram = GramMatrix([t 0.0; 0.0 t], monomial_vector([x, y]))
        set = Sets.ConvexPolySet(2, gram, nothing)
        cref = @constraint(model, set ⊆ HalfSpace(a, β))
        λ = last(all_variables(model))
        @test lower_bound(λ) == 0.0
        @test polynomial_constraint(cref) ==
              λ * (t * (x^2 + y^2) - 1) - a[1] * x - a[2] * y + β
    end

    @testset "Moment matrix cache" begin
        @polyvar x y
        calls = Tuple{Int,Int}[]
        monos = monomials([x, y], 0:2)
        matrix = SetProg._moment_matrix(monos) do i, j
            push!(calls, (i, j))
            return Float64((monos[i] * monos[j])(x => 2, y => 3))
        end
        @test length(calls) == length(monomials([x, y], 0:4))
        v = [mono(x => 2, y => 3) for mono in monos]
        @test Matrix(matrix.Q) ≈ v * v'
    end

    @testset "Polynomial membership" begin
        for degree in (2, 4), point in ([0.5, 0.25], [2.0, 0.0]), affine in (false, true)
            @polyvar x y
            monos = monomials([x, y], div(degree, 2))
            coefficients = degree == 2 ? [1.0, 1.0] : [1.0, 0.0, 1.0]
            gram = GramMatrix(Matrix(Diagonal(coefficients)), monos)
            set = Sets.ConvexPolySet(degree, gram, nothing)
            model =
                Model(optimizer_with_attributes(Hypatia.Optimizer, MOI.Silent() => true))
            @variable(model, a[1:2])
            member = affine ? [a[1] + 0.1, a[2] + 0.2] : a
            cref = @constraint(model, membership, member in set)
            @test JuMP.name(cref) == "membership"
            @test num_variables(model) == (degree == 2 ? 5 : 14)
            fix.(a, point; force = true)
            optimize!(model)
            if point[1] < 1
                @test termination_status(model) in (MOI.OPTIMAL, MOI.ALMOST_OPTIMAL)
                @test value.(a) ≈ point
            else
                @test termination_status(model) in (MOI.INFEASIBLE, MOI.ALMOST_INFEASIBLE)
            end
        end
    end

    @testset "Set constraint moment matrix" begin
        @polyvar x y
        model = Model(optimizer_with_attributes(Hypatia.Optimizer, MOI.Silent() => true))
        @variable(model, S, PolySet(degree = 2, symmetric = true, variables = [x, y]))
        superset =
            Sets.PolySet(2, GramMatrix(Matrix{Float64}(I, 2, 2), monomial_vector([x, y])))
        cref = @constraint(model, S ⊆ superset)
        @objective(model, Max, L1_heuristic(volume(S), ones(2)))
        optimize!(model)
        @test termination_status(model) in (MOI.OPTIMAL, MOI.ALMOST_OPTIMAL)
        @test objective_value(model) ≈ 8 / 3 atol = 1e-6
        matrix = SumOfSquares.moment_matrix(cref)
        transformed = SumOfSquares.moment_matrix(SetProg.transformed_constraint(cref))
        @test Matrix(matrix.Q) ≈ Matrix(transformed.Q)
        @test minimum(eigvals(Symmetric(Matrix(matrix.Q)))) ≥ -1e-6
    end
end

end
