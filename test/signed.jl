module TestSigned
using Test
using SetProg, JuMP, LinearAlgebra, SumOfSquares, DynamicPolynomials
include("utilities.jl")

mock_model() = direct_model(bridged_mock(mock -> nothing))
function fitted_model()
    return direct_model(
        bridged_mock(
            mock -> MOI.Utilities.mock_optimize!(
                mock,
                collect(1.0:MOI.get(mock, MOI.NumberOfVariables())),
            ),
        ),
    )
end

@testset "Signed affine quadratics" begin
    model = fitted_model()
    @variable(model, S, Ellipsoid(psd = false, dimension = 2))
    optimize!(model)
    @test num_variables(model) == 6
    @test num_constraints(model; count_variable_in_set_constraints = true) == 0
    @test S.variable isa SetProg.Sets.ShiftedEllipsoid{VariableRef}
    @test SetProg.Sets.space_variables(S.variable) === nothing
    @test SetProg.Sets.perspective_variable(S.variable) === nothing
    fitted = value(S)
    @test SetProg.sublevel_eval(fitted, [1.0, 2.0]) == 56
    A = [0.5 0.0; 0.0 0.25]
    mapped = SetProg.apply_map(SetProg.Sets.LinearPreImage(S.variable, A), nothing)
    @test SetProg.sublevel_eval(value(mapped), [1.0, 2.0]) ==
          SetProg.sublevel_eval(fitted, A * [1.0, 2.0])
end

@testset "Positivity can be omitted" begin
    for psd in (true, false)
        model = mock_model()
        @variable(model, S, Ellipsoid(symmetric = true, dimension = 2, psd = psd))
        optimize!(model)
        # The PSD requirement on S contributes a 2-by-2 cone only when requested.
        @test num_constraints(
            model,
            Vector{VariableRef},
            SumOfSquares.PositiveSemidefinite2x2ConeTriangle,
        ) == Int(psd)

    end
    model = mock_model()
    @variable(model, S, Ellipsoid(symmetric = true, dimension = 3, psd = false))
    unit = SetProg.Sets.Ellipsoid(Symmetric(Matrix{Float64}(I, 3, 3)))
    @constraint(model, SetProg.Polyhedra.polar(unit) ⊆ S)
    optimize!(model)
    @test SetProg.data(model).space == SetProg.DualSpace
    @test num_constraints(model, Vector{AffExpr}, MOI.PositiveSemidefiniteConeTriangle) == 1
    @test num_constraints(
        model,
        Vector{VariableRef},
        MOI.PositiveSemidefiniteConeTriangle,
    ) == 0
    model = mock_model()
    @variable(
        model,
        S,
        Ellipsoid(
            symmetric = true,
            psd = false,
            superset = SetProg.Sets.Ellipsoid(Symmetric(Matrix{Float64}(I, 2, 2))),
        )
    )
    optimize!(model)
    # Omitting intrinsic positivity must retain the requested reference bound.
    @test num_constraints(
        model,
        Vector{VariableRef},
        SumOfSquares.PositiveSemidefinite2x2ConeTriangle,
    ) == 1
    @polyvar x[1:2]
    reference = SetProg.Sets.PolySet(
        4,
        SumOfSquares.GramMatrix(Matrix{Float64}(I, 3, 3), SumOfSquares.MP.monomials(x, 2)),
    )
    model = mock_model()
    @variable(
        model,
        S,
        PolySet(symmetric = true, psd = false, degree = 4, superset = reference)
    )
    optimize!(model)
    @test num_constraints(
        model,
        Vector{VariableRef},
        MOI.PositiveSemidefiniteConeTriangle,
    ) == 1
    for psd in (true, false)
        model = mock_model()
        @variable(model, S, PolySet(symmetric = true, degree = 4, variables = x, psd = psd))
        # No additional constraints are necessary to choose primal modeling.
        optimize!(model)
        @test num_constraints(
            model,
            Vector{VariableRef},
            MOI.PositiveSemidefiniteConeTriangle,
        ) == Int(psd)
    end
end

@testset "Signed affine polynomials and preimages" begin
    @polyvar x[1:2]
    model = fitted_model()
    @variable(model, S, PolySet(psd = false, degree = 0:4, variables = x))
    optimize!(model)
    @test num_variables(model) == 21
    @test num_constraints(model; count_variable_in_set_constraints = true) == 0
    @test length(S.variable.p.basis) == 6
    MP = SumOfSquares.MP
    p = MP.polynomial(S.variable.p)
    @test MP.mindegree(p) == 0
    @test MP.maxdegree(p) == 4
    mapped =
        SetProg.apply_map(SetProg.Sets.LinearPreImage(S.variable, [0.5 0.0; 0.0 0.5]), x)
    @test mapped.degree == 0:4
    @test all(iszero, MP.coefficients(MP.polynomial(mapped.p) - MP.subs(p, x => 0.5x)))
    fitted = value(S)
    @test SetProg.sublevel_eval(fitted, [1.0, 2.0]) ==
          MP.polynomial(fitted.p)(x => [1.0, 2.0])
end
@testset "Explicit polynomial degrees are independent of positivity" begin
    @polyvar x[1:2]
    MP = SumOfSquares.MP
    for psd in (true, false), degree in (4, 0:4)
        model = mock_model()
        family = PolySet(; degree, psd, symmetric = degree isa Int, variables = x)
        @variable(model, S, family)
        optimize!(model)
        @test S.variable.degree == degree
        @test MP.mindegree(MP.polynomial(S.variable.p)) == minimum(degree)
        @test MP.maxdegree(MP.polynomial(S.variable.p)) == maximum(degree)
        @test length(S.variable.p.basis) == (degree isa Int ? 3 : 6)
        @test num_constraints(
            model,
            Vector{VariableRef},
            MOI.PositiveSemidefiniteConeTriangle,
        ) == Int(psd)
        @test num_constraints(model, AffExpr, MOI.EqualTo{Float64}) == 0
    end
    # These selections require coefficients outside the requested degrees to
    # vanish, even though the enclosing Gram basis can produce them.
    for (degree, symmetric, excluded) in
        ((0:2:4, false, 6), (0:4, true, 6), (1:4, false, 1), ([0, 4], false, 9))
        model = mock_model()
        @variable(model, S, PolySet(; degree, symmetric, psd = false, variables = x))
        optimize!(model)
        @test S.variable.degree == degree
        @test num_constraints(model, AffExpr, MOI.EqualTo{Float64}) == excluded
    end
    model = mock_model()
    @variable(model, S, PolySet(degree = 2:4, psd = false, variables = x))
    optimize!(model)
    @test length(S.variable.p.basis) == 5
    @test MP.mindegree(MP.polynomial(S.variable.p)) == 2
    @test num_constraints(model; count_variable_in_set_constraints = true) == 0
    for degree in (Int[], -2:4, 0:3)
        @test_throws ArgumentError PolySet(; degree, variables = x)
    end
end

@testset "Nonhomogeneous convex polynomial sets" begin
    @polyvar x[1:2]
    model = mock_model()
    @variable(model, S, PolySet(degree = 0:4, convex = true, variables = x))
    optimize!(model)
    @test S.variable isa SetProg.Sets.ConvexPolySet
    @test S.variable.degree == 0:4
    mapped = SetProg.apply_map(
        SetProg.Sets.LinearPreImage(S.variable, 0.5 * Matrix{Float64}(I, 2, 2)),
        x,
    )
    @test mapped.degree == 0:4
    # A lifted moment certificate for a fixed, nonhomogeneous convex set.
    gram =
        SumOfSquares.GramMatrix(Matrix{Float64}(I, 6, 6), SumOfSquares.MP.monomials(x, 0:2))
    fixed = SetProg.Sets.ConvexPolySet(0:4, gram, nothing)
    model = mock_model()
    @variable(model, point[1:2])
    @constraint(model, point in fixed)
    optimize!(model)
    @test num_variables(model) == 14
end

@testset "Fixed polynomial sublevel sets" begin
    @polyvar x[1:2]
    for p in (sum(x .^ 2)^2, -(1 + sum(x .^ 2))^2)
        fixed = SetProg.Sets.PolySet(p)
        @test fixed.degree == (iszero(p(x => zeros(2))) ? 4 : 0:4)
        @test SetProg.Sets.space_variables(fixed) == x
        @test SetProg.sublevel_eval(fixed, [2.0, 3.0]) ≈ p(x => [2.0, 3.0])
    end
    for degree in (4, 0:4)
        model = fitted_model()
        reference = SetProg.Sets.PolySet(degree, sum(x .^ 2)^2)
        @variable(
            model,
            S,
            PolySet(symmetric = true, degree = degree, superset = reference)
        )
        optimize!(model)
        @test SetProg.sublevel_eval(value(S), [1.0, 2.0]) ≈
              value(SetProg.sublevel_eval(S.variable, [1.0, 2.0]))
        @test num_constraints(
            model,
            Vector{VariableRef},
            MOI.PositiveSemidefiniteConeTriangle,
        ) == 1
        mapped = SetProg.apply_map(
            SetProg.Sets.LinearPreImage(S.variable, 0.5 * Matrix{Float64}(I, 2, 2)),
            x,
        )
        @test SetProg.sublevel_eval(value(mapped), [2.0, 4.0]) ≈
              SetProg.sublevel_eval(value(S), [1.0, 2.0])
    end
end

end
