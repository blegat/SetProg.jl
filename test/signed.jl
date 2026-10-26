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
    @variable(model, S, PolySet(psd = false, degree = 4, variables = x))
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
    @test all(iszero, MP.coefficients(MP.polynomial(mapped.p) - MP.subs(p, x => 0.5x)))
    fitted = value(S)
    @test SetProg.sublevel_eval(fitted, [1.0, 2.0]) ==
          MP.polynomial(fitted.p)(x => [1.0, 2.0])
end
end
