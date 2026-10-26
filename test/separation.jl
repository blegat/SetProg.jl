module TestSeparation
using Test
using SetProg, JuMP, LinearAlgebra, SumOfSquares, DynamicPolynomials
include("utilities.jl")

mock_model() = direct_model(bridged_mock(mock -> nothing))
const Shifted = SetProg.Sets.ShiftedEllipsoid

@testset "Affine quadratic separation" begin
    model = mock_model()
    @variable(model, S, Ellipsoid(psd = false, dimension = 2))
    initial = Shifted(Symmetric([1.0 0; 0 1]), zeros(2), -16.0)
    unsafe = Shifted(Symmetric([-1.0 0; 0 -1]), zeros(2), 36.0)
    @variable(model, margin >= 0)
    gamma = SProcedureMultiplier(; base_name = "gamma")
    inclusion =
        @constraint(model, initial ⊆ S, S_procedure_scaling = gamma, margin = margin)
    separation = @constraint(model, [S, unsafe] in Disjoint(), margin = margin)
    @constraint(model, [0.5 0.0; 0.0 0.5] * S ⊆ S)
    optimize!(model)
    @test num_variables(model) == 9 # six barrier coefficients, two multipliers, margin
    @test num_constraints(model, Vector{AffExpr}, MOI.PositiveSemidefiniteConeTriangle) == 3
    @test num_constraints(model, VariableRef, MOI.GreaterThan{Float64}) == 3
    @test S.variable isa Shifted{VariableRef}
    @test has_lower_bound(gamma.variable)
    @test SetProg.transformed_constraint(inclusion) !== nothing
    @test SetProg.transformed_constraint(separation) !== nothing
    @test SetProg.sublevel_eval(initial, [4.0, 0.0]) == 1
    @test SetProg.sublevel_eval(initial, [0.0, 0.0]) == -15
    @test_throws ArgumentError SProcedureMultiplier(; degree = 1)
    @test_throws ArgumentError SProcedureMultiplier(; degree = -2)
    @test_throws ErrorException @constraint(model, [S] in Disjoint())
end

@testset "Positivity can be omitted" begin
    for psd in (true, false)
        model = mock_model()
        @variable(model, S, Ellipsoid(symmetric = true, dimension = 2, psd = psd))
        @constraint(
            model,
            [S, Shifted(Symmetric([-1.0 0; 0 -1]), zeros(2), 36.0)] in Disjoint()
        )
        optimize!(model)
        # The PSD requirement on S contributes a 2-by-2 cone only when requested.
        @test num_constraints(
            model,
            Vector{VariableRef},
            SumOfSquares.PositiveSemidefinite2x2ConeTriangle,
        ) == Int(psd)
        @test num_constraints(
            model,
            Vector{AffExpr},
            MOI.PositiveSemidefiniteConeTriangle,
        ) == 1
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

@testset "Polynomial separation and affine preimages" begin
    @polyvar x[1:2]
    model = mock_model()
    @variable(model, S, PolySet(psd = false, degree = 4, variables = x))
    initial = Shifted(Symmetric([1.0 0; 0 1]), zeros(2), -16.0)
    gamma = SProcedureMultiplier(; degree = 2)
    @constraint(model, initial ⊆ S, S_procedure_scaling = gamma)
    @constraint(model, [S, initial] in Disjoint(), S_procedure_scaling = 1)
    @constraint(model, [0.5 0.0; 0.0 0.5] * S ⊆ S)
    optimize!(model)
    @test gamma.variable isa SumOfSquares.GramMatrix
    @test length(S.variable.p.basis) == 6
    MP = SumOfSquares.MP
    mapped =
        SetProg.apply_map(SetProg.Sets.LinearPreImage(S.variable, [0.5 0.0; 0.0 0.5]), x)
    p = MP.polynomial(S.variable.p)
    @test all(iszero, MP.coefficients(MP.polynomial(mapped.p) - MP.subs(p, x => 0.5x)))
    @test MP.maxdegree(MP.polynomial(gamma.variable)) == 2
end
@testset "SOS multipliers for quadratic sublevels" begin
    model = mock_model()
    @variable(model, S, Ellipsoid(psd = false, dimension = 2))
    initial = Shifted(Symmetric([1.0 0; 0 1]), zeros(2), -16.0)
    gamma = SProcedureMultiplier(; degree = 2)
    @constraint(model, initial ⊆ S, S_procedure_scaling = gamma)
    optimize!(model)
    @test gamma.variable isa SumOfSquares.GramMatrix
    @test SumOfSquares.MP.maxdegree(SumOfSquares.MP.polynomial(gamma.variable)) == 2
end

end
