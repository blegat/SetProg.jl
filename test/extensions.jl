using Test

@testset "SumOfSquares extension" begin
    project = dirname(Base.active_project())
    core_then_sos = raw"""
    using Test, LinearAlgebra, SetProg
    const MOI = SetProg.MOI
    polynomial_packages = (
        "SumOfSquares", "DynamicPolynomials", "MultivariatePolynomials",
        "MultivariateBases", "MultivariateMoments",
    )
    function check_core_only()
        @test Base.get_extension(SetProg, :SetProgSumOfSquaresExt) === nothing
        @test all(pkg.name ∉ polynomial_packages for pkg in keys(Base.loaded_modules))
    end
    check_core_only()
    err = ArgumentError("Polynomial sets and operations require SumOfSquares. Load it with `using SumOfSquares`.")
    @test_throws err PolySet(degree = 4, symmetric = true)
    @test_throws err SetProg.new_polyvars(2)
    for n in 1:3
        mock = MOI.Utilities.MockOptimizer(MOI.Utilities.Model{Float64}())
        MOI.Utilities.set_mock_optimize!(mock, m -> MOI.Utilities.mock_optimize!(
            m, zeros(MOI.get(m, MOI.NumberOfVariables())),
        ))
        model = direct_model(mock)
        @variable(model, E, Ellipsoid(symmetric = true, dimension = n))
        @constraint(model, ones(n) in E)
        @objective(model, Min, L1_heuristic(volume(E), ones(n)))
        optimize!(model)
        @test termination_status(model) == MOI.OPTIMAL
        @test value(E) isa SetProg.Sets.Ellipsoid
        @test SetProg.data(model).perspective_polyvar === nothing
        @test isempty(SetProg.data(model).spaces.polyvars)
        @test num_constraints(model, Vector{AffExpr}, MOI.PositiveSemidefiniteConeTriangle) +
            num_constraints(model, Vector{VariableRef}, MOI.PositiveSemidefiniteConeTriangle) == 1
    end
    model = Model()
    @variable(model, P, Polytope(symmetric = true, dimension = 2))
    SetProg.set_space(SetProg.data(model), model)
    SetProg.create_spaces(SetProg.data(model))
    SetProg.load(model, P)
    @test P.variable isa SetProg.Sets.PolarPoint
    mock = MOI.Utilities.MockOptimizer(MOI.Utilities.Model{Float64}())
    MOI.Utilities.set_mock_optimize!(mock, m -> MOI.Utilities.mock_optimize!(
        m, zeros(MOI.get(m, MOI.NumberOfVariables())),
    ))
    model = direct_model(mock)
    @variable(model, S, Ellipsoid(psd = false, dimension = 2))
    initial = SetProg.Sets.ShiftedEllipsoid(Symmetric([1.0 0; 0 1]), zeros(2), -16.0)
    inclusion = @constraint(model, initial ⊆ S)
    @constraint(model, [S, initial] in Disjoint())
    optimize!(model)
    @test num_variables(model) == 8
    @test num_constraints(model, Vector{AffExpr}, MOI.PositiveSemidefiniteConeTriangle) == 2
    @test size(SetProg.certificate_matrix(inclusion)) == (3, 3)
    @test value(S) isa SetProg.Sets.ShiftedEllipsoid
    @test SetProg.sublevel_eval(value(S), [1.0, 1.0]) == 1
    @test SetProg.data(model).perspective_polyvar === nothing
    @test isempty(SetProg.data(model).spaces.polyvars)
    check_core_only()
    using DynamicPolynomials
    @test Base.get_extension(SetProg, :SetProgSumOfSquaresExt) === nothing
    @test_throws err PolySet(degree = 4)
    using SumOfSquares
    @test Base.get_extension(SetProg, :SetProgSumOfSquaresExt) !== nothing
    @polyvar x[1:2]
    @test PolySet(degree = 4, variables = x).variables == x
    @test PolySet(degree = 4).basis == Monomial
    @test length(SetProg.new_polyvars(2)) == 2
    """
    sos_then_core = raw"""
    using Test, SumOfSquares, SetProg
    @test Base.get_extension(SetProg, :SetProgSumOfSquaresExt) !== nothing
    @test PolySet(degree = 2).degree == 2
    """
    for code in (core_then_sos, sos_then_core)
        @test success(
            pipeline(
                `$(Base.julia_cmd()) --startup-file=no --project=$project -e $code`,
                stdout = stdout,
                stderr = stderr,
            ),
        )
    end
end
