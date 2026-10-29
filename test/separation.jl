module TestSeparation
using Test
using SetProg, JuMP, LinearAlgebra, SumOfSquares, DynamicPolynomials
import GLPK
import Hypatia
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

@testset "Disjoint homogeneous ellipsoids" begin
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
end

@testset "Polynomial separation" begin
    @polyvar x[1:2]
    model = mock_model()
    @variable(model, S, PolySet(psd = false, degree = 0:4, variables = x))
    initial = Shifted(Symmetric([1.0 0; 0 1]), zeros(2), -16.0)
    gamma = SProcedureMultiplier(; degree = 2)
    @constraint(model, initial ⊆ S, S_procedure_scaling = gamma)
    @constraint(model, [S, initial] in Disjoint(), S_procedure_scaling = 1)
    @constraint(model, [0.5 0.0; 0.0 0.5] * S ⊆ S)
    optimize!(model)
    @test gamma.variable isa SumOfSquares.GramMatrix
    @test length(S.variable.p.basis) == 6
    MP = SumOfSquares.MP
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

@testset "Polytope separation in vertex coordinates" begin
    lib = Polyhedra.DefaultLibrary{Float64}(
        optimizer_with_attributes(
            GLPK.Optimizer,
            MOI.Silent() => true,
            "presolve" => GLPK.GLP_ON,
        ),
    )
    partition = polyhedron(vrep([[1.0, 0], [-1.0, 0], [0.0, 1], [0.0, -1], [0.0, 0]]), lib)
    for reverse in (false, true)
        model = direct_model(
            bridged_mock(
                mock -> MOI.Utilities.mock_optimize!(
                    mock,
                    zeros(MOI.get(mock, MOI.NumberOfVariables())),
                ),
            ),
        )
        @variable(model, S, Polytope(symmetric = true, piecewise = partition))
        reference = polar(
            SetProg.Sets.Piecewise(
                [SetProg.Sets.PolarPoint(collect(h.a)) for h in halfspaces(partition)],
                partition,
            ),
        )
        @constraint(model, reference ⊆ S) # Select vertex coordinates.
        initial =
            reverse ? SetProg.Sets.Ellipsoid(Symmetric(Matrix{Float64}(I, 2, 2) / 16)) :
            Shifted(Symmetric(Matrix{Float64}(I, 2, 2)), zeros(2), -16.0)
        unsafe = Shifted(Symmetric(-Matrix{Float64}(I, 2, 2)), zeros(2), 36.0)
        c0 = @constraint(model, initial ⊆ S)
        cu =
            reverse ? @constraint(model, [unsafe, S] in Disjoint()) :
            @constraint(model, [S, unsafe] in Disjoint())
        optimize!(model)
        @test SetProg.data(model).space == SetProg.DualSpace
        @test num_constraints(model, Vector{AffExpr}, MOI.SecondOrderCone) == 4
        @test num_constraints(
            model,
            Vector{AffExpr},
            MOI.PositiveSemidefiniteConeTriangle,
        ) == 0
        @test length(SetProg.transformed_constraint(cu)) == 4
        @test size(SetProg.certificate_matrix(cu)) == (12, 12)
        @test eigmin(SetProg.certificate_matrix(cu)) == 5
        @test size(SetProg.certificate_matrix(c0)) == (4, 4)
        @test eigmin(SetProg.certificate_matrix(c0)) == -4
    end
    pieces =
        [SetProg.Sets.PolarPoint([2.0, 3.0] .* sign.(h.a)) for h in halfspaces(partition)]
    box = polar(SetProg.Sets.Piecewise(pieces, partition))
    @test SetProg.sublevel_eval(box, [1.0, 1.0]) == 0.5
    @test SetProg.sublevel_eval(box, [-2.0, 3.0]) == 1
    @test SetProg.sublevel_eval(box, [0.0, 0.0]) == 0
    # The polar also includes zero even when all listed vertices are shifted.
    # Merely putting the listed vertices in the safe ellipsoid is insufficient.
    shifted = polar(
        SetProg.Sets.Piecewise(
            [
                SetProg.Sets.PolarPoint([10.0, 0.0] + 0.1sign.(h.a)) for
                h in halfspaces(partition)
            ],
            partition,
        ),
    )
    unsafe = Shifted(Symmetric(-Matrix{Float64}(I, 2, 2)), [10.0, 0.0], -99.0)
    model =
        direct_model(bridged_mock(mock -> MOI.Utilities.mock_optimize!(mock, Float64[])))
    certificate = @constraint(model, [shifted, unsafe] in Disjoint(), margin = 0.1)
    optimize!(model)
    @test length(certificate) == 5
    @test eigmin(SetProg.certificate_matrix(certificate)) ≈ -9.1
    invalid = polyhedron(vrep([[0.0, -1], [2.0, -1], [2.0, 1], [0.0, 1]]), lib)
    inner = polar(box)
    bad_fan = SetProg.Sets.Piecewise(inner.sets, invalid, inner.pieces, inner.graph)
    @test_throws ArgumentError SetProg._check_normal_fan(bad_fan)
    @test SetProg.sublevel_eval(inner, [1.0, 1.0]) == 5
    @test SetProg.sublevel_eval(first(inner.sets), [0.0, 0.0]) == 0
    # With no other preference, the same declarations keep primal gauge
    # coordinates: initial containment is SOC, unsafe separation is linear.
    model = direct_model(
        bridged_mock(
            mock -> MOI.Utilities.mock_optimize!(
                mock,
                zeros(MOI.get(mock, MOI.NumberOfVariables())),
            ),
        ),
    )
    @variable(model, S, Polytope(symmetric = true, piecewise = partition))
    initial = Shifted(Symmetric(Matrix{Float64}(I, 2, 2)), zeros(2), -16.0)
    unsafe_shifted = Shifted(Symmetric(-Matrix{Float64}(I, 2, 2)), [1.0, 0.0], 35.0)
    c0 = @constraint(model, initial ⊆ S)
    cu = @constraint(model, [unsafe_shifted, S] in Disjoint())
    optimize!(model)
    @test SetProg.data(model).space == SetProg.PrimalSpace
    @test length(SetProg.transformed_constraint(c0)) == 4
    @test eigmin(SetProg.certificate_matrix(c0)) == 1
    @test eigmin(SetProg.certificate_matrix(cu)) ≈ -1.2
    # The earlier far-away safe ellipsoid excludes zero, hence cannot contain
    # any origin-containing gauge sublevel set.
    model =
        direct_model(bridged_mock(mock -> MOI.Utilities.mock_optimize!(mock, Float64[])))
    impossible = @constraint(model, [inner, unsafe] in Disjoint(), margin = 0.1)
    optimize!(model)
    @test eigmin(SetProg.certificate_matrix(impossible)) == -1
end

@testset "Quadratic separation on cones" begin
    lib = Polyhedra.DefaultLibrary{Float64}(GLPK.Optimizer)
    cone = polyhedron(hrep([-1.0;;], [0.0]), lib)
    for sign in (-1, 1)
        model = Model(Hypatia.Optimizer)
        set_silent(model)
        source = Shifted(Symmetric(zeros(1, 1)), zeros(1), 0.0)
        target = Shifted(Symmetric(zeros(1, 1)), [sign / 2], 0.0)
        constraint = @constraint(
            model,
            SetProg.Sets.Intersection(source, cone) ⊆ target,
            S_procedure_scaling = 1
        )
        optimize!(model)
        @test is_solved_and_feasible(model) == (sign == -1)
        if sign == -1
            @test eigmin(Symmetric(SetProg.certificate_matrix(constraint))) >= -1e-6
        end
    end
end

@testset "Quadratic domain interfaces" begin
    lib = Polyhedra.DefaultLibrary{Float64}(GLPK.Optimizer)
    cone = polyhedron(hrep([-1.0;;], [0.0]), lib)
    zero = Shifted(Symmetric(zeros(1, 1)), zeros(1), 0.0)
    source = SetProg.Sets.Intersection(zero, cone)
    model = Model(Hypatia.Optimizer)
    set_silent(model)
    @variable(model, S, Ellipsoid(psd = false, dimension = 1))
    c = @constraint(model, source ⊆ S, S_procedure_scaling = 1)
    optimize!(model)
    @test is_solved_and_feasible(model)
    @test SetProg.data(model).space == SetProg.PrimalSpace
    @test eigmin(Symmetric(SetProg.certificate_matrix(c))) >= -1e-6
    unsafe = Shifted(Symmetric(zeros(1, 1)), [0.5], 1.0)
    for sets in ([source, unsafe], [unsafe, source])
        model = Model(Hypatia.Optimizer)
        set_silent(model)
        c = @constraint(model, sets in Disjoint(), margin = 0)
        optimize!(model)
        @test is_solved_and_feasible(model)
        @test eigmin(Symmetric(SetProg.certificate_matrix(c))) >= -1e-6
    end
    translated = polyhedron(hrep([-1.0;;], [1.0]), lib)
    model = Model()
    @test_throws ArgumentError @constraint(
        model,
        SetProg.Sets.Intersection(zero, translated) ⊆ unsafe,
        S_procedure_scaling = 1
    )
end

@testset "Piecewise quadratic separation" begin
    lib = Polyhedra.DefaultLibrary{Float64}(GLPK.Optimizer)
    partition = polyhedron(hrep([1.0 1.0; 1 -1; -1 1; -1 -1], ones(4)), lib)
    pieces = [
        SetProg.Sets.Ellipsoid(Symmetric([h.a[1] > 0 ? 1.0 : 100.0 0.0; 0.0 1.0])) for
        h in halfspaces(partition)
    ]
    set = SetProg.Sets.Piecewise(pieces, partition)
    initial = Shifted(Symmetric(Matrix{Float64}(I, 2, 2)), [-0.5, 0.0], 0.24)
    unsafe = Shifted(Symmetric(-Matrix{Float64}(I, 2, 2)), [0.45, 0.0], 1.2375)
    for reverse in (false, true)
        model = Model(Hypatia.Optimizer)
        set_silent(model)
        c0 = @constraint(model, initial ⊆ set, margin = 1e-3)
        cu = @constraint(
            model,
            (reverse ? [unsafe, set] : [set, unsafe]) in Disjoint(),
            margin = 1e-3
        )
        outer = Shifted(Symmetric(-Matrix{Float64}(I, 2, 2)), zeros(2), -1.0)
        bound = @constraint(model, set ⊆ outer, S_procedure_scaling = ones(4))
        optimize!(model)
        @test is_solved_and_feasible(model)
        @test length(c0) == length(cu) == length(bound) == 4
        for constraint in (c0, cu, bound)
            @test eigmin(Symmetric(SetProg.certificate_matrix(constraint))) >= -1e-6
        end
    end
    # The inactive pieces fail these tests globally. Cone restrictions are
    # essential, rather than simply imposing every quadratic everywhere.
    for (source, target, disjoint) in
        ((initial, pieces[end], false), (pieces[1], unsafe, true))
        source, target = SetProg._as_shifted(source), SetProg._as_shifted(target)
        model = Model(Hypatia.Optimizer)
        set_silent(model)
        if disjoint
            @constraint(model, [source, target] in Disjoint(), margin = 1e-3)
        else
            @constraint(model, source ⊆ target, margin = 1e-3)
        end
        optimize!(model)
        @test !is_solved_and_feasible(model)
    end
    @test SetProg.sublevel_eval(set, [0.5, 0.0]) ≈ 0.25
    @test SetProg.sublevel_eval(set, [-0.5, 0.0]) ≈ 25.0
    @test_throws DimensionMismatch SetProg._piecewise_scalings([1], 4)
end

end
