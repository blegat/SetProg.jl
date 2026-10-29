module TestLpSum

using Test, LinearAlgebra, SetProg
import GLPK, Hypatia
const Sets = SetProg.Sets
const E = Sets.Ellipsoid(Symmetric([2.0 0.25; 0.25 1.0]))
const F = Sets.Ellipsoid(Symmetric([1.0 0.0; 0.0 3.0]))
const BALL = Sets.Ellipsoid(Symmetric(Matrix{Float64}(I, 2, 2)))

@testset "Lₚ sum metadata and polarity" begin
    for operation in (LpSum, PolarLpSum)
        set = operation(E, F, 2)
        @test Sets.dimension(set) == 2
        @test Sets.space_variables(set) === nothing
        @test Sets.perspective_variable(set) === nothing
        @test !SetProg.needs_polyvars(set)
        @test !SetProg.need_variablify(set)
        opposite = polar(set)
        @test opposite isa (operation === LpSum ? PolarLpSum : LpSum)
        @test polar(opposite).first === E
        @test polar(opposite).second === F
        @test polar(opposite).p == 2
        for p in (-1, 0, 0.5, Inf, NaN)
            @test_throws ArgumentError operation(E, F, p)
        end
        @test_throws DimensionMismatch operation(E, Sets.PolarPoint([1.0]), 2)
    end
end

@testset "Powered gauges and support functions" begin
    for p in (1, 1.5, 2, 4), x in ([1.0, 0.0], [0.5, -2.0], zeros(2))
        sum = LpSum(E, F, p)
        inverse_sum = PolarLpSum(E, F, p)
        @test support_function(x, sum)^p ≈ dot(x, E.Q \ x)^(p / 2) + dot(x, F.Q \ x)^(p / 2)
        @test SetProg.sublevel_eval(inverse_sum, x) ≈
              dot(x, E.Q * x)^(p / 2) + dot(x, F.Q * x)^(p / 2)
        # Polarity changes support addition into gauge addition.
        @test support_function(x, polar(inverse_sum))^p ≈
              SetProg.sublevel_eval(inverse_sum, x)
        nested = PolarLpSum(inverse_sum, E, p)
        @test SetProg.sublevel_eval(nested, x) ≈
              2dot(x, E.Q * x)^(p / 2) + dot(x, F.Q * x)^(p / 2)
    end
    x = [2.0, -1.0]
    @test support_function(x, LpSum(BALL, BALL, 1)) ≈ 2norm(x)
    @test support_function(x, LpSum(BALL, BALL, 2)) ≈ sqrt(2) * norm(x)
    @test SetProg.sublevel_eval(PolarLpSum(E, F, 2), x) ≈ dot(x, (E.Q + F.Q) * x)
    @test SetProg.sublevel_eval(PolarLpSum(PolarLpSum(BALL, BALL, 1), BALL, 2), x) ≈
          5dot(x, x)
end

const LIB = Polyhedra.DefaultLibrary{Float64}(
    optimizer_with_attributes(GLPK.Optimizer, "presolve" => GLPK.GLP_ON),
)
const BOX = polyhedron(hrep([1.0 0.0; 0.0 1.0; -1.0 0.0; 0.0 -1.0], ones(4)), LIB)
const DIAMOND = polyhedron(hrep([1.0 1.0; 1.0 -1.0; -1.0 1.0; -1.0 -1.0], ones(4)), LIB)
linear_pieces(partition) =
    Sets.Piecewise([Sets.PolarPoint(h.a) for h in halfspaces(partition)], partition)

@testset "Polyhedral and piecewise quadratic gauges" begin
    box = linear_pieces(BOX)
    diamond = linear_pieces(DIAMOND)
    piecewise = Sets.Piecewise(
        [
            Sets.Ellipsoid(
                Symmetric(iszero(h.a[1]) ? [1.0 0.0; 0.0 2.0] : [2.0 0.0; 0.0 1.0]),
            ) for h in halfspaces(BOX)
        ],
        BOX,
    )
    for x in ([2.0, 0.5], [0.5, -2.0], [-1.0, 1.0])
        @test SetProg.sublevel_eval(PolarLpSum(BALL, box, 2), x) ≈
              dot(x, x) + maximum(abs, x)^2
        @test SetProg.sublevel_eval(PolarLpSum(diamond, box, 1), x) ≈
              sum(abs, x) + maximum(abs, x)
        @test support_function(x, LpSum(box, diamond, 1)) ≈ sum(abs, x) + maximum(abs, x)
        @test SetProg.sublevel_eval(PolarLpSum(BALL, piecewise, 2), x) ≈
              2dot(x, x) + maximum(abs2, x)
    end
    orthant = polyhedron(hrep(-Matrix{Float64}(I, 2, 2), zeros(2)), LIB)
    cost = Sets.Intersection(Sets.PolarPoint([1.0, 2.0]), orthant)
    @test SetProg.sublevel_eval(PolarLpSum(cost, box, 1), [2.0, 1.0]) == 6
    @test isinf(SetProg.sublevel_eval(PolarLpSum(cost, box, 1), [-2.0, 1.0]))
    # A halfspace's gauge is zero, rather than negative, in recession directions.
    @test SetProg.sublevel_eval(
        PolarLpSum(Sets.PolarPoint([1.0, 0.0]), box, 2),
        [-2.0, 1.0],
    ) == 4
end

@testset "Set-variable expressions preserve space inference" begin
    for space in (SetProg.PrimalSpace, SetProg.DualSpace)
        mock = MOI.Utilities.MockOptimizer(MOI.Utilities.Model{Float64}())
        MOI.Utilities.set_mock_optimize!(
            mock,
            m -> MOI.Utilities.mock_optimize!(
                m,
                [
                    i % 3 == 1 ? 2.0 : (i % 3 == 2 ? 0.0 : 3.0) for
                    i = 1:MOI.get(m, MOI.NumberOfVariables())
                ],
            ),
        )
        model = direct_model(mock)
        @variable(model, S, Ellipsoid(symmetric = true, psd = false))
        @variable(model, T, Ellipsoid(symmetric = true, psd = false))
        expression = PolarLpSum(BALL, PolarLpSum(S, T, 2), 2)
        d = SetProg.data(model)
        d.spaces = SetProg.Spaces()
        d.space = space
        index = SetProg.create_spaces(expression, d.spaces)
        @test SetProg.space_dimension(d.spaces, index) == 2
        @test SetProg.space_dimension(d.spaces, S.space_index) == 2
        @test SetProg.space_dimension(d.spaces, T.space_index) == 2
        @test SetProg.need_variablify(expression)
        @test !SetProg.needs_polyvars(expression)
        @test SetProg.needs_polyvars(Ellipsoid(symmetric = false))
        @test !SetProg.needs_polyvars(Polytope())
        @test !SetProg.needs_polyvars(Domain(Polytope(), MOI.Nonnegatives))
        for variable in (S, T)
            SetProg.load(model, variable)
        end
        resolved = SetProg.variablify(expression)
        @test !SetProg.need_variablify(resolved)
        optimize!(model; ignore_optimize_hook = true)
        gauge_squared = space == SetProg.PrimalSpace ? 5.0 : 5 / 6
        @test SetProg.sublevel_eval(value(expression), [1.0, 1.0]) ≈ 2 + 2gauge_squared
        support_squared = space == SetProg.PrimalSpace ? 5 / 3 : 10.0
        @test support_function([1.0, 1.0], value(LpSum(S, T, 2)))^2 ≈ support_squared
        matrix = space == SetProg.PrimalSpace ? S.variable.Q : polar(S.variable).Q
        mapped = Sets.LinearPreImage(BALL, matrix)
        @test SetProg.sublevel_eval(PolarLpSum(value(mapped), BALL, 2), [1.0, 1.0]) ≈ 15
        @test value(LpSum(S, T, 2)) isa LpSum
        @test SetProg.variablify(LpSum(S, T, 2)) isa LpSum
        outer_polar = polar(LpSum(S, PolarLpSum(BALL, T, 2), 2))
        @test SetProg.need_variablify(outer_polar)
        @test SetProg.variablify(outer_polar) isa PolarLpSum
        SetProg.clear_spaces(outer_polar)
        @test S.space_index === nothing
        @test T.space_index === nothing
        d.spaces = SetProg.Spaces()
        index = SetProg.create_spaces(outer_polar, d.spaces)
        @test SetProg.space_dimension(d.spaces, index) == 2
    end
end

@testset "Quadratic Lₚ inclusion preserves the Schur complement" begin
    model = Model()
    @variable(model, S[1:2, 1:2], Symmetric)
    @variable(model, T[1:2, 1:2], Symmetric)
    @variable(model, Y[1:1, 1:2])
    source = Sets.Ellipsoid(S)
    destination = polar(Sets.Ellipsoid(T))
    A = [0.5 0.25; -0.5 0.75]
    B = [1.0; 0.5;;]
    Q = Symmetric([2.0 0.25; 0.25 1.0])
    R = Symmetric([3.0;;])
    next = Sets.LinearPreImage(destination, A * S + B * Y)
    state = Sets.LinearPreImage(Sets.Ellipsoid(Q), S)
    input = Sets.LinearPreImage(Sets.Ellipsoid(R), Y)
    target = PolarLpSum(next, PolarLpSum(state, input, 2), 2)
    margin = 1e-4
    c = @constraint(model, source ⊆ target, psd_margin = margin)
    expected =
        [
            S S*A'+(B*Y)' S Y';
            A*S+B*Y T zeros(2, 2) zeros(2, 1);
            S zeros(2, 2) inv(Matrix(Q)) zeros(2, 1);
            Y zeros(1, 2) zeros(1, 2) inv(Matrix(R))
        ] - margin * I
    canonical = vec(convert.(AffExpr, expected))
    @test all(isequal_canonical.(constraint_object(c).func, canonical))
    @test num_variables(model) == 8
    @test num_constraints(model; count_variable_in_set_constraints = true) == 1
    @test constraint_object(c).set == MOI.PositiveSemidefiniteConeSquare(7)
    @test Sets.dimension(state) == 2
    @test Sets.space_variables(state) === nothing
    @test Sets.perspective_variable(state) === nothing
    @test_throws ArgumentError @constraint(model, source ⊆ target, psd_margin = -1)
    @test_throws MethodError @constraint(model, source ⊆ PolarLpSum(E, F, 4))
    invalid = PolarLpSum(Sets.LinearPreImage(E, ones(1, 2)), F, 2)
    @test_throws DimensionMismatch @constraint(model, source ⊆ invalid)

    # Direct operands use identity maps; 2.0 and 2 have the same certificate.
    fixed = PolarLpSum(E, polar(F), 2.0)
    c = @constraint(model, Sets.Ellipsoid(Symmetric(10.0 * Matrix(I, 2, 2))) ⊆ fixed)
    @test constraint_object(c).set == MOI.PositiveSemidefiniteConeSquare(6)
    @test num_variables(model) == 8
    mapped = Sets.LinearPreImage(E, [1.0 0.0; 0.0 2.0])
    @test SetProg.sublevel_eval(PolarLpSum(mapped, F, 2), [1.0, 2.0]) ≈
          dot([1.0, 4.0], E.Q * [1.0, 4.0]) + dot([1.0, 2.0], F.Q * [1.0, 2.0])
    @test value(mapped).A == mapped.A
end

# Polynomial operands use the SOS extension; the tests above run without it.
@testset "Polyhedral ray certificates" begin
    for p in (1, 2)
        model = Model(GLPK.Optimizer)
        # GLPK handles the LP; for p=2 inspect the generated SOC model.
        p == 2 && (model = Model())
        S = @variable(model, variable_type = Polytope(symmetric = true, piecewise = BOX))
        cost = linear_pieces(DIAMOND)
        target = PolarLpSum(Sets.LinearPreImage(S, 0.5I(2)), cost, p)
        edge = @constraint(model, S ⊆ target)
        @objective(model, Max, L1_heuristic(volume(S)))
        if p == 1
            set_silent(model)
            optimize!(model)
            @test is_solved_and_feasible(model)
            value_set = value(S)
            @test minimum(eigvals(Symmetric(SetProg.certificate_matrix(edge)))) >= -1e-7
            for x in ([1.0, 0.3], [-0.2, 1.0], [-1.0, -0.8])
                @test SetProg.sublevel_eval(value_set, x) + 1e-7 >=
                      SetProg.sublevel_eval(value(target), x)
            end
        else
            d = SetProg.data(model)
            SetProg.set_space(d, model)
            SetProg.synchronize_perspective(d)
            SetProg.clear_spaces(d)
            SetProg.create_spaces(d)
            SetProg.load(model, d)
            @test num_constraints(model, Vector{AffExpr}, MOI.SecondOrderCone) > 0
        end
    end
end

@testset "Ray certificate margins and dispatch" begin
    source = linear_pieces(BOX)
    cost = linear_pieces(BOX)
    for p in (1, 2)
        model = Model()
        constraints = @constraint(model, source ⊆ PolarLpSum(cost, cost, p), margin = 0.1)
        @test !isempty(constraints)
        @test_throws ArgumentError @constraint(
            model,
            source ⊆ PolarLpSum(cost, cost, p),
            margin = -1
        )
    end
    model = Model()
    @test_throws MethodError @constraint(model, source ⊆ PolarLpSum(cost, cost, 1.5))
    # Several variable terms use independent nonnegative gauge epigraphs.
    a = @variable(model, [1:2])
    terms = SetProg._ray_lp_terms(
        PolarLpSum(Sets.PolarPoint(a), Sets.PolarPoint(-a), 1),
        [1.0, 0.0],
        1,
    )
    constraints = SetProg._add_ray_lp_constraint(model, Val(1), 2.0, terms, 0, [1.0, 0.0])
    @test length(constraints) == 5
end

using SumOfSquares
@testset "Defining polynomial degree and operation exponent are independent" begin
    x = SetProg.new_polyvars(2)
    MP = SumOfSquares.MP
    Q = Matrix(Diagonal([1.0, 2.0, 1.0]))
    quartic = Sets.PolySet(4, SumOfSquares.GramMatrix(Q, MP.monomials(x, 2)))
    quadratic = Sets.PolySet(2, SumOfSquares.GramMatrix(Matrix{Float64}(I, 2, 2), x))
    convex = Sets.ConvexPolySet(4, quartic.p, nothing)
    @test value(convex) === convex
    @test SetProg.sublevel_eval(PolarLpSum(BALL, convex, 2), [1.0, 2.0]) ≈ 10
    for point in ([1.0, 2.0], [-0.5, 1.0])
        @test SetProg.sublevel_eval(PolarLpSum(quartic, quartic, 4), point) ≈
              2dot(point, point)^2
        @test SetProg.sublevel_eval(PolarLpSum(BALL, quartic, 2), point) ≈
              2dot(point, point)
        @test Sets.space_variables(PolarLpSum(BALL, quartic, 2)) == x
        @test SetProg.sublevel_eval(value(PolarLpSum(quartic, quartic, 4)), point) ≈
              2dot(point, point)^2
        @test support_function(point, LpSum(polar(quartic), polar(quadratic), 2))^2 ≈
              2dot(point, point)
    end
    @test SetProg.needs_polyvars(PolarLpSum(BALL, quartic, 2))
    model = Model()
    @variable(model, S, PolySet(symmetric = true, degree = 4))
    @test SetProg.needs_polyvars(PolarLpSum(quartic, S, 4))
    # The sum preserves the coordinate chart of perspective representations.
    householder(z) = SetProg.perspective_dual_polyset(
        2,
        quadratic.p,
        SetProg.InteriorPoint(zeros(2)),
        z,
        x,
    ).set
    z = SetProg.new_perspective_variable()
    H = householder(z)
    @test Sets.perspective_variable(PolarLpSum(H, BALL, 2)) == z
    @test Sets.perspective_variable(PolarLpSum(BALL, H, 2)) == z
    @test Sets.perspective_variable(LpSum(H, H, 2)) == z
    @test_throws ArgumentError Sets.perspective_variable(
        LpSum(H, householder(SetProg.new_perspective_variable()), 2),
    )
    @test_throws ArgumentError PolarLpSum(
        quartic,
        Sets.PolySet(
            4,
            SumOfSquares.GramMatrix(Q, MP.monomials(SetProg.new_polyvars(2), 2)),
        ),
        4,
    )
    nonhomogeneous = Sets.PolySet(
        0:4,
        SumOfSquares.GramMatrix(Matrix{Float64}(I, 6, 6), MP.monomials(x, 0:2)),
    )
    @test_throws MethodError SetProg.sublevel_eval(
        PolarLpSum(BALL, nonhomogeneous, 2),
        ones(2),
    )
end

@testset "Polynomial Lₚ inclusion certificates" begin
    x = SetProg.new_polyvars(2)
    MP = SumOfSquares.MP
    gram = SumOfSquares.GramMatrix(Matrix(Diagonal([2.0, 4.0, 2.0])), MP.monomials(x, 2))
    source = Sets.PolySet(4, gram)
    A = [0.5 0.0; 0.0 0.25]
    K = [0.25 0.5]
    input = Sets.Ellipsoid(Symmetric([2.0;;]))
    target = PolarLpSum(
        Sets.LinearPreImage(source, A),
        PolarLpSum(BALL, Sets.LinearPreImage(input, K), 4),
        4,
    )
    ext = Base.get_extension(SetProg, :SetProgSumOfSquaresExt)
    expected =
        MP.polynomial(source.p) - MP.polynomial(source.p)(x => A*x) - sum(x .^ 2)^2 -
        (2only(K*x)^2)^2
    @test all(
        iszero,
        MP.coefficients(
            ext._powered_gauge_polynomial(source, x, 4) - MP.polynomial(source.p),
        ),
    )
    @test ext._powered_gauge_polynomial(target, x, 4) ≈ MP.polynomial(source.p) - expected
    actual = JuMP.build_constraint(error, source, SetProg.PowerSet(target); margin = 0.1)
    reference = JuMP.build_constraint(error, expected - 0.1sum(x .^ 2)^2, SOSCone())
    model = Model()
    actual_ref = JuMP.add_constraint(model, actual)
    reference_ref = JuMP.add_constraint(model, reference)
    @test JuMP.constant.(constraint_object(actual_ref).func) ≈
          JuMP.constant.(constraint_object(reference_ref).func)
    @test_throws ArgumentError JuMP.build_constraint(
        error,
        source,
        SetProg.PowerSet(target);
        margin = -1,
    )
    @test_throws MethodError JuMP.build_constraint(
        error,
        source,
        SetProg.PowerSet(PolarLpSum(BALL, source, 2)),
    )
    mixed = PolarLpSum(BALL, PolarLpSum(BALL, BALL, 2), 4)
    @test_throws MethodError JuMP.build_constraint(error, source, SetProg.PowerSet(mixed))
    convex = Sets.ConvexPolySet(4, gram, nothing)
    @test ext._powered_gauge_polynomial(convex, x, 4) == MP.polynomial(source.p)
    model = Model()
    @variable(model, S, PolySet(symmetric = true, psd = false, degree = 4))
    mapped = Sets.LinearPreImage(S, A)
    expression = PolarLpSum(mapped, BALL, 4)
    c = @constraint(model, source ⊆ expression)
    d = SetProg.data(model)
    SetProg.set_space(d, model)
    @test d.space == SetProg.PrimalSpace
    SetProg.clear_spaces(d)
    SetProg.create_spaces(d)
    @test SetProg.space_dimension(d.spaces, S.space_index) == 2
    @test SetProg.needs_polyvars(mapped)
    @test !SetProg.needs_polyvars(Sets.LinearPreImage(BALL, A))
    SetProg.load(model, S)
    @test !SetProg.need_variablify(SetProg.variablify(expression))
    @test SetProg.variablify(constraint_object(c)) isa JuMP.AbstractConstraint
    SetProg.clear_spaces(expression)
    @test S.space_index === nothing
end

@testset "Piecewise quadratic Bellman inclusions" begin
    model = Model(() -> Hypatia.Optimizer(verbose = false))
    @variable(model, S, Ellipsoid(symmetric = true, piecewise = BOX))
    @variable(model, T, Ellipsoid(symmetric = true, piecewise = DIAMOND))
    A = [0.3 0.1; -0.05 0.2]
    # Different partitions, nested sums and a rectangular input map all
    # contribute to the common refinement in the source coordinates.
    input = Sets.LinearPreImage(Sets.Ellipsoid(Symmetric([1.0;;])), [0.1 -0.2])
    stage = PolarLpSum(BALL, input, 2)
    c = @constraint(
        model,
        S ⊆ PolarLpSum(Sets.LinearPreImage(T, A), stage, 2),
        margin = 0.01
    )
    d = @constraint(model, T ⊆ PolarLpSum(Sets.LinearPreImage(S, A), stage, 2))
    @objective(model, Max, L1_heuristic(volume(S)) + L1_heuristic(volume(T)))
    optimize!(model)
    @test is_solved_and_feasible(model; allow_almost = true)
    @test SetProg.data(model).space == SetProg.PrimalSpace
    fitted = value.([S, T])
    for cref in (c, d)
        @test eigmin(Symmetric(SetProg.certificate_matrix(cref))) >= -1e-6
    end
    for t in range(0, 2pi; length = 101)
        x = [cos(t), sin(t)]
        for (i, j, margin) in ((1, 2, 0.01), (2, 1, 0.0))
            @test SetProg.sublevel_eval(fitted[i], x) + 1e-6 >=
                  SetProg.sublevel_eval(fitted[j], A*x) +
                  dot(x, x) +
                  only([0.1 -0.2]*x)^2 +
                  margin * dot(x, x)
        end
    end
    @test_throws ArgumentError JuMP.add_constraint(
        model,
        SetProg.InclusionConstraint(fitted[1], PolarLpSum(BALL, BALL, 2), (; margin = -1)),
    )
    @test_throws MethodError SetProg._quadratic_lp_regions(PolarLpSum(BALL, BALL, 1))
    # Each residual is indefinite globally, but nonnegative on its cone:
    # 0.5(x1^2-x2^2) on |x1| >= |x2|, and the reverse on the other cones.
    source = Sets.Piecewise(
        [
            Sets.Ellipsoid(
                Symmetric(iszero(h.a[1]) ? [1.0 0.0; 0.0 2.0] : [2.0 0.0; 0.0 1.0]),
            ) for h in halfspaces(BOX)
        ],
        BOX,
    )
    stage = PolarLpSum(
        Sets.Ellipsoid(Symmetric(0.75Matrix{Float64}(I, 2, 2))),
        Sets.Ellipsoid(Symmetric(0.75Matrix{Float64}(I, 2, 2))),
        2,
    )
    check = Model(() -> Hypatia.Optimizer(verbose = false))
    ref = @constraint(check, source ⊆ stage)
    optimize!(check)
    @test is_solved_and_feasible(check)
    @test eigmin(Symmetric(SetProg.certificate_matrix(ref))) >= -1e-6
    check = Model(() -> Hypatia.Optimizer(verbose = false))
    @constraint(check, Sets.Piecewise(fill(BALL, 4), BOX) ⊆ stage)
    optimize!(check)
    @test termination_status(check) == MOI.INFEASIBLE
    # Preimage cones must remain valid for singular maps, including zero.
    for M in ([0.3 0.0; 0.0 0.0], zeros(2, 2))
        check = Model(() -> Hypatia.Optimizer(verbose = false))
        ref = @constraint(
            check,
            fitted[1] ⊆ PolarLpSum(Sets.LinearPreImage(fitted[2], M), BALL, 2)
        )
        optimize!(check)
        @test is_solved_and_feasible(check)
        @test eigmin(Symmetric(SetProg.certificate_matrix(ref))) >= -1e-6
    end
end

end
