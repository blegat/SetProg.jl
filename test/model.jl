using LinearAlgebra
using Test

using SetProg, SetProg.Sets
using Polyhedra
using DynamicPolynomials
using JuMP
import GLPK

include("utilities.jl")

const _lp_solver = optimizer_with_attributes(
    GLPK.Optimizer,
    MOI.Silent() => true,
    "presolve" => GLPK.GLP_ON,
)
const _lib = Polyhedra.DefaultLibrary{Float64}(_lp_solver)

function _square(n)
    interval = HalfSpace([1.0], 1.0) ∩ HalfSpace([-1.0], 1.0)
    h = interval
    for _ in 2:n
        h = h * interval
    end
    return polyhedron(h, _lib)
end

# Model whose optimizer does not return any solution, only useful to check
# how the set program is reformulated
_mock_model() = direct_model(bridged_mock(mock -> nothing))

function _num_constraints(model, F, S)
    return MOI.get(backend(model), MOI.NumberOfConstraints{F,S}())
end

@testset "Polytope" begin
    @test_throws DimensionMismatch Polytope(dimension = 3, piecewise = _square(2))
    @test Polytope(symmetric = true, piecewise = _square(2)).dimension == 2
    @testset "Non-symmetric" begin
        model = _mock_model()
        @variable(model, S, Polytope(dimension = 2))
        @constraint(model, S ⊆ _square(2))
        err = ErrorException("Non-symmetric polytope not supported yet")
        @test_throws err SetProg.optimize!(model)
    end
    @testset "Discrete-time controlled invariant" begin
        # See `docs/src/examples/discrete_controlled.jl`
        A = [1.0 0.5]
        E = [1.0 0.0]
        □ = _square(2)
        ◇ = polyhedron(convexhull([1.0, 0], [0, 1], [-1, 0], [0, -1]), _lib)
        for (piecewise, obj) in [(□, 8 / 3), (◇, 4√2 / 3)]
            model = Model(_lp_solver)
            @variable(model, S, Polytope(symmetric = true, piecewise = piecewise))
            @constraint(model, S ⊆ □)
            @constraint(model, A * S ⊆ E * S)
            @objective(model, Max, L1_heuristic(volume(S)))
            optimize!(model)
            @test termination_status(model) == MOI.OPTIMAL
            @test objective_value(model) ≈ obj
            sol = value(S)
            @test sol isa Sets.Polar{Float64,<:Sets.Piecewise{Float64,Sets.PolarPoint{Float64}}}
            @test Sets.dimension(sol) == 2
            @test length(Sets.polar(sol).sets) == 4
        end
    end
end

@testset "Continuous-time controlled invariant" begin
    # See `docs/src/examples/continuous_controlled.jl`
    A = [
        0.0 1.0 0.0
        0.0 0.0 1.0
        0.0 0.0 0.0
    ]
    E = [
        1.0 0.0 0.0
        0.0 1.0 0.0
    ]
    C = A[1:2, :]
    □_3 = _square(3)
    dirs = [[-1 + √3, -1 + √3], [-1, 1]]
    function build(family)
        model = _mock_model()
        @variable(model, S, family)
        @constraint(model, S ⊆ □_3)
        x = boundary_point(S, :x)
        @test tangent_cone(S, x) isa SetProg.TangentCone
        @constraint(model, C * x in E * tangent_cone(S, x))
        S_2 = project(S, 1:2)
        @variable(model, γ)
        for point in dirs
            @constraint(model, γ * point in S_2)
        end
        @objective(model, Max, γ)
        SetProg.optimize!(model)
        return model
    end
    @testset "Ellipsoid" begin
        model = build(Ellipsoid(symmetric = true))
        # One per facet of the cube and one per direction
        @test _num_constraints(model, MOI.ScalarAffineFunction{Float64}, MOI.LessThan{Float64}) == 6
        # For the tangent cone constraint
        @test _num_constraints(model, MOI.VectorAffineFunction{Float64}, MOI.PositiveSemidefiniteConeTriangle) == 2
    end
    @testset "Piecewise ellipsoid" begin
        model = build(Ellipsoid(symmetric = true, piecewise = polar(□_3)))
        # One PSD variable per piece
        @test _num_constraints(model, MOI.VectorOfVariables, MOI.PositiveSemidefiniteConeTriangle) == 8
    end
    @testset "PolySet" begin
        model = build(PolySet(symmetric = true, degree = 4, convex = true))
        @test _num_constraints(model, MOI.ScalarAffineFunction{Float64}, MOI.LessThan{Float64}) == 6
    end
end

@testset "Membership" begin
    # Vertex representation with a line and a ray
    P = polyhedron(vrep([[1.0, 0.0]], [Line([0.0, 1.0])], [Ray([1.0, 1.0])]))
    @testset "Line and ray in primal $name" for (name, family, metric) in [
        ("Ellipsoid", Ellipsoid(symmetric = true, dimension = 2), nth_root),
        (
            "PolySet",
            PolySet(symmetric = true, degree = 2, dimension = 2, convex = true),
            vol -> L1_heuristic(vol, [1.0, 1.0]),
        ),
    ]
        model = _mock_model()
        @variable(model, S, family)
        @constraint(model, P ⊆ S)
        @objective(model, Min, metric(volume(S)))
        SetProg.optimize!(model)
        # The line and the ray must be in the recession cone
        @test _num_constraints(model, MOI.ScalarAffineFunction{Float64}, MOI.EqualTo{Float64}) == 2
        # The point must be in the set
        @test _num_constraints(model, MOI.ScalarAffineFunction{Float64}, MOI.LessThan{Float64}) == 1
    end
    @testset "Point in dual $name" for (name, family, metric) in [
        ("Ellipsoid", Ellipsoid(symmetric = true, dimension = 2), nth_root),
        (
            "PolySet",
            PolySet(symmetric = true, degree = 4, dimension = 2, convex = true),
            vol -> L1_heuristic(vol, [1.0, 1.0]),
        ),
    ]
        model = _mock_model()
        @variable(model, S, family)
        @constraint(model, S ⊆ _square(2))
        @constraint(model, [0.5, 0.5] in S)
        @objective(model, Max, metric(volume(S)))
        SetProg.optimize!(model)
        # One per facet of the square and one for the point
        @test _num_constraints(model, MOI.ScalarAffineFunction{Float64}, MOI.LessThan{Float64}) == 4
    end
    @testset "Point in dual piecewise ellipsoid" begin
        model = _mock_model()
        □ = _square(2)
        @variable(model, S, Ellipsoid(symmetric = true, piecewise = polar(□)))
        @constraint(model, S ⊆ □)
        @constraint(model, [0.5, 0.5] in S)
        @objective(model, Max, L1_heuristic(volume(S)))
        SetProg.optimize!(model)
        # One PSD variable per piece
        @test _num_constraints(model, MOI.VectorOfVariables, SetProg.SumOfSquares.PositiveSemidefinite2x2ConeTriangle) == 4
    end
end

@testset "Affine objective" begin
    model = _mock_model()
    @variable(model, S, Ellipsoid(symmetric = true, dimension = 2))
    @variable(model, T, Ellipsoid(symmetric = true, dimension = 2))
    @constraint(model, S ⊆ _square(2))
    @constraint(model, T ⊆ _square(2))
    l1(set) = L1_heuristic(volume(set), [1.0, 1.0])
    obj = l1(S) + l1(T)
    @test obj isa SetProg.AffineExpression
    @test length((obj + l1(S)).terms) == 3
    obj = l1(T) + obj
    @test length(obj.terms) == 3
    @objective(model, Max, obj)
    SetProg.optimize!(model)
    @test JuMP.objective_sense(model) == MOI.MAX_SENSE
    @test JuMP.objective_function(model) isa JuMP.AffExpr
end

@testset "Printing" begin
    model = _mock_model()
    @variable(model, S, Ellipsoid(symmetric = true, dimension = 2))
    □ = _square(2)
    c = @constraint(model, S ⊆ □, base_name = "c")
    @test JuMP.name(c) == "c"
    @test JuMP.is_valid(model, c)
    @test !JuMP.is_valid(_mock_model(), c)
    @test JuMP.constraint_object(c) isa SetProg.InclusionConstraint
    @test sprint(print, S) == "S"
    @test startswith(sprint(print, c), "c : S ⊆ HalfSpace")
    @test startswith(
        JuMP.constraint_string(MIME"text/latex"(), JuMP.constraint_object(c)),
        "S \\subseteq HalfSpace",
    )
    m = @constraint(model, [0.5, 0.5] in S)
    @test sprint(print, m) == "[0.5, 0.5] ∈ S"
    @test JuMP.constraint_string(MIME"text/latex"(), JuMP.constraint_object(m)) ==
          "[0.5, 0.5] \\in S"
    @test sprint(show, L1_heuristic(volume(S))) == "L1-heuristic(S)"
    @test sprint(show, nth_root(volume(S))) == "volume^(1/n)(S)"
end

@testset "Incompatible spaces" begin
    model = _mock_model()
    @variable(model, S, Ellipsoid(symmetric = true, dimension = 2))
    □ = _square(2)
    # Requires the dual space
    @constraint(model, S ⊆ □)
    # Requires the primal space
    @constraint(model, □ ⊆ S)
    err = ErrorException(
        "Incompatible constraints/objective, some require to do the modeling in the primal space and some in the dual space.",
    )
    @test_throws err SetProg.optimize!(model)
end

@testset "Macros" begin
    model = _mock_model()
    @variable(model, S, Ellipsoid(symmetric = true, dimension = 2))
    □ = _square(2)
    @test_throws ErrorException @macroexpand @constraint(model, S ⊂ □)
    @test_throws ErrorException @macroexpand @constraint(model, □ ⊃ S)
    c = @constraint(model, □ ⊇ S)
    @test JuMP.constraint_object(c).subset === S
    cs = @constraint(model, [S, S] .⊆ [□, □])
    @test length(cs) == 2
    @test all(c -> JuMP.constraint_object(c).subset === S, cs)
end

@testset "Log volume of PolySet" begin
    model = _mock_model()
    @variable(model, S, PolySet(symmetric = true, degree = 4, dimension = 2, convex = true))
    @constraint(model, S ⊆ _square(2))
    @objective(model, Max, log(volume(S)))
    SetProg.optimize!(model)
    @test JuMP.objective_sense(model) == MOI.MAX_SENSE
    @test _num_constraints(model, MOI.VectorAffineFunction{Float64}, MOI.LogDetConeTriangle) == 1
    @test _num_constraints(model, MOI.VectorOfVariables, MOI.RootDetConeTriangle) == 0
end

@testset "Constant superset" begin
    □ = _square(2)
    @polyvar x[1:2]
    unit_ellipsoid = Sets.Ellipsoid(Symmetric(Matrix(1.0I, 2, 2)))
    # (x₁² + x₂²)² = x₁⁴ + 2x₁²x₂² + x₂⁴
    unit_quartic = Sets.PolySet(
        4,
        SetProg.GramMatrix(Matrix(Diagonal([1.0, 2.0, 1.0])), monomials(x, 2)),
    )
    unit_piecewise = Sets.Piecewise([unit_ellipsoid for _ in 1:4], □)
    unit_convex_quartic = Sets.ConvexPolySet(4, unit_quartic.p, nothing)
    # The gauge of the square, linear on the cone over each of its facets
    unit_polytope = Sets.Piecewise([Sets.PolarPoint(h.a / h.β) for h in halfspaces(□)], □)
    @testset "$name" for (name, family, unit, F, S) in [
        (
            "Ellipsoid",
            Ellipsoid(symmetric = true, dimension = 2),
            unit_ellipsoid,
            MOI.VectorAffineFunction{Float64},
            SetProg.SumOfSquares.PositiveSemidefinite2x2ConeTriangle,
        ),
        (
            "PolySet",
            PolySet(symmetric = true, degree = 4, variables = x),
            unit_quartic,
            MOI.VectorAffineFunction{Float64},
            SetProg.SumOfSquares.SOSPolynomialSet,
        ),
        (
            "Convex PolySet",
            PolySet(symmetric = true, degree = 4, convex = true, variables = x),
            unit_convex_quartic,
            MOI.VectorAffineFunction{Float64},
            SetProg.SumOfSquares.SOSPolynomialSet,
        ),
        (
            "Piecewise",
            Ellipsoid(symmetric = true, piecewise = □),
            unit_piecewise,
            MOI.VectorAffineFunction{Float64},
            SetProg.SumOfSquares.PositiveSemidefinite2x2ConeTriangle,
        ),
        (
            "Polytope",
            Polytope(symmetric = true, piecewise = □),
            unit_polytope,
            MOI.VectorAffineFunction{Float64},
            MOI.Zeros,
        ),
    ]
        model = _mock_model()
        @variable(model, V, family)
        c = @constraint(model, V ⊆ unit)
        @constraint(model, [1.0 0.5; 0.0 0.5] * V ⊆ V)
        SetProg.optimize!(model)
        # The gauge function of `V` is compared with the one of `unit`
        @test SetProg.data(model).space == SetProg.PrimalSpace
        types = MOI.get(backend(model), MOI.ListOfConstraintTypesPresent())
        @test any(((f, s),) -> f == F && s <: S, types)
    end
    @testset "Ellipsoid in dual space" begin
        model = _mock_model()
        @variable(model, V, Ellipsoid(symmetric = true, dimension = 2))
        # Forces the dual space
        @constraint(model, V ⊆ □)
        @constraint(model, V ⊆ Sets.Ellipsoid(Symmetric(2 * unit_ellipsoid.Q)))
        SetProg.optimize!(model)
        @test SetProg.data(model).space == SetProg.DualSpace
        # The polar of `V`, of matrix `Q`, must contain the polar of the
        # superset, of matrix `inv(2I) = I / 2`, so `I / 2 - Q` is PSD
        F = MOI.VectorAffineFunction{Float64}
        S = SetProg.SumOfSquares.PositiveSemidefinite2x2ConeTriangle
        @test _num_constraints(model, F, S) == 1
    end
end
