using LinearAlgebra
using Test

using SetProg, SetProg.Sets
using Polyhedra
using MultivariatePolynomials

using JuMP
const MOIT = MOI.Test

function invariant_square_test(optimizer, config::MOIT.TestConfig,
                               inner::Bool, variable::SetProg.AbstractVariable,
                               metric::Function, objective_value, set_test)
    model = _model(optimizer)

    □ = polyhedron(HalfSpace([1, 0], 1.0) ∩ HalfSpace([-1, 0], 1) ∩ HalfSpace([0, 1], 1) ∩ HalfSpace([0, -1], 1))

    @variable(model, ◯, variable)
    if inner
        cref = @constraint(model, ◯ ⊆ □)
    else
        cref = @constraint(model, □ ⊆ ◯)
    end

    A = [0.0 -1.0
         1.0  0.0]

    if variable.symmetric
        @constraint(model, A * ◯ ⊆ ◯)
    else
        @constraint(model, A * ◯ ⊆ ◯, S_procedure_scaling = 1.0)
    end

    @objective(model, inner ? MOI.MAX_SENSE : MOI.MIN_SENSE,
               metric(volume(◯)))

    SetProg.optimize!(model)
    @test JuMP.termination_status(model) == MOI.OPTIMAL
    @test JuMP.objective_sense(model) == MOI.MAX_SENSE
    @test JuMP.objective_value(model) ≈ objective_value atol=config.atol rtol=config.rtol
    return set_test(JuMP.value(◯))
end

function maximal_invariant_ell_homogeneous_test(optimizer, config)
    invariant_square_test(
        optimizer, config, true, Ellipsoid(symmetric=true, dimension=2),
        nth_root, 1.0,
        ◯ -> begin
            @test ◯ isa Sets.Polar{Float64, Sets.EllipsoidAtOrigin{Float64}}
            @test Sets.polar(◯).Q ≈ Symmetric([1.0 0.0; 0.0 1.0]) atol=config.atol rtol=config.rtol
        end)
end

function maximal_convex_invariant_quad_homogeneous_test(optimizer, config)
    invariant_square_test(
        optimizer, config, true,
        PolySet(degree=2, convex=true, symmetric=true),
        nth_root, 2.0,
        ◯ -> begin
            @test ◯ isa Sets.Polar{Float64, Sets.ConvexPolynomialSublevelSetAtOrigin{Float64}}
            x, y = Sets.space_variables(◯)
            ◯_polar = Sets.polar(◯)
            @test ◯_polar.p ≈ x^2 + y^2 atol=config.atol rtol=config.rtol
        end)
end

function minimal_invariant_ell_homogeneous_test(optimizer, config)
    invariant_square_test(
        optimizer, config, false, Ellipsoid(symmetric=true, dimension=2),
        nth_root, 0.5,
        ◯ -> begin
            @test ◯ isa Sets.EllipsoidAtOrigin{Float64}
            @test ◯.Q ≈ Symmetric([0.5 0.0; 0.0 0.5]) atol=config.atol rtol=config.rtol
        end)
end

function minimal_invariant_quad_homogeneous_test(optimizer, config)
    invariant_square_test(
        optimizer, config, false,
        PolySet(degree=2, symmetric=true),
        set -> L1_heuristic(set, ones(2)), 4/3,
        ◯ -> begin
            @test ◯ isa Sets.PolynomialSublevelSetAtOrigin{Float64}
            x, y = Sets.space_variables(◯)
            @test ◯.p ≈ 0.5x^2 + 0.5y^2 atol=config.atol rtol=config.rtol
        end)
end

function minimal_convex_invariant_quad_homogeneous_test(optimizer, config)
    invariant_square_test(
        optimizer, config, false,
        PolySet(degree=2, convex=true, symmetric=true),
        nth_root, 1.0,
        ◯ -> begin
            @test ◯ isa Sets.ConvexPolynomialSublevelSetAtOrigin{Float64}
            x, y = Sets.space_variables(◯)
            @test ◯.p ≈ 0.5x^2 + 0.5y^2 atol=config.atol rtol=config.rtol
        end)
end

const invariant_tests = Dict("maximal_invariant_ell_homogeneous_test"  => maximal_invariant_ell_homogeneous_test,
                             "maximal_convex_invariant_quad_homogeneous_test" => maximal_convex_invariant_quad_homogeneous_test,
                             "minimal_invariant_ell_homogeneous_test"  => minimal_invariant_ell_homogeneous_test,
                             "minimal_invariant_quad_homogeneous_test" => minimal_invariant_quad_homogeneous_test,
                             "minimal_convex_invariant_quad_homogeneous_test" => minimal_convex_invariant_quad_homogeneous_test)

@test_suite invariant
