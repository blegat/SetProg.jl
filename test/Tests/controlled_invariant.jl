using LinearAlgebra
using Test

using SetProg, SetProg.Sets
using Polyhedra
using MultivariatePolynomials

using JuMP
const MOIT = MOI.Test

function ci_square_test(optimizer, config::MOIT.TestConfig,
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

    Δt = 0.5
    A = [1.0 Δt]
    E = [1.0 0.0]

    if variable.symmetric
        @constraint(model, A * ◯ ⊆ E * ◯)
    else
        @constraint(model, A * ◯ ⊆ E * ◯, S_procedure_scaling = 1.0)
    end

    @objective(model, inner ? MOI.MAX_SENSE : MOI.MIN_SENSE,
               metric(volume(◯)))

    SetProg.optimize!(model)
    @test JuMP.termination_status(model) == MOI.OPTIMAL
    @test JuMP.objective_sense(model) == MOI.MAX_SENSE
    @test JuMP.objective_value(model) ≈ objective_value atol=config.atol rtol=config.rtol
    set_test(JuMP.value(◯))
end

function ci_ell_homogeneous_test(optimizer, config)
    ci_square_test(optimizer, config, true,
                   Ellipsoid(symmetric=true, dimension=2), nth_root, √15/4,
                   ◯ -> begin
                       @test ◯ isa Sets.Polar{Float64, Sets.EllipsoidAtOrigin{Float64}}
                       @test Sets.polar(◯).Q ≈ Symmetric([1.0 -1/4; -1/4 1.0]) atol=config.atol rtol=config.rtol
                   end)
end

function ci_ell_nonhomogeneous_test(optimizer, config)
    ci_square_test(optimizer, config, true,
                   Ellipsoid(point=SetProg.InteriorPoint([0.0, 0.0])),
                   nth_root, √15/4,
                   ◯ -> begin
                       @test ◯ isa Sets.PerspectiveDual{Float64, Sets.Householder{Float64, Sets.ShiftedEllipsoid{Float64}, Float64}}
                       z = Sets.perspective_variable(◯)
                       x, y = Sets.space_variables(◯)
                       ◯_dual = Sets.perspective_dual(◯)
                       @test ◯_dual.p ≈ -z^2 + x^2 - x*y/2 + y^2 atol=config.atol rtol=config.rtol
                       @test ◯_dual.set.Q ≈ Symmetric([1.0 -1/4; -1/4 1.0]) atol=config.atol rtol=config.rtol
                       @test ◯_dual.set.b ≈ [0.0, 0.0] atol=config.atol rtol=config.rtol
                       @test ◯_dual.set.β ≈ -1.0 atol=config.atol rtol=config.rtol
                       @test Sets._householder(◯_dual.h) ≈ [-1.0 0.0 0.0
                                                             0.0 1.0 0.0
                                                             0.0 0.0 1.0] atol=config.atol rtol=config.rtol
                   end)
end

function ci_quad_nonhomogeneous_test(optimizer, config)
    ci_square_test(optimizer, config, true,
                   PolySet(degree=2, convex=true, point=SetProg.InteriorPoint([0.0, 0.0])),
                   set -> L1_heuristic(set, [1.0, 1.0]), 8/3,
                   ◯ -> begin
                       @test ◯ isa Sets.PerspectiveDual{Float64, Sets.Householder{Float64, Sets.ConvexPolynomialSet{Float64}, Float64}}
                       z = Sets.perspective_variable(◯)
                       x, y = Sets.space_variables(◯)
                       ◯_dual = Sets.perspective_dual(◯)
                       # The coefficient of `x*y` does not influence the volume
                       # and with the values of the other parameters, it should
                       # simply be in the interval [-2, -0.5].
                       α = coefficient(◯_dual.p, x*y)
                       @test α ≥ -2 - 2config.atol - config.rtol
                       @test α ≤ -0.5 + 0.5config.atol + config.rtol
                       @test ◯_dual.p ≈ -z^2 + x^2 + α*x*y + y^2 atol=config.atol rtol=config.rtol
                   end)
end

function ci_quartic_homogeneous_test(optimizer, config)
    ci_square_test(optimizer, config, true,
                   PolySet(symmetric=true, degree=4, convex=true),
                   set -> L1_heuristic(set, [1.0, 1.0]),
                   64/15,
                   ◯ -> begin
                       @test ◯ isa Sets.Polar{Float64, Sets.ConvexPolynomialSublevelSetAtOrigin{Float64}}
                       @test Sets.polar(◯).degree == 4
                       x, y = variables(Sets.polar(◯).p)
                       α = coefficient(Sets.polar(◯).p, x^3*y) / 2
                       q = x^4 + 2α*x^3*y + 6x^2*y^2 + 2α*x*y^3 + y^4
                       @test all(eigvals(Matrix(Sets.polar(◯).p.Q)) .≥ -config.atol)
                       @test polynomial(Sets.polar(◯).p) ≈ q atol=config.atol rtol=config.rtol
                       convexity_proof = Sets.convexity_proof(◯)
                       @test convexity_proof.n == 4

                       hess = 6 * [2.0, α, 2.0, α, 2.0,
                                   2.0, 2.0, α, α, 2.0]
                       Hess = SetProg.SumOfSquares.MultivariateMoments.SymMatrix(hess, 4)
                       @test all(eigvals(Matrix(Hess)) .≥ -config.atol)
                       @test convexity_proof.Q ≈ hess atol=config.atol rtol=config.rtol
                   end)
end

const ci_tests = Dict("ci_ell_homogeneous" => ci_ell_homogeneous_test,
                      "ci_ell_nonhomogeneous" => ci_ell_nonhomogeneous_test,
                      "ci_quad_nonhomogeneous" => ci_quad_nonhomogeneous_test,
                      "ci_quartic_homogeneous" => ci_quartic_homogeneous_test)

@test_suite ci
