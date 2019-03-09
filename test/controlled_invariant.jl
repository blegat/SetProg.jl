using LinearAlgebra
using Test

using SetProg, SetProg.Sets
using Polyhedra
using MultivariatePolynomials

using JuMP
const MOI = JuMP.MOI
const MOIT = MOI.Test

function mock(mock_optimize!::Function)
    mock = MOI.Utilities.MockOptimizer(JuMP._MOIModel{Float64}())
    MOI.Utilities.set_mock_optimize!(mock, mock_optimize!)
    return mock
end

function ci_square_test(optimizer::MOI.AbstractOptimizer, config::MOIT.TestConfig,
                        inner::Bool, variable::SetProg.AbstractVariable,
                        metric::Function, objective_value, set_test)
    MOI.empty!(optimizer)
    model = JuMP.direct_model(optimizer)

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
                       @test ◯_dual.p ≈ -z^2 + x^2 - 0.98661919361x*y + y^2 atol=config.atol rtol=config.rtol
                   end)
end

const ci_quartic_α = -1.4529635030551264
const ci_quartic_β =  0.25613759221607024
const ci_quartic_γ = -8.717781018330758
const ci_quartic_hess = [12.0, ci_quartic_γ, 12.0, ci_quartic_γ, 12.0,
                         12.0, 12.0, ci_quartic_γ, ci_quartic_γ, 12.0]

const ci_quartic_obj = 0.25369938382997853

function ci_quartic_homogeneous_test(optimizer, config)
    α = 2.905927006110253
    ci_square_test(optimizer, config, true,
                   PolySet(symmetric=true, degree=4, convex=true),
                   set -> L1_heuristic(set, [1.0, 1.0]),
                   64/15,
                   ◯ -> begin
                       @test ◯ isa Sets.Polar{Float64, Sets.ConvexPolynomialSublevelSetAtOrigin{Float64}}
                       @test Sets.polar(◯).degree == 4
                       x, y = variables(Sets.polar(◯).p)
                       q = x^4 - α*x^3*y + 6x^2*y^2 - α*x*y^3 + y^4
                       @test polynomial(Sets.polar(◯).p) ≈ q atol=config.atol rtol=config.rtol
                       convexity_proof = Sets.convexity_proof(◯)
                       @test convexity_proof.n == 4
                       @test convexity_proof.Q ≈ ci_quartic_hess atol=config.atol rtol=config.rtol
                   end)
end

@testset "Controlled invariant" begin
    config = MOIT.TestConfig()
    @testset "Ellipsoid" begin
        # Q = [1 0
        #      0 1]
        # t = √det(Q) = 1
        Q = [1.0, -1/4, 1.0]
        t = √15/4
        @testset "Homogeneous" begin
            ci_ell_homogeneous_test(mock(mock -> MOI.Utilities.mock_optimize!(mock, [Q; t])),
                                    config)
        end
        @testset "Non-homogeneous" begin
            @testset "Ellipsoid" begin
                β = -1.0
                b = [0.0, 0.0]
                ci_ell_nonhomogeneous_test(mock(mock -> MOI.Utilities.mock_optimize!(mock, [Q; β; b; zeros(3); t])),
                                           config)
            end
            @testset "PolySet" begin
                β = -1.0
                b = [0.0, 0.0]
                Q = [1.0, -0.4933095968, 1.0]

                ci_quad_nonhomogeneous_test(mock(mock -> begin
                                                     # β+1 b[1] b[2]
                                                     #  .  Q[1] Q[2]
                                                     #  .   .   Q[3]
                                                     MOI.Utilities.mock_optimize!(mock, [β+1; b[1]; Q[1]; b[2]; Q[2]; Q[3]; 2Q; 0.0; 0.0; 0.24331])
                                                 end),
                                            config)
            end
        end
    end
    @testset "Quartic" begin
        # The PSD matrix for the variable  is 3 x 3 so 3 * (3+1) / 2 = 6
        # The PSD matrix for the convexity is 6 x 6 so 6 * (6+1) / 2 = 21
        # entries
        # 1 variable for t
        # hence 28 variables
        sol = [1.0; ci_quartic_α; 6 - 2ci_quartic_β; ci_quartic_β; ci_quartic_α; 1.0;
               ci_quartic_hess; ci_quartic_obj]
        ci_quartic_homogeneous_test(mock(mock -> MOI.Utilities.mock_optimize!(mock, sol)),
                                    config)
    end
end

const citests = Dict("ci_ell_homogeneous" => ci_ell_homogeneous_test,
                     "ci_ell_nonhomogeneous" => ci_ell_nonhomogeneous_test,
                     "ci_quad_nonhomogeneous" => ci_quad_nonhomogeneous_test,
                     "ci_quartic_homogeneous" => ci_quartic_homogeneous_test)

MOIT.@moitestset ci
