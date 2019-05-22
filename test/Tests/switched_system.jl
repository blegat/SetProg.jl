using LinearAlgebra
using Test

using SetProg, SetProg.Sets
using Polyhedra
using MultivariatePolynomials
using DynamicPolynomials

using JuMP
const MOIT = MOI.Test

function switched_system_test(optimizer, config::MOIT.TestConfig,
                              variable::SetProg.AbstractVariable, γ,
                              feasible, objective_value, set_test, dual_test)
    model = _model(optimizer)

    @variable(model, ◯, variable)

    # See https://github.com/blegat/SwitchOnSafety.jl/blob/master/examples/LPJ17e43.ipynb
    A1 = [1  0
          1  0] / γ
    A2 = [0  1
          0 -1] / γ

    cref1 = @constraint(model, A1 * ◯ ⊆ ◯)
    cref2 = @constraint(model, A2 * ◯ ⊆ ◯)

	if objective_value !== nothing
		@objective(model, MOI.MAX_SENSE, L1_heuristic(volume(◯), ones(2)))
	end

    SetProg.optimize!(model)
	if objective_value === nothing
		@test JuMP.objective_sense(model) == MOI.FEASIBILITY_SENSE
	else
		@test JuMP.objective_sense(model) == MOI.MIN_SENSE
	end
    if feasible
        @test JuMP.termination_status(model) == MOI.OPTIMAL
        @test JuMP.primal_status(model) == MOI.FEASIBLE_POINT
		if objective_value !== nothing
			@test JuMP.objective_value(model) ≈ objective_value atol=config.atol rtol=config.rtol
		end
        return set_test(JuMP.value(◯))
    else
        @test JuMP.termination_status(model) == MOI.INFEASIBLE
        @test JuMP.dual_status(model) == MOI.INFEASIBILITY_CERTIFICATE
        return dual_test(cref1, cref2)
    end
end

function feasible_switched_system_ell_test(optimizer, config, ε=1e-3)
    Q = [1.0 0.0
         0.0 1.0]
    superset = SetProg.Sets.EllipsoidAtOrigin(Symmetric(Q))
    switched_system_test(
        optimizer, config,
        Ellipsoid(symmetric=true, superset=superset),
        √2 + ε, true, 4/3,
        ◯ -> begin
            @test ◯ isa Sets.EllipsoidAtOrigin{Float64}
            @test ◯.Q ≈ Q atol=config.atol rtol=config.rtol
        end,
        (cref1, cref2) -> begin end)
end
function infeasible_switched_system_ell_test(optimizer, config, ε=1e-3)
    Q = [1.0 0.0
         0.0 1.0]
    superset = SetProg.Sets.EllipsoidAtOrigin(Symmetric(Q))
    switched_system_test(
        optimizer, config,
        Ellipsoid(symmetric=true, superset=superset),
        √2 - ε, false, NaN,
        ◯ -> begin end,
        (cref1, cref2) -> begin
            function ispsd(q)
                Q = [q[1] q[2]; q[2] q[3]]
                return all(eigvals(Q) .≥ -config.atol)
            end
            @test ispsd(JuMP.dual(cref1))
            @test ispsd(JuMP.dual(cref2))
        end)
end

function superset(x, d)
    q = SetProg.SumOfSquares.GramMatrix(SetProg.SumOfSquares.SOSDecomposition(x.^d))
    return SetProg.Sets.PolynomialSublevelSetAtOrigin(2d, q)
end

function feasible_switched_system_quad_test(optimizer, config, ε=1e-3)
    @polyvar x[1:2]
    switched_system_test(
        optimizer, config,
        PolySet(symmetric=true, degree=2, superset=superset(x, 1)),
        √2 + ε, true, 8/3,
        ◯ -> begin
            @test ◯ isa Sets.PolynomialSublevelSetAtOrigin{Float64}
            @test polynomial(◯.p) ≈ x[1]^2 + x[2]^2 atol=config.atol rtol=config.rtol
        end,
        (cref1, cref2) -> begin end)
end
function infeasible_switched_system_quad_test(optimizer, config, ε=1e-3)
    @polyvar x[1:2]
    switched_system_test(
        optimizer, config,
        PolySet(symmetric=true, degree=2, superset=superset(x, 1)),
        √2 - ε, false, NaN,
        ◯ -> begin end,
        (cref1, cref2) -> begin
            function ispsd(M)
                return all(eigvals(Matrix(M.Q)) .≥ -config.atol)
            end
            @test ispsd(SetProg.SumOfSquares.moment_matrix(cref1))
            @test ispsd(SetProg.SumOfSquares.moment_matrix(cref2))
        end)
end

function feasible_switched_system_quartic_test(optimizer, config, ε=1e-2)
    @polyvar x[1:2]
    switched_system_test(
        optimizer, config,
        PolySet(symmetric=true, degree=4, superset=superset(x, 2)),
        1.0 + ε, true, 10.001105454190741,
        ◯ -> begin
            @test ◯ isa Sets.PolynomialSublevelSetAtOrigin{Float64}
			α = 11.814054544955727
            @test polynomial(◯.p) ≈ (α+1) * x[1]^4 - 2α * x[1]^2*x[2]^2 + (α+1) * x[2]^4 atol=config.atol rtol=config.rtol
        end,
        (cref1, cref2) -> begin end)
end
function infeasible_switched_system_quartic_test(optimizer, config, ε=2e-1)
    @polyvar x[1:2]
    switched_system_test(
        optimizer, config,
        PolySet(symmetric=true, degree=4, superset=superset(x, 2)),
        1.0 - ε, false, nothing,
        ◯ -> begin end,
        (cref1, cref2) -> begin
            function ispsd(M)
                return all(eigvals(Matrix(M.Q)) .≥ -config.atol)
            end
            @test ispsd(SetProg.SumOfSquares.moment_matrix(cref1))
            @test ispsd(SetProg.SumOfSquares.moment_matrix(cref2))
        end)
end

const switched_system_tests = Dict("feasible_switched_system_ell"     => feasible_switched_system_ell_test,
                                   "infeasible_switched_system_ell"   => infeasible_switched_system_ell_test,
                                   "feasible_switched_system_quad"    => feasible_switched_system_quad_test,
                                   "infeasible_switched_system_quad"  => infeasible_switched_system_quad_test,
                                   "feasible_switched_system_quartic"    => feasible_switched_system_quartic_test,
                                   "infeasible_switched_system_quartic"  => infeasible_switched_system_quartic_test)

@test_suite switched_system
