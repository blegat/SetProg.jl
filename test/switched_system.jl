using LinearAlgebra
using Test

using SetProg, SetProg.Sets
using Polyhedra
using MultivariatePolynomials
import DynamicPolynomials

using JuMP
const MOIT = MOI.Test

@testset "Switched System" begin
    config = MOIT.TestConfig()
    @testset "Ellipsoid" begin
        @testset "Feasible" begin
            optimize!(mock) = MOIU.mock_optimize!(
                mock, MOI.OPTIMAL, (MOI.FEASIBLE_POINT, zeros(3)))
            Tests.feasible_switched_system_ell_test(bridged_mock(optimize!), config)
        end
        @testset "Infeasible" begin
            optimize!(mock) = MOIU.mock_optimize!(
                mock, MOI.INFEASIBLE, tuple(),
                (MOI.VectorOfVariables, SetProg.SumOfSquares.PositiveSemidefinite2x2ConeTriangle) => [[0.00282108, 0.0, 0.00282108]],
                (MOI.VectorAffineFunction{Float64}, SetProg.SumOfSquares.PositiveSemidefinite2x2ConeTriangle) => [[0.00117187, 0.0, 2.82044], [2.82044, 0.0, 0.00117187]]
            )
            Tests.infeasible_switched_system_ell_test(bridged_mock(optimize!), config)
        end
    end
    @testset "Quadratic" begin
        @testset "Feasible" begin
            optimize!(mock) = MOIU.mock_optimize!(
                mock, MOI.OPTIMAL, (MOI.FEASIBLE_POINT, zeros(3)))
            Tests.feasible_switched_system_quad_test(bridged_mock(optimize!), config)
        end
        # Blocked by a bug in MockOptimizer which does not implement support for MomentMatrixAttribute correctly
#        @testset "Infeasible" begin
#            function optimize!(mock)
#                MOI.set(mock, MOI.TerminationStatus(), MOI.INFEASIBLE)
#                MOI.set(mock, MOI.PrimalStatus(), MOI.NO_SOLUTION)
#                MOI.set(mock, MOI.DualStatus(), MOI.INFEASIBILITY_CERTIFICATE)
#                DynamicPolynomials.@polyvar x[1:2]
#                for (F, S) in MOI.get(mock, MOI.ListOfConstraints())
#                    if F == MOI.VectorAffineFunction{Float64}
#                        cis = MOI.get(mock, MOI.ListOfConstraintIndices{F, S}())
#                        function _moment_matrix(q)
#                            SetProg.SumOfSquares.build_moment_matrix(q, monovec(x))
#                        end
#                        MOI.set(mock, SetProg.SumOfSquares.MomentMatrixAttribute(),
#                                cis[1], _moment_matrix([0.00336272, 0.0, 5.64559]))
#                        MOI.set(mock, SetProg.SumOfSquares.MomentMatrixAttribute(),
#                                cis[2], _moment_matrix([5.64559, 0.0, 0.00336272]))
#                    end
#                end
#            end
#            Tests.infeasible_switched_system_quad_test(bridged_mock(optimize!), config)
#        end
    end
    @testset "Quartic" begin
        @testset "Feasible" begin
			α = 11.814054544955727
            optimize!(mock) = MOIU.mock_optimize!(
                mock, MOI.OPTIMAL, (MOI.FEASIBLE_POINT, [α, 0.0, 0.0, -α, -0.0, α]))
            Tests.feasible_switched_system_quartic_test(bridged_mock(optimize!), config)
        end
        @testset "Infeasible" begin
            # Blocked for the same reason than quad
        end
    end
end
