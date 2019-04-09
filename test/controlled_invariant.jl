using LinearAlgebra
using Test

using SetProg, SetProg.Sets
using Polyhedra
using MultivariatePolynomials

using JuMP
const MOIT = MOI.Test

const ci_quartic_α = -1.4529635030551264
const ci_quartic_β =  0.25613759221607024
const ci_quartic_γ = -8.717781018330758
const ci_quartic_hess = [12.0, ci_quartic_γ, 12.0, ci_quartic_γ, 12.0,
                         12.0, 12.0, ci_quartic_γ, ci_quartic_γ, 12.0]

const ci_quartic_obj = 0.25369938382997853

@testset "Controlled invariant" begin
    config = MOIT.TestConfig()
    @testset "Ellipsoid" begin
        # Q = [1 0
        #      0 1]
        # t = √det(Q) = 1
        Q = [1.0, -1/4, 1.0]
        t = √15/4
        @testset "Homogeneous" begin
            Tests.ci_ell_homogeneous_test(bridged_mock(mock -> MOI.Utilities.mock_optimize!(mock, [Q; t])),
                                          config)
        end
        @testset "Non-homogeneous" begin
            @testset "Ellipsoid" begin
                β = -1.0
                b = [0.0, 0.0]
                # The constraint (0)z² + (b[2])zx₁ + (-Q[1,2] - 0.25 Q[2,2])x₁² is SOS
                # use the certificate monomials [x₁] so it needs a 1x1 PSD matrix so it creates 1 variable whose value is 0
                sos = 0.0
                Tests.ci_ell_nonhomogeneous_test(bridged_mock(mock -> MOI.Utilities.mock_optimize!(mock, [Q; β; b; sos; t])),
                                                 config)
            end
            @testset "PolySet" begin
                β = -1.0
                b = [0.0, 0.0]
                Q = [1.0, -0.4933095968, 1.0]

                Tests.ci_quad_nonhomogeneous_test(bridged_mock(mock -> begin
                         # β+1 b[1] b[2]
                         #  .  Q[1] Q[2]
                         #  .   .   Q[3]
                         MOI.Utilities.mock_optimize!(mock, [β+1; b[1]; Q[1]; b[2]; Q[2]; Q[3]; 2Q; 0.24331])
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
        Tests.ci_quartic_homogeneous_test(bridged_mock(mock -> MOI.Utilities.mock_optimize!(mock, sol)),
                                          config)
    end
end
