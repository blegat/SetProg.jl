using LinearAlgebra
using Test

using SetProg, SetProg.Sets
using Polyhedra
using MultivariatePolynomials

using JuMP
const MOIT = MOI.Test

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
                Tests.ci_ell_nonhomogeneous_test(bridged_mock(mock -> MOI.Utilities.mock_optimize!(mock, [Q; β; b; t])),
                                                 config)
            end
            @testset "PolySet" begin
                β = -1.0
                b = [0.0, 0.0]
                Q = [1.0, -0.5, 1.0]

                Tests.ci_quad_nonhomogeneous_test(bridged_mock(mock -> begin
                         # β+1 b[1] b[2]
                         #  .  Q[1] Q[2]
                         #  .   .   Q[3]
                         MOI.Utilities.mock_optimize!(mock, [β+1; b[1]; Q[1]; b[2]; Q[2]; Q[3]; 2Q])
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
        ci_quartic_α = -1/8
        ci_quartic_β =  1/4
        ci_quartic_hess = 6 * [2.0, ci_quartic_α, 2.0, ci_quartic_α, 2.0,
                               2.0, 2.0, ci_quartic_α, ci_quartic_α, 2.0]
        sol = [1.0; ci_quartic_α; 6 - 2ci_quartic_β; ci_quartic_β; ci_quartic_α; 1.0;
               ci_quartic_hess]
        Tests.ci_quartic_homogeneous_test(bridged_mock(mock -> MOI.Utilities.mock_optimize!(mock, sol)),
                                          config)
    end
end
