using LinearAlgebra
using Test

using SetProg, SetProg.Sets
using Polyhedra
using MultivariatePolynomials

using JuMP
const MOIT = MOI.Test

const quartic_inner_poly = [3.1518541833100864, -0.1617384194869734]
const quartic_inner_obj = 6.447419478140056
const quartic_inner_α = 5.6567546886722795
const quartic_inner_convexity = [12.0, 0.0, quartic_inner_α, 0.0, quartic_inner_poly[1]+2quartic_inner_poly[2],
                                 quartic_inner_α, 8.48516455194103, 0.0, 0.0, 12.0]

const quartic_outer_β = 0.30177574048813055
const quartic_outer_γ = 0.5936049698923986
const quartic_outer_λ = -0.09857757888257276
const quartic_outer_obj = 1.611854896946893
const quartic_outer_α = 0.7928996242545062
const quartic_outer_convexity = [3.621308885857567, 0.0, quartic_outer_α, 0.0, 0.08578956499151169,
                                 quartic_outer_α, 1.5, 0.0, 0.0, 3.6212933687704307]

@testset "Square" begin
    config = MOIT.TestConfig()
    @testset "Ellipsoid" begin
        @testset "John" begin
            # Q = [1 0
            #      0 1]
            # t = √det(Q) = 1
            Q = [1.0, 0.0, 1.0]
            t = 1.0
            @testset "Homogeneous" begin
                Tests.john_homogeneous_square_test(bridged_mock(mock -> MOI.Utilities.mock_optimize!(mock, [Q; t])),
                                                   config)
            end
            @testset "Non-homogeneous" begin
                @testset "Ellipsoid" begin
                    β = -1.0
                    b = [0.0, 0.0]
                    Tests.john_nonhomogeneous_ell_square_test(bridged_mock(mock -> MOI.Utilities.mock_optimize!(mock, [Q; β; b; t])),
                                                              config)
                end
                @testset "PolySet" begin
                    β = 1.0
                    b = [0.0, 0.0]
                    Tests.john_nonhomogeneous_quad_square_test(bridged_mock(mock -> begin
                        # β-1 b[1] b[2]
                        #  .  Q[1] Q[2]
                        #  .   .   Q[3]
                        MOI.Utilities.mock_optimize!(mock, [β-1; b[1]; Q[1]; b[2]; Q[2]; Q[3]; 2Q])
                    end), config)
                end
            end
        end
        @testset "Löwner" begin
            # Q = [√2  0
            #       0 √2]
            # t = √det(Q) = 2                                                                  Q11  Q12  Q22  t
            Tests.löwner_homogeneous_square_test(bridged_mock(mock -> MOI.Utilities.mock_optimize!(mock, [0.5, 0.0, 0.5, 0.5])),
                                                 config)
        end
    end
    @testset "Quartic" begin
        @testset "Inner" begin
            # The PSD matrix for the variable  is 3 x 3 so 3 * (3+1) / 2 = 6
            # The PSD matrix for the convexity is 4 x 4 so 4 * (5+1) / 2 = 10
            # entries
            # 1 variable for t
            # hence 17 variables
            sol = [1.0; 0.0; quartic_inner_poly; 0.0; 1.0;
                   quartic_inner_convexity; quartic_inner_obj]
            Tests.quartic_inner_homogeneous_square_test(bridged_mock(mock -> MOI.Utilities.mock_optimize!(mock, sol)),
                                                        config)
        end
        @testset "Outer" begin
            sol = [quartic_outer_β; 0.0; quartic_outer_γ; quartic_outer_λ; 0.0; quartic_outer_β;
                   quartic_outer_convexity; quartic_outer_obj]
            Tests.quartic_outer_homogeneous_square_test(bridged_mock(mock -> MOI.Utilities.mock_optimize!(mock, sol)),
                                                        config)
        end
    end
end
