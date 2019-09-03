using LinearAlgebra
using Test

using SetProg, SetProg.Sets
using Polyhedra
using MultivariatePolynomials

using JuMP
const MOIT = MOI.Test

@testset "Invariant" begin
    config = MOIT.TestConfig()
    @testset "Maximal" begin
        # Q = [1 0
        #      0 1]
        # t = √det(Q) = 1
        Q = [1.0, 0.0, 1.0]
        t = 1.0
        @testset "Homogeneous" begin
            @testset "Ellipsoid" begin
                Tests.maximal_invariant_ell_homogeneous_test(
                    bridged_mock(mock -> MOI.Utilities.mock_optimize!(mock, [Q; t])),
                    config)
            end
            @testset "Convex Quadratic" begin
                Tests.maximal_convex_invariant_quad_homogeneous_test(
                    bridged_mock(mock -> MOI.Utilities.mock_optimize!(mock, [1.0, 1.0, 0.0, 2.0, 0.0, 2.0, 2.0])),
                    config)
            end
        end
    end
    @testset "Minimal" begin
        # Q = [0.5 0
        #      0   0.5]
        # t = √det(Q) = 0.5
        Q = [0.5, 0.0, 0.5]
        t = 0.5
        @testset "Homogeneous" begin
            @testset "Ellipsoid" begin
                Tests.minimal_invariant_ell_homogeneous_test(
                    bridged_mock(mock -> MOI.Utilities.mock_optimize!(mock, [Q; t])),
                    config)
            end
            @testset "Quadratic" begin
                @testset "Non-convex" begin
                    Tests.minimal_invariant_quad_homogeneous_test(
                        bridged_mock(mock -> MOI.Utilities.mock_optimize!(mock, [0.5, 0.5, 0.0])),
                        config)
                end
                @testset "Convex" begin
                    Tests.minimal_convex_invariant_quad_homogeneous_test(
                        bridged_mock(mock -> MOI.Utilities.mock_optimize!(mock, [[0.5, 0.5, 0.0]; 2Q; 2t])),
                        config)
                end
            end
        end
    end
end
