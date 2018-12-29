using LinearAlgebra
using Test

using SetProg, SetProg.Sets
using Polyhedra
using MultivariatePolynomials

using JuMP
const MOI = JuMP.MOI
const MOIT = MOI.Test

function mock(mock_optimize!::Function)
    mock = MOI.Utilities.MockOptimizer(JuMP.JuMPMOIModel{Float64}())
    MOI.Utilities.set_mock_optimize!(mock, mock_optimize!)
    return mock
end

function square_test(optimizer::MOI.AbstractOptimizer, config, inner::Bool,
                     variable::SetProg.AbstractVariable, metric::Function,
                     objective_value, set_test)
    MOI.empty!(optimizer)
    model = JuMP.direct_model(optimizer)

    □ = polyhedron(HalfSpace([1, 0], 1.0) ∩ HalfSpace([-1, 0], 1) ∩ HalfSpace([0, 1], 1) ∩ HalfSpace([0, -1], 1))

    @variable(model, ◯, variable)
    if inner
        cref = @constraint(model, ◯ ⊆ □)
    else
        cref = @constraint(model, □ ⊆ ◯)
    end
    @objective(model, inner ? MOI.MAX_SENSE : MOI.MIN_SENSE,
               metric(volume(◯)))

    SetProg.optimize!(model)
    @test JuMP.termination_status(model) == MOI.OPTIMAL
    @test JuMP.objective_sense(model) == MOI.MAX_SENSE
    @test JuMP.objective_value(model) ≈ objective_value atol=config.atol rtol=config.rtol
    set_test(JuMP.value(◯))
end

function john_homogeneous_square_test(optimizer, config)
    square_test(optimizer, config, true,
                Ellipsoid(symmetric=true, dimension=2),
                nth_root, 1.0,
                ◯ -> begin
                    @test ◯ isa Sets.Polar{Float64, Sets.EllipsoidAtOrigin{Float64}}
                    @test Sets.polar(◯).Q ≈ Symmetric([1.0 0.0; 0.0 1.0]) atol=config.atol rtol=config.rtol
                end)
end

function john_nonhomogeneous_ell_square_test(optimizer, config)
    square_test(optimizer, config, true,
                Ellipsoid(point=SetProg.InteriorPoint([0.0, 0.0])),
                nth_root, 1.0,
                ◯ -> begin
                    @test ◯ isa Sets.PerspectiveDual{Float64, Sets.PerspectiveInteriorEllipsoid{Float64, Float64}}
                    z = Sets.perspective_variable(◯)
                    x, y = Sets.space_variables(◯)
                    ◯_dual = Sets.perspective_dual(◯)
                    @test ◯_dual.p ≈ -z^2 + x^2 + y^2 atol=config.atol rtol=config.rtol
                    @test ◯_dual.Q ≈ Symmetric([1.0 0.0; 0.0 1.0]) atol=config.atol rtol=config.rtol
                    @test ◯_dual.b ≈ [0.0, 0.0] atol=config.atol rtol=config.rtol
                    @test ◯_dual.β ≈ -1.0 atol=config.atol rtol=config.rtol
                    @test ◯_dual.H ≈ [-1.0 0.0 0.0
                                      0.0 1.0 0.0
                                      0.0 0.0 1.0] atol=config.atol rtol=config.rtol
                end)
end

function john_nonhomogeneous_quad_square_test(optimizer, config)
    square_test(optimizer, config, true,
                PolySet(degree=2, convex=true, point=SetProg.InteriorPoint([0.0, 0.0])),
                set -> L1_heuristic(set, [1.0, 1.0]),
                8/3,
                ◯ -> begin
                    @test ◯ isa Sets.PerspectiveDual{Float64, Sets.PerspectiveConvexPolynomialSet{Float64, Float64}}
                    z = Sets.perspective_variable(◯)
                    x, y = Sets.space_variables(◯)
                    ◯_dual = Sets.perspective_dual(◯)
                    @test ◯_dual.p ≈ -z^2 + x^2 + y^2 atol=config.atol rtol=config.rtol
                end)
end

function löwner_homogeneous_square_test(optimizer, config)
    square_test(optimizer, config, false,
                Ellipsoid(symmetric=true, dimension=2),
                nth_root, 0.5,
                ◯ -> begin
                    @test ◯ isa Sets.EllipsoidAtOrigin
                    @test ◯.Q ≈ Symmetric([0.5 0.0; 0.0 0.5]) atol=config.atol rtol=config.rtol
                end)
end

function quartic_inner_homogeneous_square_test(optimizer, config)
    square_test(optimizer, config, true,
                PolySet(symmetric=true, degree=4, dimension=2, convex=true),
                nth_root, 1.0,
                ◯ -> begin
                    @test ◯ isa Sets.Polar{Float64, Sets.ConvexPolynomialSublevelSetAtOrigin{Float64}}
                    ◯_polar = Sets.polar(◯)
                    @test ◯_polar.degree == 4
                    x, y = variables(◯_polar.p)
                    @test polynomial(◯_polar.p) == x^4 + 2x^3*y + 3x^2*y^2 + 2x*y^3 + y^4
                    convexity_proof = Sets.convexity_proof(◯)
                    @test convexity_proof.n == 6
                    @test convexity_proof.Q == ones(21)
                end)
end

function quartic_outer_homogeneous_square_test(optimizer, config)
    square_test(optimizer, config, false,
                PolySet(symmetric=true, degree=4, dimension=2, convex=true),
                nth_root,
                1.0,
                ◯ -> begin
                    @test ◯ isa Sets.ConvexPolynomialSublevelSetAtOrigin{Float64}
                    @test ◯.degree == 4
                    x, y = variables(◯.p)
                    @test polynomial(◯.p) == x^4 + 2x^3*y + 3x^2*y^2 + 2x*y^3 + y^4
                    convexity_proof = Sets.convexity_proof(◯)
                    @test convexity_proof.n == 6
                    @test convexity_proof.Q == ones(21)
                end)
end

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
                john_homogeneous_square_test(mock(mock -> MOI.Utilities.mock_optimize!(mock, [Q; t])),
                                             config)
            end
            @testset "Non-homogeneous" begin
                @testset "Ellipsoid" begin
                    β = -1.0
                    b = [0.0, 0.0]
                    john_nonhomogeneous_ell_square_test(mock(mock -> MOI.Utilities.mock_optimize!(mock, [Q; β; b; t])),
                                                        config)
                end
                @testset "PolySet" begin
                    β = 1.0
                    b = [0.0, 0.0]
                    john_nonhomogeneous_quad_square_test(mock(mock -> begin
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
            löwner_homogeneous_square_test(mock(mock -> MOI.Utilities.mock_optimize!(mock, [0.5, 0.0, 0.5, 0.5])),
                                           config)
        end
    end
    @testset "Quartic" begin
        @testset "Inner" begin
            # The PSD matrix for the variable  is 3 x 3 so 3 * (3+1) / 2 = 6
            # The PSD matrix for the convexity is 6 x 6 so 6 * (6+1) / 2 = 21
            # entries
            # 1 variable for t
            # hence 28 variables
            quartic_inner_homogeneous_square_test(mock(mock -> MOI.Utilities.mock_optimize!(mock, ones(28))),
                                                  config)
        end
        @testset "Outer" begin
            quartic_outer_homogeneous_square_test(mock(mock -> MOI.Utilities.mock_optimize!(mock, ones(28))),
                                                  config)
        end
    end
end

const squaretests = Dict("john_homogeneous_square" => john_homogeneous_square_test,
                         "john_nonhomogeneous_ell_square" => john_nonhomogeneous_ell_square_test,
                         "john_nonhomogeneous_quad_square" => john_nonhomogeneous_quad_square_test,
                         "löwner_homogeneous_square" => löwner_homogeneous_square_test,
                         "quartic_inner_homogeneous" => quartic_inner_homogeneous_square_test,
                         "quartic_outer_homogeneous" => quartic_outer_homogeneous_square_test)

MOIT.@moitestset square
