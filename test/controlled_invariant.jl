using LinearAlgebra
using Test

using SetProg, SetProg.Sets
using Polyhedra
using MultivariatePolynomials

using JuMP
const MOI = JuMP.MOI

@testset "Controlled invariant" begin
    □ = polyhedron(HalfSpace([1, 0], 1.0) ∩ HalfSpace([-1, 0], 1) ∩ HalfSpace([0, 1], 1) ∩ HalfSpace([0, -1], 1))
    Δt = 0.5
    A = [1.0 Δt]
    E = [1.0 0.0]
    function test(inner::Bool, variable::SetProg.AbstractVariable,
                  metric::Function, mock_optimize!::Function, objective_value,
                  set_test)
        mock = MOI.Utilities.MockOptimizer(JuMP.JuMPMOIModel{Float64}())
        model = JuMP.direct_model(mock);
        MOI.Utilities.set_mock_optimize!(mock, mock_optimize!)

        @variable(model, ◯, variable)
        if inner
            cref = @constraint(model, ◯ ⊆ □)
        else
            cref = @constraint(model, □ ⊆ ◯)
        end
        @constraint(model, A * ◯ ⊆ E * ◯)
        @objective(model, inner ? MOI.MAX_SENSE : MOI.MIN_SENSE,
                   metric(volume(◯)))

        SetProg.optimize!(model)
        @test JuMP.termination_status(model) == MOI.OPTIMAL
        @test JuMP.objective_sense(model) == MOI.MAX_SENSE
        @test JuMP.objective_value(model) ≈ objective_value rtol=1e-5
        set_test(JuMP.value(◯))
    end
    @testset "Ellipsoid" begin
        # Q = [1 0
        #      0 1]
        # t = √det(Q) = 1
        Q = [1.0, -1/4, 1.0]
        t = √15/4
        @testset "Homogeneous" begin
            test(true,
                 Ellipsoid(symmetric=true, dimension=2), nth_root,
                 mock -> MOI.Utilities.mock_optimize!(mock, [Q; t]),
                 t,
                 ◯ -> begin
                     @test ◯ isa Sets.Polar{Float64, Sets.EllipsoidAtOrigin{Float64}}
                     @test Sets.polar(◯).Q == Symmetric([1.0 -1/4; -1/4 1.0])
                 end)
        end
        @testset "Non-homogeneous" begin
            @testset "Ellipsoid" begin
                β = -1.0
                b = [0.0, 0.0]
                test(true,
                     Ellipsoid(point=SetProg.InteriorPoint([0.0, 0.0])),
                     nth_root,
                     mock -> MOI.Utilities.mock_optimize!(mock, [Q; β; b; zeros(3); t]),
                     t,
                     ◯ -> begin
                         @test ◯ isa Sets.PerspectiveDual{Float64, Sets.PerspectiveInteriorEllipsoid{Float64,Float64}}
                         z = Sets.perspective_variable(◯)
                         x, y = Sets.space_variables(◯)
                         ◯_dual = Sets.perspective_dual(◯)
                         @test ◯_dual.p == -z^2 + x^2 - x*y/2 + y^2
                         @test ◯_dual.Q == Symmetric([1.0 -1/4; -1/4 1.0])
                         @test ◯_dual.b == [0.0, 0.0]
                         @test ◯_dual.β == β
                         @test ◯_dual.H ≈ [-1.0 0.0 0.0
                                            0.0 1.0 0.0
                                            0.0 0.0 1.0]
                     end)
            end
            @testset "PolySet" begin
                β = -1.0
                b = [0.0, 0.0]
                Q = [1.0, -1/2, 1.0]

                test(true,
                     PolySet(degree=2, convex=true, point=SetProg.InteriorPoint([0.0, 0.0])),
                     set -> L1_heuristic(set, [1.0, 1.0]),
                     mock -> begin
                         # β+1 b[1] b[2]
                         #  .  Q[1] Q[2]
                         #  .   .   Q[3]
                         MOI.Utilities.mock_optimize!(mock, [β+1; b[1]; Q[1]; b[2]; Q[2]; Q[3]; 2Q; 0.0; 0.0; 0.24331])
                     end,
                     8/3,
                     ◯ -> begin
                         @test ◯ isa Sets.PerspectiveDual{Float64, Sets.PerspectiveConvexPolynomialSet{Float64,Float64}}
                         z = Sets.perspective_variable(◯)
                         x, y = Sets.space_variables(◯)
                         ◯_dual = Sets.perspective_dual(◯)
                         @test ◯_dual.p == -z^2 + x^2 - x*y + y^2
                     end)
            end
        end
    end
    @testset "Quartic" begin
        # The PSD matrix for the variable  is 3 x 3 so 3 * (3+1) / 2 = 6
        # The PSD matrix for the convexity is 6 x 6 so 6 * (6+1) / 2 = 21
        # entries
        # 1 variable for t
        # hence 28 variables
        α = -8.7398984946
        β = -1.4566
        ϵ =  0.2521560726
        δ =  6 - 2ϵ
        ξ =  0.258304355
        η = -8.7397089279
        sol = [1.0; β; δ; ϵ; β; 1.0; 0.0; 0.0; 12.0; 0.0; α; 12.0; 0.0; α; 12.0;
               12.0; 0.0; 12.0; η; η; 12.0; zeros(6); ξ]
        test(true,
             PolySet(symmetric=true, degree=4, convex=true),
             set -> L1_heuristic(set, [1.0, 1.0]),
             mock -> MOI.Utilities.mock_optimize!(mock, sol),
             64/15,
             ◯ -> begin
                 @test ◯ isa Sets.Polar{Float64, Sets.ConvexPolynomialSublevelSetAtOrigin{Float64}}
                 @test Sets.polar(◯).degree == 4
                 x, y = variables(Sets.polar(◯).p)
                 q = x^4 - 2.9132x^3*y + 6x^2*y^2 - 2.9132x*y^3 + y^4
                 @test polynomial(Sets.polar(◯).p) ≈ q rtol=1e-4
                 convexity_proof = Sets.convexity_proof(◯)
                 @test convexity_proof.n == 6
                 hessian = [0.0; 0.0; 12.0; 0.0; α; 12.0; 0.0; α; 12.0; 12.0;
                            0.0; 12.0; η; η; 12.0; zeros(6)]
                 @test convexity_proof.Q == hessian
             end)
    end
end
