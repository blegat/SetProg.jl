using LinearAlgebra
using Test

using SetProg, SetProg.Sets
using Polyhedra
using MultivariatePolynomials

import DynamicPolynomials

using JuMP

const □ = polyhedron(HalfSpace([1, 0], 1.0) ∩ HalfSpace([-1, 0], 1) ∩ HalfSpace([0, 1], 1) ∩ HalfSpace([0, -1], 1))
const ◇ = polyhedron(convexhull([1.0, 0], [0, 1], [-1, 0], [0, -1]), lib)

function square_test(optimizer, config::MOI.Test.Config,
                     inner::Bool, variable::SetProg.AbstractVariable,
                     metric::Function, objective_value, set_test)
    model = _model(optimizer)

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
                    @test ◯ isa Sets.Polar{Float64, Sets.Ellipsoid{Float64}}
                    @test Sets.polar(◯).Q ≈ Symmetric([1.0 0.0; 0.0 1.0]) atol=config.atol rtol=config.rtol
                end)
end

function john_nonhomogeneous_ell_square_test(optimizer, config)
    square_test(optimizer, config, true,
                Ellipsoid(point=SetProg.InteriorPoint([0.0, 0.0])),
                nth_root, 1.0,
                ◯ -> begin
                    @test ◯ isa Sets.PerspectiveDual{Float64, Sets.Householder{Float64, Sets.ShiftedEllipsoid{Float64}, Float64}}
                    z = Sets.perspective_variable(◯)
                    x, y = Sets.space_variables(◯)
                    ◯_dual = Sets.perspective_dual(◯)
                    @test ◯_dual.p ≈ -z^2 + x^2 + y^2 atol=config.atol rtol=config.rtol
                    @test Sets._householder(◯_dual.h) ≈ [-1.0 0.0 0.0
                                                          0.0 1.0 0.0
                                                          0.0 0.0 1.0] atol=config.atol rtol=config.rtol
                    @test ◯_dual.set.Q ≈ Symmetric([1.0 0.0; 0.0 1.0]) atol=config.atol rtol=config.rtol
                    @test ◯_dual.set.b ≈ [0.0, 0.0] atol=config.atol rtol=config.rtol
                    @test ◯_dual.set.β ≈ -1.0 atol=config.atol rtol=config.rtol
                end)
end

function john_nonhomogeneous_quad_square_test(optimizer, config)
    square_test(optimizer, config, true,
                PolySet(degree=2, convex=true, point=SetProg.InteriorPoint([0.0, 0.0])),
                set -> L1_heuristic(set, [1.0, 1.0]),
                8/3,
                ◯ -> begin
                    @test ◯ isa Sets.PerspectiveDual{Float64, Sets.Householder{Float64, Sets.ConvexPolynomialSet{Float64, SetProg.Sets.MonoBasis, Float64}, Float64}}
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
                    @test ◯ isa Sets.Ellipsoid
                    @test ◯.Q ≈ Symmetric([0.5 0.0
                                           0.0 0.5]) atol=config.atol rtol=config.rtol
                end)
end

function piecewise_semiell_inner_homogeneous_◇_square_test(optimizer, config)
    square_test(optimizer, config, true,
                Ellipsoid(symmetric=true, piecewise=◇),
                set -> L1_heuristic(set), 4,
                p◯ -> begin
                    @test p◯ isa Sets.Polar
                    ◯ = p◯.set
                    @test ◯ isa Sets.Piecewise{Float64, Sets.Ellipsoid{Float64}}
                    @test length(◯.sets) == 4
                    Q1 = Symmetric([ 1.0  1.0
                                     1.0  1.0])
                    Q2 = Symmetric([ 1.0 -1.0
                                    -1.0  1.0])
                    @test ◯.sets[1].Q ≈ Q1 atol=config.atol rtol=config.rtol
                    @test ◯.sets[2].Q ≈ Q2 atol=config.atol rtol=config.rtol
                    @test ◯.sets[3].Q ≈ Q2 atol=config.atol rtol=config.rtol
                    @test ◯.sets[4].Q ≈ Q1 atol=config.atol rtol=config.rtol
                end)
end

function piecewise_semiell_inner_homogeneous_□_square_test(optimizer, config)
    square_test(optimizer, config, true,
                Ellipsoid(symmetric=true, piecewise=□),
                set -> L1_heuristic(set), 8/3,
                p◯ -> begin
                    @test p◯ isa Sets.Polar
                    ◯ = p◯.set
                    @test ◯ isa Sets.Piecewise{Float64, Sets.Ellipsoid{Float64}}
                    @test length(◯.sets) == 4
                    Q = Symmetric([1.0 0.0
                                   0.0 1.0])
                    @test ◯.sets[1].Q ≈ Q atol=config.atol rtol=config.rtol
                    @test ◯.sets[2].Q ≈ Q atol=config.atol rtol=config.rtol
                    @test ◯.sets[3].Q ≈ Q atol=config.atol rtol=config.rtol
                    @test ◯.sets[4].Q ≈ Q atol=config.atol rtol=config.rtol
                end)
end

const quartic_inner_poly = [3.1518541833100864, -0.1617384194869734]
const quartic_inner_obj = 6.447419478140056
const quartic_inner_α = 5.6567546886722795
const quartic_inner_convexity = [12.0, 0.0, quartic_inner_α, 0.0, quartic_inner_poly[1]+2quartic_inner_poly[2],
                                 quartic_inner_α, 8.48516455194103, 0.0, 0.0, 12.0]

function quartic_inner_homogeneous_square_test(optimizer, config)
    square_test(optimizer, config, true,
                PolySet(symmetric=true, degree=4, dimension=2, convex=true),
                nth_root, quartic_inner_obj,
                ◯ -> begin
                    @test ◯ isa Sets.Polar{Float64, Sets.ConvexPolySet{Float64, SetProg.Sets.MonoBasis, Float64}}
                    ◯_polar = Sets.polar(◯)
                    @test ◯_polar.degree == 4
                    x, y = variables(◯_polar.p)
                    @test polynomial(◯_polar.p) ≈ x^4 + quartic_inner_convexity[5]*x^2*y^2 + y^4 atol=config.atol rtol=config.rtol
                    convexity_proof = Sets.convexity_proof(◯)
                    @test convexity_proof.n == 4
                    @test convexity_proof.Q ≈ quartic_inner_convexity atol=config.atol rtol=config.rtol
                end)
end

const quartic_outer_β = 0.30177574048813055
const quartic_outer_γ = 0.5936049698923986
const quartic_outer_λ = -0.09857757888257276
const quartic_outer_obj = 1.611854896946893
const quartic_outer_α = 0.7928996242545062
const quartic_outer_convexity = [3.621308885857567, 0.0, quartic_outer_α, 0.0, 0.08578956499151169,
                                 quartic_outer_α, 1.5, 0.0, 0.0, 3.6212933687704307]

function quartic_outer_homogeneous_square_test(optimizer, config)
    square_test(optimizer, config, false,
                PolySet(symmetric=true, degree=4, dimension=2, convex=true),
                nth_root,
                quartic_outer_obj,
                ◯ -> begin
                    @test ◯ isa Sets.ConvexPolySet{Float64, SetProg.Sets.MonoBasis, Float64}
                    @test ◯.degree == 4
                    x, y = variables(◯.p)
                    @test polynomial(◯.p) ≈ quartic_outer_β*x^4 + (quartic_outer_γ+2quartic_outer_λ)*x^2*y^2 + quartic_outer_β*y^4 atol=config.atol rtol=config.rtol
                    convexity_proof = Sets.convexity_proof(◯)
                    @test convexity_proof.n == 4
                    @test convexity_proof.Q ≈ quartic_outer_convexity atol=config.atol rtol=config.rtol
                end)
end

const square_tests = Dict(
    "john_homogeneous_square" =>
     john_homogeneous_square_test,
    "john_nonhomogeneous_ell_square" =>
     john_nonhomogeneous_ell_square_test,
    "john_nonhomogeneous_quad_square" =>
     john_nonhomogeneous_quad_square_test,
    "löwner_homogeneous_square" =>
     löwner_homogeneous_square_test,
    "piecewise_semiell_inner_homogeneous_◇_square" =>
     piecewise_semiell_inner_homogeneous_◇_square_test,
    "piecewise_semiell_inner_homogeneous_□_square" =>
     piecewise_semiell_inner_homogeneous_□_square_test,
    "quartic_inner_homogeneous_square" =>
     quartic_inner_homogeneous_square_test,
    "quartic_outer_homogeneous_square" =>
     quartic_outer_homogeneous_square_test
)

@test_suite square
