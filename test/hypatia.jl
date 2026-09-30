using LinearAlgebra
using Test

using SetProg, SetProg.Sets
using Polyhedra
using JuMP
import Hypatia
import GLPK

const _□ = polyhedron(
    HalfSpace([1, 0], 1.0) ∩ HalfSpace([-1, 0], 1) ∩ HalfSpace([0, 1], 1) ∩
    HalfSpace([0, -1], 1),
)

# Cones received by Hypatia, `nothing` for the `Zeros` of equality constraints
_hypatia_cones(model) = typeof.(JuMP.unsafe_backend(model).moi_cones)

function _test_volume(variable, inner::Bool, metric, cone, obj, set_test; num_cones = 1)
    model = Model(optimizer_with_attributes(Hypatia.Optimizer, MOI.Silent() => true))
    @variable(model, ◯, variable)
    if inner
        @constraint(model, ◯ ⊆ _□)
    else
        @constraint(model, _□ ⊆ ◯)
    end
    @objective(model, inner ? MOI.MAX_SENSE : MOI.MIN_SENSE, metric(volume(◯)))
    optimize!(model)
    # Hypatia may only reach `ALMOST_OPTIMAL` with its default tolerances
    @test termination_status(model) in [MOI.OPTIMAL, MOI.ALMOST_OPTIMAL]
    @test primal_status(model) in [MOI.FEASIBLE_POINT, MOI.NEARLY_FEASIBLE_POINT]
    @test JuMP.objective_sense(model) == MOI.MAX_SENSE
    @test objective_value(model) ≈ obj atol = 1e-5
    set_test(value(◯))
    cones = _hypatia_cones(model)
    # The cone is given as is to Hypatia instead of being bridged into a
    # `PositiveSemidefiniteConeTriangle` with either a `GeometricMeanCone` or
    # an `ExponentialCone`
    @test count(isequal(MOI.Scaled{cone}), cones) == num_cones
    # The only PSD cones are the ones added explicitly by the model
    psd = MOI.PositiveSemidefiniteConeTriangle
    num_psd =
        num_constraints(model, Vector{VariableRef}, psd) +
        num_constraints(model, Vector{AffExpr}, psd)
    @test count(isequal(MOI.Scaled{psd}), cones) == num_psd
    @test !any(C -> C <: Union{MOI.GeometricMeanCone,MOI.ExponentialCone}, cones)
    return
end

@testset "Hypatia volume" begin
    @testset "$name" for (name, metric, cone, john, löwner) in [
        ("nth_root", nth_root, MOI.RootDetConeTriangle, 1.0, 0.5),
        ("log", log, MOI.LogDetConeTriangle, 0.0, 2log(0.5)),
    ]
        @testset "John homogeneous" begin
            _test_volume(
                Ellipsoid(symmetric = true, dimension = 2),
                true,
                metric,
                cone,
                john,
                ◯ -> begin
                    @test ◯ isa Sets.Polar{Float64,Sets.Ellipsoid{Float64}}
                    @test Sets.polar(◯).Q ≈ I atol = 1e-5
                end,
            )
        end
        @testset "John non-homogeneous" begin
            _test_volume(
                Ellipsoid(point = SetProg.InteriorPoint([0.0, 0.0])),
                true,
                metric,
                cone,
                john,
                ◯ -> begin
                    @test ◯ isa Sets.PerspectiveDual
                    ell = Sets.perspective_dual(◯).set
                    @test ell.Q ≈ I atol = 1e-5
                    @test ell.b ≈ zeros(2) atol = 1e-5
                    @test ell.β ≈ -1 atol = 1e-5
                end,
            )
        end
        @testset "Löwner homogeneous" begin
            _test_volume(
                Ellipsoid(symmetric = true, dimension = 2),
                false,
                metric,
                cone,
                löwner,
                ◯ -> begin
                    @test ◯ isa Sets.Ellipsoid{Float64}
                    @test ◯.Q ≈ I / 2 atol = 1e-5
                end,
            )
        end
    end
end

@testset "Hypatia volume of piecewise semi-ellipsoids" begin
    # The pieces need an LP solver to detect their linearity
    lp_solver = optimizer_with_attributes(
        GLPK.Optimizer,
        MOI.Silent() => true,
        "presolve" => GLPK.GLP_ON,
    )
    lib = Polyhedra.DefaultLibrary{Float64}(lp_solver)
    □ = polyhedron(hrep(_□), lib)
    ◇ = polyhedron(convexhull([1.0, 0], [0, 1], [-1, 0], [0, -1]), lib)
    # The volume heuristic is the sum of the heuristic of each of the 4 pieces.
    # By symmetry of the square, each piece is the John (resp. Löwner) ellipsoid.
    @testset "$name $pieces_name" for (name, metric, cone, john, löwner) in [
            ("nth_root", nth_root, MOI.RootDetConeTriangle, 1.0, 0.5),
            ("log", log, MOI.LogDetConeTriangle, 0.0, 2log(0.5)),
        ],
        (pieces_name, pieces) in [("□", □), ("◇", ◇)]

        @testset "John" begin
            _test_volume(
                Ellipsoid(symmetric = true, piecewise = pieces),
                true,
                metric,
                cone,
                4john,
                ◯ -> begin
                    @test ◯ isa Sets.Polar{Float64,<:Sets.Piecewise}
                    for piece in Sets.polar(◯).sets
                        @test piece.Q ≈ I atol = 1e-3
                    end
                end;
                num_cones = 4,
            )
        end
        @testset "Löwner" begin
            _test_volume(
                Ellipsoid(symmetric = true, piecewise = pieces),
                false,
                metric,
                cone,
                4löwner,
                ◯ -> begin
                    @test ◯ isa Sets.Piecewise
                    for piece in ◯.sets
                        @test piece.Q ≈ I / 2 atol = 1e-3
                    end
                end;
                num_cones = 4,
            )
        end
    end
end

@testset "Log volume" begin
    model = Model()
    @variable(model, S, Ellipsoid(symmetric = true, dimension = 2))
    @test sprint(show, log(volume(S))) == "log(volume(S))"
    @test copy(log(volume(S))) isa SetProg.LogVolume
end
