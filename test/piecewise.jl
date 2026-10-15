module TestPiecewiseCones

using Test
using LinearAlgebra
using SetProg
import GLPK

const Sets = SetProg.Sets
const SIGNS = ([1.0, 1.0], [-1.0, 1.0], [1.0, -1.0], [-1.0, -1.0])
const COEFFICIENTS = [sign .* [2.0, 3.0] for sign in SIGNS]
const CONES =
    vrep.(Ref(zeros(1, 2)), permutedims.([Matrix(Diagonal(sign)) for sign in SIGNS]))
const PIECES = Sets.PolarPoint.(COEFFICIENTS)

function check_membership(set, point, feasible)
    model = Model(GLPK.Optimizer)
    set_silent(model)
    @constraint(model, point in set)
    optimize!(model)
    @test termination_status(model) == (feasible ? MOI.OPTIMAL : MOI.INFEASIBLE)
end

@testset "Linear pieces on explicit cones" begin
    set = Sets.piecewise_from_cones(PIECES, CONES)
    @test set.polytope === nothing
    @test Sets.dimension(set) == 2
    @test Sets.space_variables(set) === nothing
    @test Sets.perspective_variable(set) === nothing
    @test !SetProg.needs_polyvars(set)
    @test length(set.pieces) == length(CONES)
    @test all(isempty, set.graph)
    explicit = Sets.Piecewise(PIECES, nothing, set.pieces)
    @test [piece.a for piece in explicit.sets] == COEFFICIENTS
    @test all(isempty, explicit.graph)
    @test_throws DimensionMismatch Sets.Piecewise(PIECES[1:1], nothing, set.pieces)
    @test_throws ArgumentError Sets.Piecewise(PIECES, nothing, nothing)
    @test_throws ArgumentError Sets.piecewise_from_cones(
        Sets.PolarPoint{Float64}[],
        typeof(first(CONES))[],
    )
    for (piece, coefficient, cone, sign) in zip(set.sets, COEFFICIENTS, set.pieces, SIGNS)
        @test piece.a == coefficient
        @test sign in cone
        @test -sign ∉ cone
    end
    for sign in SIGNS
        check_membership(set, sign .* [0.2, 0.1], true)
        check_membership(set, sign .* [0.3, 0.2], false)
    end
    fitted = value(set)
    @test fitted.polytope === nothing
    @test [piece.a for piece in fitted.sets] == COEFFICIENTS
    @test Sets.dimension(fitted) == 2

    integer_set = Sets.piecewise_from_cones(
        [Sets.PolarPoint([2, 3])],
        [vrep(zeros(Int, 1, 2), [1 0; 0 1])],
    )
    @test Sets.dimension(integer_set) == 2
    @test integer_set.sets[1].a == [2, 3]

    # PCC stores generators as columns; Polyhedra matrices store them as rows.
    matrices = [[1.0 1.0; 0.0 1.0]]
    cones = vrep.(Ref(zeros(1, 2)), permutedims.(matrices))
    triangular = Sets.piecewise_from_cones([Sets.PolarPoint([1.0, 1.0])], cones)
    @test [1.0, 0.0] in only(triangular.pieces)
    @test [0.0, 1.0] ∉ only(triangular.pieces)

end

@testset "Linear preimages of explicit cone partitions" begin
    set = Sets.piecewise_from_cones(PIECES, CONES)
    A = [2.0 0.0; 0.0 1.0]
    preimage = SetProg.apply_map(Sets.LinearPreImage(set, A), nothing)
    @test preimage.polytope === nothing
    @test Sets.dimension(preimage) == 2
    @test [piece.a for piece in preimage.sets] == [A' * a for a in COEFFICIENTS]
    for sign in SIGNS
        check_membership(preimage, A \ (sign .* [0.2, 0.1]), true)
        check_membership(preimage, A \ (sign .* [0.3, 0.2]), false)
    end

    # The dimension comes from the preimage cones, even for rectangular maps.
    B = reshape([1.0, 0.0], 2, 1)
    preimage = SetProg.apply_map(Sets.LinearPreImage(set, B), nothing)
    @test Sets.dimension(preimage) == 1
    check_membership(preimage, [0.4], true)
    check_membership(preimage, [-0.6], false)
end

@testset "Partitions from a polytope still map their parent" begin
    partition = polyhedron(
        hrep([1.0 0.0; 0.0 1.0; -1.0 0.0; 0.0 -1.0], ones(4)),
        Polyhedra.DefaultLibrary{Float64}(GLPK.Optimizer),
    )
    set = Sets.piecewise_from_polytope(
        [
            Sets.PolarPoint(collect(a)) for
            a in eachrow([1.0 0.0; 0.0 1.0; -1.0 0.0; 0.0 -1.0])
        ],
        partition,
    )
    @test Sets.dimension(set) == 2
    explicit = Sets.Piecewise(set.sets, partition, set.pieces)
    @test explicit.graph == set.graph
    @test explicit.pieces === set.pieces
    @test Sets.Piecewise(set.sets, partition).graph == set.graph
    A = [2.0 0.0; 0.0 1.0]
    preimage = SetProg.apply_map(Sets.LinearPreImage(set, A), nothing)
    @test preimage.polytope !== nothing
    @test [0.5, 0.0] in preimage.polytope
    @test [0.75, 0.0] ∉ preimage.polytope
    check_membership(preimage, [0.4, 0.0], true)
    check_membership(preimage, [0.6, 0.0], false)
end

end
