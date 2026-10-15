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
    @test all(isempty, explicit.graph)
    @test explicit.polytope === partition
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

@testset "Polytope variables on explicit cones" begin
    @test Polytope(piecewise = CONES).dimension == 2
    @test_throws DimensionMismatch Polytope(dimension = 3, piecewise = CONES)
    @test_throws DimensionMismatch Polytope(
        piecewise = [CONES[1], vrep(zeros(1, 3), Matrix{Float64}(I, 3, 3))],
    )
    @test_throws ArgumentError Polytope(piecewise = typeof(CONES)())

    # The upper-half-plane cones suffice when symmetry supplies the others.
    for symmetric in (false, true), space in (SetProg.PrimalSpace, SetProg.DualSpace)
        cones = symmetric ? CONES[1:2] : CONES
        family = Polytope(symmetric = symmetric, piecewise = cones)
        model = Model(GLPK.Optimizer)
        set_silent(model)
        result = SetProg.variable_set(model, family, space, 2, nothing)
        set = space == SetProg.PrimalSpace ? result : polar(result)
        for (piece, expected) in zip(
            set.sets,
            symmetric ? vcat(COEFFICIENTS[1:2], [-c for c in COEFFICIENTS[1:2]]) :
            COEFFICIENTS,
        )
            @constraint(model, piece.a .== expected)
        end
        optimize!(model)
        @test termination_status(model) == MOI.OPTIMAL
        @test set.polytope === nothing
        @test length(set.sets) == 4
        fitted = value(set)
        for (piece, cone) in zip(fitted.sets, fitted.pieces), ray in rays(cone)
            r = Polyhedra.coord(ray)
            @test dot(piece.a, r) ≈ 2abs(r[1]) + 3abs(r[2])
        end
    end

    # Positive coefficients can still violate continuity across a boundary.
    model = Model(GLPK.Optimizer)
    set_silent(model)
    set = SetProg.variable_set(
        model,
        Polytope(piecewise = CONES),
        SetProg.PrimalSpace,
        2,
        nothing,
    )
    for (i, piece) in enumerate(set.sets)
        @constraint(model, piece.a .== (i == 1 ? [4.0, 3.0] : COEFFICIENTS[i]))
    end
    optimize!(model)
    @test termination_status(model) == MOI.INFEASIBLE

    # A globally linear function is convex, but must also be nonnegative.
    model = Model(GLPK.Optimizer)
    set_silent(model)
    set = SetProg.variable_set(
        model,
        Polytope(piecewise = CONES),
        SetProg.PrimalSpace,
        2,
        nothing,
    )
    for piece in set.sets
        @constraint(model, piece.a .== [1.0, 0.0])
    end
    optimize!(model)
    @test termination_status(model) == MOI.INFEASIBLE

    # These positive, continuous pieces have a concave kink at (1, 1).
    generators = [
        [1.0, 0.0],
        [1.0, 1.0],
        [0.0, 1.0],
        [-1.0, 1.0],
        [-1.0, 0.0],
        [-1.0, -1.0],
        [0.0, -1.0],
        [1.0, -1.0],
    ]
    values = [1.0, 3.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0]
    cones = [
        vrep(zeros(1, 2), permutedims(hcat(generators[i], generators[mod1(i + 1, 8)])))
        for i = 1:8
    ]
    model = Model(GLPK.Optimizer)
    set_silent(model)
    set = SetProg.variable_set(
        model,
        Polytope(piecewise = cones),
        SetProg.PrimalSpace,
        2,
        nothing,
    )
    for (i, piece) in enumerate(set.sets)
        j = mod1(i + 1, 8)
        @constraint(model, dot(piece.a, generators[i]) == values[i])
        @constraint(model, dot(piece.a, generators[j]) == values[j])
    end
    optimize!(model)
    @test termination_status(model) == MOI.INFEASIBLE

    # Line directions must constrain both signs, allowing strip gauges.
    cones = [vrep([Line([0.0, -1.0])], [Ray([sign, 0.0])]) for sign in (-1.0, 1.0)]
    model = Model(GLPK.Optimizer)
    set_silent(model)
    set = SetProg.variable_set(
        model,
        Polytope(piecewise = cones),
        SetProg.PrimalSpace,
        2,
        nothing,
    )
    for (piece, sign) in zip(set.sets, (-1.0, 1.0))
        @constraint(model, piece.a[1] == sign)
        @constraint(model, -1 <= piece.a[2] <= 1)
    end
    @objective(model, Max, sum(piece.a[2] for piece in set.sets))
    optimize!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test objective_value(model) ≈ 0
    @constraint(model, set.sets[1].a[2] == -1)
    optimize!(model)
    @test termination_status(model) == MOI.INFEASIBLE

    # Exercise JuMP's set-variable construction and automatic dimension inference.
    for cones in (
        CONES,
        vrep.(Ref(zeros(Int, 1, 2)), [Matrix(Diagonal(Int.(sign))) for sign in SIGNS]),
        [hrep(-Matrix(Diagonal(sign)), zeros(2)) for sign in SIGNS],
    )
        model = Model(GLPK.Optimizer)
        set_silent(model)
        @variable(model, S, Polytope(piecewise = cones))
        @constraint(model, [0.25, 0.25] in S)
        optimize!(model)
        @test termination_status(model) == MOI.OPTIMAL
        @test Sets.dimension(value(S)) == 2
    end
end

end
