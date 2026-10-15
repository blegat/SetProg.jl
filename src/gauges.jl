set_space(space::Space, ::InclusionConstraint{<:Sets.AbstractEllipsoid,<:SetVariableRef}) =
    set_space(space, PrimalSpace)

"""
    gauge_value(set, x)

Evaluate the homogeneous gauge of a fitted set. Ellipsoid and polynomial
representations use their degree-two and degree-`degree` sublevel functions.
"""
gauge_value(set::Sets.Ellipsoid, x::AbstractVector) = dot(x, set.Q * x)
gauge_value(set::Sets.PolarPoint, x::AbstractVector) = dot(set.a, x)
function gauge_value(set::Sets.Piecewise, x::AbstractVector)
    i = findfirst(piece -> x in piece, set.pieces)
    i === nothing && throw(ArgumentError("the point lies in none of the pieces"))
    return gauge_value(set.sets[i], x)
end
(set::Sets.AbstractSet)(x::AbstractVector{<:Real}) = gauge_value(set, x)
Base.Matrix(set::Sets.Ellipsoid) = Matrix(set.Q)

# Linear pieces are shared by ordinary and polar polyhedral families.
const _LinearPolytope = Sets.Piecewise{<:Any,<:Sets.PolarPoint}

function _ray_cone(C::AbstractMatrix)
    n = size(C, 1)
    rep = vrep(
        [zeros(n)],
        Line{Float64,Vector{Float64}}[],
        [Ray(collect(c)) for c in eachcol(C)],
    )
    return polyhedron(rep)
end

function _linear_piecewise(coefficients, cones)
    pieces = _ray_cone.(cones)
    sets = [Sets.PolarPoint(collect(a)) for a in coefficients]
    T = eltype(first(sets).a)
    return Sets.Piecewise{T,eltype(sets),Float64,Nothing,eltype(pieces)}(
        sets,
        nothing,
        pieces,
        [Tuple{Int,Vector{Float64}}[] for _ in sets],
    )
end

function _polytope_variable_set(model, family::Polytope, n)
    partition = family.piecewise
    floor = something(family.min_scale, 1)
    transform = family.transform === nothing ? Matrix{Float64}(I, n, n) : family.transform
    if partition isa Symbol
        signs =
            partition === :orthant ? [ones(n)] :
            vec([Float64[s...] for s in Iterators.product(ntuple(_ -> (-1, 1), n)...)])
        w = @variable(model, [1:n], lower_bound = floor)
        cones = [transform * Diagonal(sign) for sign in signs]
        coefficients = [inv(transform') * (sign .* w) for sign in signs]
    else
        cones = [transform * Matrix{Float64}(C) for C in partition]
        a = @variable(model, [1:length(cones), 1:n])
        coefficients = [a[i, :] for i in eachindex(cones)]
        if family.symmetric
            cones = vcat(cones, [-C for C in cones])
            coefficients = vcat(coefficients, [-c for c in coefficients])
        end
        # A piece must dominate every other piece on its cone. On shared
        # boundaries this also enforces continuity, without an adjacency graph.
        for (i, C) in enumerate(cones), ray in eachcol(C)
            for b in coefficients
                @constraint(model, dot(coefficients[i] - b, ray) >= 0)
            end
            @constraint(model, dot(coefficients[i], ray) >= floor * norm(ray, Inf))
        end
    end
    set = _linear_piecewise(coefficients, cones)
    result = family.polar ? Polyhedra.polar(set) : set
    if partition === :orthant
        result = Sets.Intersection(result, _ray_cone(transform))
    end
    return result
end

function _domain_invariant(A, source, destination)
    return all(r -> A * Polyhedra.coord(r) in destination, rays(source)) && all(
        l ->
            A * Polyhedra.coord(l) in destination && -A * Polyhedra.coord(l) in destination,
        lines(source),
    )
end

function JuMP.build_constraint(
    error::Function,
    subset::Sets.LinearImage{<:Sets.Intersection},
    sup::PowerSet{<:Sets.Intersection},
)
    _domain_invariant(subset.A, subset.set.domain, sup.set.domain) ||
        throw(ArgumentError("the linear map must preserve the polytope domain"))
    return JuMP.build_constraint(
        error,
        Sets.LinearImage(subset.set.set, subset.A),
        PowerSet(sup.set.set),
    )
end

function _cone_directions(cone)
    return Iterators.flatten((
        (Polyhedra.coord(r) for r in rays(cone)),
        (Polyhedra.coord(l) for l in lines(cone)),
        (-Polyhedra.coord(l) for l in lines(cone)),
    ))
end

function _linear_inclusion!(model, source::_LinearPolytope, destination::_LinearPolytope, A)
    for (s, cone) in zip(source.sets, source.pieces), ray in _cone_directions(cone)
        for d in destination.sets
            @constraint(model, dot(s.a, ray) - dot(d.a, A * ray) >= 0)
        end
    end
    return nothing
end

function JuMP.add_constraint(
    model::JuMP.Model,
    constraint::InclusionConstraint{
        <:Sets.LinearImage{<:_LinearPolytope},
        <:_LinearPolytope,
    },
    ::String = "",
)
    _linear_inclusion!(model, constraint.subset.set, constraint.supset, constraint.subset.A)
    return nothing
end

function JuMP.add_constraint(
    model::JuMP.Model,
    constraint::InclusionConstraint{<:_LinearPolytope,<:_LinearPolytope},
    ::String = "",
)
    if constraint.subset.polytope !== nothing && constraint.supset.polytope !== nothing
        return _add_piecewise_inclusion(model, constraint.subset, constraint.supset)
    end
    _linear_inclusion!(
        model,
        constraint.subset,
        constraint.supset,
        Matrix{Float64}(
            I,
            Sets.dimension(constraint.subset),
            Sets.dimension(constraint.subset),
        ),
    )
    return nothing
end

# Polar inclusion reverses the two sets and transposes the dynamics. In
# particular it needs no inverse, so singular dynamics are supported too.
function JuMP.build_constraint(
    error::Function,
    subset::Sets.LinearImage{<:Sets.PolarOf{<:_LinearPolytope}},
    sup::PowerSet{<:Sets.PolarOf{<:_LinearPolytope}},
)
    return JuMP.build_constraint(
        error,
        Sets.LinearImage(Polyhedra.polar(sup.set), subset.A'),
        PowerSet(Polyhedra.polar(subset.set)),
    )
end

function JuMP.add_constraint(
    model::JuMP.Model,
    constraint::InclusionConstraint{<:_LinearPolytope,<:Polyhedra.HalfSpace},
    ::String = "",
)
    h = constraint.supset
    h.β > 0 || throw(ArgumentError("a polar membership bound must be positive"))
    for (s, cone) in zip(constraint.subset.sets, constraint.subset.pieces),
        ray in _cone_directions(cone)

        @constraint(model, h.β * dot(s.a, ray) - dot(h.a, ray) >= 0)
    end
    return nothing
end

gauge_value(set::_LinearPolytope, x::AbstractVector) =
    maximum(dot(s.a, x) for s in set.sets)
function gauge_value(set::Sets.PolarOf{<:_LinearPolytope}, x::AbstractVector)
    core = Polyhedra.polar(set)
    return max(
        0,
        maximum(
            dot(x, Polyhedra.coord(r)) / dot(s.a, Polyhedra.coord(r)) for
            (s, cone) in zip(core.sets, core.pieces) for r in rays(cone)
        ),
    )
end
gauge_value(set::Sets.Intersection, x::AbstractVector) =
    x in set.domain ? gauge_value(set.set, x) : Inf

function JuMP.build_constraint(
    error::Function,
    subset::Sets.Intersection,
    sup::PowerSet{<:Sets.Intersection},
)
    n = Sets.dimension(subset)
    return JuMP.build_constraint(
        error,
        Sets.LinearImage(subset, Matrix{Float64}(I, n, n)),
        sup,
    )
end

function JuMP.build_constraint(
    error::Function,
    subset::Sets.Intersection,
    sup::PowerSet{<:Polyhedra.HalfSpace},
)
    h = sup.set
    if subset.set isa Sets.PolarOf{<:_LinearPolytope}
        # The support of a downward closed polar, restricted to the positive
        # orthant, is evaluated at the positive part of the normal.
        h = HalfSpace(max.(h.a, 0), h.β)
    end
    return JuMP.build_constraint(error, subset.set, PowerSet(h))
end

"""
    Sets.Piecewise(coefficients, cones)

A convex homogeneous piecewise-linear sublevel set. Each coefficient vector
is the linear piece on the cone generated by the corresponding matrix's
columns. Coefficients must already satisfy continuity and convexity.
"""
Sets.Piecewise(
    coefficients::AbstractVector{<:AbstractVector},
    cones::AbstractVector{<:AbstractMatrix},
) = _linear_piecewise(coefficients, cones)

function _polytope_image_inclusion(error, subset, superset)
    B = superset.A
    size(B, 1) == size(B, 2) && rank(B) == size(B, 1) || throw(
        ArgumentError(
            "the target image must be invertible for this polytope representation",
        ),
    )
    return JuMP.build_constraint(
        error,
        Sets.LinearImage(subset.set, B \ subset.A),
        PowerSet(superset.set),
    )
end
function JuMP.build_constraint(
    error::Function,
    subset::Sets.LinearImage{<:_LinearPolytope},
    sup::PowerSet{<:Sets.LinearImage{<:_LinearPolytope}},
)
    return _polytope_image_inclusion(error, subset, sup.set)
end
function JuMP.build_constraint(
    error::Function,
    subset::Sets.LinearImage{<:Sets.Intersection},
    sup::PowerSet{<:Sets.LinearImage{<:Sets.Intersection}},
)
    return _polytope_image_inclusion(error, subset, sup.set)
end

function JuMP.build_constraint(
    error::Function,
    subset::Sets.LinearImage{<:_LinearPolytope},
    sup::PowerSet{<:_LinearPolytope};
    kws...,
)
    return InclusionConstraint(subset, sup.set, kws)
end
