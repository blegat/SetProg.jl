"""
    homogeneity_degree(family)

The degree of the homogeneous gauge used to represent a set family.
"""
homogeneity_degree(::Union{Ellipsoid,Sets.Ellipsoid}) = 2
homogeneity_degree(
    ::Union{Polytope,Sets.PolarPoint,Sets.Piecewise{<:Any,<:Sets.PolarPoint}},
) = 1
homogeneity_degree(family::Union{PolySet,Sets.PolySet,Sets.ConvexPolySet}) = family.degree
homogeneity_degree(set::Union{Sets.Polar,Sets.Intersection}) = homogeneity_degree(set.set)
function homogeneity_degree(families::AbstractVector{<:AbstractVariable})
    isempty(families) && throw(ArgumentError("at least one set family is required"))
    degrees = unique(homogeneity_degree.(families))
    length(degrees) == 1 ||
        throw(ArgumentError("every node must have the same gauge degree"))
    return only(degrees)
end

"""
    with_dimension(family, dimension; variables = nothing)

Fill the dimension of a homogeneous gauge family, validating explicit data.
`variables` supplies a shared polynomial coordinate system when needed.
"""
function with_dimension(family::AbstractVariable, n::Integer; variables = nothing)
    n > 0 || throw(ArgumentError("the dimension must be positive"))
    family.dimension === nothing || family.dimension == n || throw(DimensionMismatch())
    return _with_dimension(family, Int(n), variables)
end
function _with_dimension(family::Ellipsoid, n, vars)
    family.symmetric || throw(ArgumentError("a homogeneous ellipsoid must be symmetric"))
    family.point === nothing && family.superset === nothing ||
        throw(ArgumentError("point and superset are not supported for homogeneous gauges"))
    return Ellipsoid(symmetric = true, dimension = n, piecewise = family.piecewise)
end
function _with_dimension(family::Polytope, n, vars)
    family.piecewise !== nothing ||
        throw(ArgumentError("a polytope gauge requires a partition"))
    return Polytope(
        symmetric = family.symmetric,
        dimension = n,
        piecewise = family.piecewise,
        transform = family.transform,
        polar = family.polar,
        min_scale = family.min_scale,
    )
end
function _with_dimension(family::PolySet, n, vars)
    family.symmetric ||
        throw(ArgumentError("a homogeneous polynomial set must be symmetric"))
    family.point === nothing && family.superset === nothing ||
        throw(ArgumentError("point and superset are not supported for homogeneous gauges"))
    x =
        family.variables === nothing ? (vars === nothing ? new_polyvars(n) : vars) :
        family.variables
    length(x) == n || throw(DimensionMismatch())
    return PolySet(
        symmetric = true,
        degree = family.degree,
        dimension = n,
        convex = family.convex,
        variables = x,
        basis = family.basis,
    )
end
"""
    unit_set(family)

The unit gauge used to normalize a set variable in a homogeneous certificate.
"""
function unit_set(f::Ellipsoid)
    ball = Sets.Ellipsoid(Symmetric(Matrix{Float64}(I, f.dimension, f.dimension)))
    f.piecewise === nothing && return ball
    return Sets.Piecewise([ball for _ in halfspaces(f.piecewise)], f.piecewise)
end
unit_set(f::Polytope) = Sets.Piecewise(
    [Sets.PolarPoint(h.a / h.β) for h in halfspaces(f.piecewise)],
    f.piecewise,
)
unit_set(f::PolySet) = sos_extension().unit_set(f)

"""
    add_function_variables!(model, family, dimension, node)

Create a homogeneous set variable, named for `node`, for a certificate model.
"""
function add_function_variables!(
    model::JuMP.Model,
    family::AbstractVariable,
    n::Integer,
    node::Integer,
)
    variable = with_dimension(family, n)
    return @variable(model, variable_type = variable, base_name = "S_$node")
end

"""
    add_normalization!(model, family, S; conditioning_bound = Inf)

Exclude the zero gauge by requiring `S ⊆ unit_set(S.set)`. Orthant and
polyhedral families with bounded weights impose their normalization when
creating the set. A finite ellipsoid conditioning bound also requires
`Q ≤ conditioning_bound * I`.
"""
function add_normalization!(
    model::JuMP.Model,
    family::AbstractVariable,
    S::SetVariableRef;
    conditioning_bound::Real = Inf,
)
    conditioning_bound >= 1 || throw(ArgumentError("conditioning_bound must be at least 1"))
    if S.set isa Polytope && S.set.piecewise isa Union{Symbol,AbstractVector}
        isfinite(conditioning_bound) &&
            throw(ArgumentError("conditioning bounds require ellipsoids"))
        return nothing
    end
    unit = unit_set(S.set)
    @constraint(model, S ⊆ unit)
    if isfinite(conditioning_bound)
        S.set isa Ellipsoid && S.set.piecewise === nothing ||
            throw(ArgumentError("conditioning bounds require an ellipsoid"))
        bound = Sets.Ellipsoid(
            Symmetric(
                conditioning_bound * Matrix{Float64}(I, S.set.dimension, S.set.dimension),
            ),
        )
        @constraint(model, bound ⊆ S)
    end
    return nothing
end
set_space(space::Space, ::InclusionConstraint{<:Sets.AbstractEllipsoid,<:SetVariableRef}) =
    set_space(space, PrimalSpace)

"""
    add_decrease!(model, family, source, destination, A; scale = 1)

Require `V_destination(A * x) ≤ scale * V_source(x)`, using set inclusion.
"""
function add_decrease!(
    model::JuMP.Model,
    family::AbstractVariable,
    source::SetVariableRef,
    destination::SetVariableRef,
    A::AbstractMatrix;
    scale::Real = 1,
)
    scale >= 0 || throw(ArgumentError("the scale must be nonnegative"))
    if iszero(scale)
        iszero(A) || @constraint(model, zero(JuMP.AffExpr) >= 1)
        return nothing
    end
    γ = scale^(1 / homogeneity_degree(family))
    @constraint(model, (A / γ) * source ⊆ destination)
    return nothing
end

"""
    solution_value(family, variable)

The fitted set with solver values substituted for its variables.
"""
solution_value(::AbstractVariable, variable::SetVariableRef) = JuMP.value(variable)

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

"""
    decrease_slack(family, source, destination, A; scale = 1)

Measure the fitted gauge's decrease condition. A nonnegative slack certifies
that the modeled edge inequality holds. Polynomial and piecewise-quadratic
slacks require a separate optimization and are not implemented.
"""
function decrease_slack(
    ::Ellipsoid,
    source::Sets.Ellipsoid,
    destination::Sets.Ellipsoid,
    A;
    scale = 1,
)
    return minimum(eigvals(Symmetric(Matrix(scale * source.Q - A' * destination.Q * A))))
end
function decrease_slack(family::AbstractVariable, source, destination, A; scale = 1)
    throw(ArgumentError("decrease slacks are not implemented for $(typeof(family))"))
end

"""
    supports_dual(family)

Whether the fitted homogeneous gauge has an implemented dual gauge family.
"""
supports_dual(::AbstractVariable) = false
supports_dual(f::Ellipsoid) = f.symmetric && f.piecewise === nothing
supports_dual(f::Polytope) = f.piecewise === :orthant

"""
    dual_family(family)

The family of dual homogeneous gauges. On the positive orthant this uses
blocking duality between positive simplices and boxes.
"""
dual_family(f::Ellipsoid) =
    supports_dual(f) ? f :
    throw(ArgumentError("piecewise ellipsoid duals are not implemented"))
function dual_family(f::Polytope)
    supports_dual(f) ||
        throw(ArgumentError("no dual family is implemented for this partition"))
    return Polytope(
        symmetric = f.symmetric,
        dimension = f.dimension,
        piecewise = f.piecewise,
        transform = f.transform,
        polar = !f.polar,
        min_scale = f.min_scale,
    )
end

"""
    dual_value(set)

The fitted dual gauge, represented as a set in the primal space.
"""
dual_value(set::Sets.Ellipsoid) = Sets.Ellipsoid(Symmetric(inv(set.Q)))

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

JuMP.value(set::Sets.Intersection) = Sets.Intersection(JuMP.value(set.set), set.domain)
needs_polyvars(set::Sets.Intersection) = needs_polyvars(set.set)

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

function JuMP.add_constraint(
    model::JuMP.Model,
    constraint::MembershipConstraint{<:Point,<:Sets.Intersection},
    ::String = "",
)
    point, set = constraint.member, constraint.set
    for h in halfspaces(set.domain)
        @constraint(model, dot(h.a, Polyhedra.coord(point)) <= h.β * scaling(point))
    end
    JuMP.add_constraint(model, JuMP.build_constraint(error, point, set.set))
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

dual_value(set::Sets.Intersection) = Sets.Intersection(Polyhedra.polar(set.set), set.domain)
dual_value(set::Union{_LinearPolytope,Sets.PolarOf{<:_LinearPolytope}}) =
    Polyhedra.polar(set)

function _polytope_slack(source::_LinearPolytope, destination::_LinearPolytope, A, scale)
    return minimum(
        scale * dot(s.a, Polyhedra.coord(r)) -
        gauge_value(destination, A * Polyhedra.coord(r)) for
        (s, cone) in zip(source.sets, source.pieces) for r in rays(cone)
    )
end
function decrease_slack(
    ::Polytope,
    source::_LinearPolytope,
    destination::_LinearPolytope,
    A;
    scale = 1,
)
    return _polytope_slack(source, destination, A, scale)
end
function decrease_slack(
    ::Polytope,
    source::Sets.PolarOf{<:_LinearPolytope},
    destination::Sets.PolarOf{<:_LinearPolytope},
    A;
    scale = 1,
)
    return _polytope_slack(Polyhedra.polar(destination), Polyhedra.polar(source), A', scale)
end
function decrease_slack(
    family::Polytope,
    source::Sets.Intersection,
    destination::Sets.Intersection,
    A;
    scale = 1,
)
    return decrease_slack(family, source.set, destination.set, A; scale)
end

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
