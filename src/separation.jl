"""
    Disjoint()

Use `@constraint(model, [S, T] in Disjoint())` to certify disjointness.
For quadratic and polynomial sublevel sets, the certificate is a positive
combination of their defining inequalities. `margin` subtracts a positive
reference polynomial; its default of one certifies strict separation.
An explicit zero margin permits shared boundary points.
For a polytope with a fixed normal fan and an exterior ellipsoid, separation
uses second-order cones and `margin` reduces the complementary ellipsoid's
radius.
"""
struct Disjoint <: JuMP.AbstractVectorSet end

"""
    SProcedureMultiplier(; degree = 0, base_name = "λ")

A nonnegative S-procedure multiplier, constructed when its set constraint is
loaded. Degree zero creates a nonnegative scalar; a positive even degree
creates an SOS polynomial. Read its fitted multiplier with `JuMP.value`.
"""
mutable struct SProcedureMultiplier
    degree::Int
    base_name::String
    variable::Any
    function SProcedureMultiplier(; degree::Integer = 0, base_name = "λ")
        degree >= 0 && iseven(degree) ||
            throw(ArgumentError("the multiplier degree must be nonnegative and even"))
        return new(Int(degree), base_name, nothing)
    end
end
JuMP.value(multiplier::SProcedureMultiplier) = JuMP.value(multiplier.variable)

struct DisjointConstraint{S,T,K} <: SetConstraint
    first::S
    second::T
    kws::K
end
function JuMP.build_constraint(error::Function, sets::AbstractVector, ::Disjoint; kws...)
    length(sets) == 2 || error("Disjoint expects two sets")
    return DisjointConstraint(sets[1], sets[2], kws)
end
need_variablify(c::DisjointConstraint) =
    need_variablify(c.first) || need_variablify(c.second)
variablify(c::DisjointConstraint) =
    DisjointConstraint(variablify(c.first), variablify(c.second), c.kws)
function clear_spaces(c::DisjointConstraint)
    clear_spaces(c.first)
    clear_spaces(c.second)
end
Sets.perspective_variable(c::DisjointConstraint) = synchronize_perspective(
    Sets.perspective_variable(c.first),
    Sets.perspective_variable(c.second),
)
function create_spaces(c::DisjointConstraint, spaces::Spaces)
    return merge_spaces(
        spaces,
        create_spaces(c.first, spaces),
        create_spaces(c.second, spaces),
    )
end
set_space(space::Space, ::DisjointConstraint) = set_space(space, PrimalSpace)
set_space(space::Space, ::InclusionConstraint{<:SetVariableRef,<:Sets.ShiftedEllipsoid}) =
    set_space(space, PrimalSpace)
set_space(space::Space, ::InclusionConstraint{<:Sets.ShiftedEllipsoid,<:SetVariableRef}) =
    set_space(space, PrimalSpace)
JuMP.function_string(mode, c::DisjointConstraint) = string(c.first, " ∩ ", c.second)
JuMP.in_set_string(mode, ::DisjointConstraint) = " = ∅"

_fixed_set(::Sets.ShiftedEllipsoid{<:Number}) = true
_fixed_set(::Sets.ShiftedEllipsoid) = false
function _multiplier(model, multiplier::SProcedureMultiplier, variables)
    if multiplier.variable === nothing
        multiplier.variable = if iszero(multiplier.degree)
            @variable(model, lower_bound = 0, base_name = multiplier.base_name)
        else
            sos_extension().polynomial_multiplier(
                model,
                variables,
                multiplier.degree,
                multiplier.base_name,
            )
        end
    end
    return multiplier.variable
end
function _multiplier(model, multiplier::Real, variables)
    multiplier >= 0 || throw(ArgumentError("an S-procedure multiplier must be nonnegative"))
    return multiplier
end
function _multiplier(model, multiplier::JuMP.VariableRef, variables)
    if !(JuMP.has_lower_bound(multiplier) && JuMP.lower_bound(multiplier) >= 0)
        @constraint(model, multiplier >= 0)
    end
    return multiplier
end
function _scaling(model, first, second, scaling, variables = nothing)
    fixed_first, fixed_second = _fixed_set(first), _fixed_set(second)
    if scaling === nothing
        scaling = fixed_first || fixed_second ? SProcedureMultiplier() : 1
    end
    !fixed_first &&
        !fixed_second &&
        !(scaling isa Real) &&
        throw(ArgumentError("two variable sets require a fixed S-procedure multiplier"))
    return _multiplier(model, scaling, variables), fixed_first
end

function _matrix_certificate(model, matrix)
    cref = psd_constraint(model, Symmetric(matrix))
    certificates = get!(model.ext, :SetProgCertificates) do
        Dict{Any,Any}()
    end
    certificates[cref] = matrix
    return cref
end
"""
    certificate_matrix(constraint)

Return the fitted PSD certificate of a set constraint. For an SOS certificate,
return its fitted Gram matrix. The model must have been optimized.
"""
certificate_matrix(cref::SetConstraintRef) =
    certificate_matrix(transformed_constraint(cref))
function certificate_matrix(cref::JuMP.ConstraintRef)
    certificates = get(cref.model.ext, :SetProgCertificates, Dict{Any,Any}())
    if haskey(certificates, cref)
        return JuMP.value.(certificates[cref])
    end
    object = JuMP.constraint_object(cref)
    if object.set isa MOI.LessThan
        return reshape([object.set.upper - JuMP.value(object.func)], 1, 1)
    elseif object.set isa MOI.Nonnegatives
        return Diagonal(JuMP.value.(object.func))
    elseif object.set isa MOI.SecondOrderCone
        values = JuMP.value.(object.func)
        t, x = first(values), values[2:end]
        return Symmetric([t x'; x t * Matrix{Float64}(I, length(x), length(x))])
    end
    if object.set isa MOI.PositiveSemidefiniteConeTriangle
        values = JuMP.value.(object.func)
        n = MOI.side_dimension(object.set)
        matrix = zeros(eltype(values), n, n)
        k = 0
        for j = 1:n, i = 1:j
            matrix[i, j] = matrix[j, i] = values[k+=1]
        end
        return Symmetric(matrix)
    end
    return sos_extension().polynomial_certificate_matrix(cref)
end

function JuMP.add_constraint(
    model::JuMP.Model,
    c::InclusionConstraint{<:Sets.ShiftedEllipsoid,<:Sets.ShiftedEllipsoid},
    name::String = "",
)
    return _quadratic_separation(model, c.subset, c.supset, false; c.kws...)
end
function JuMP.add_constraint(
    model::JuMP.Model,
    c::DisjointConstraint{<:Sets.ShiftedEllipsoid,<:Sets.ShiftedEllipsoid},
    name::String = "",
)
    return _quadratic_separation(model, c.first, c.second, true; c.kws...)
end
function _quadratic_separation(
    model,
    first,
    second,
    disjoint;
    S_procedure_scaling = nothing,
    margin = disjoint ? 1 : 0,
)
    if S_procedure_scaling isa SProcedureMultiplier && S_procedure_scaling.degree > 0
        return sos_extension()._polynomial_separation(
            model,
            first,
            second,
            disjoint;
            S_procedure_scaling,
            margin,
        )
    end
    λ, fixed_first = _scaling(model, first, second, S_procedure_scaling)
    A, B = Sets._perspective_cat(first), Sets._perspective_cat(second)
    matrix = if disjoint
        fixed_first ? λ * A + B : A + λ * B
    else
        fixed_first ? λ * A - B : A - λ * B
    end
    return _matrix_certificate(
        model,
        matrix - margin * Matrix{Float64}(I, size(matrix, 1), size(matrix, 2)),
    )
end

_as_shifted(set::Sets.ShiftedEllipsoid) = set
function _as_shifted(set::Sets.Ellipsoid)
    return Sets.ShiftedEllipsoid(set.Q, zeros(size(set.Q, 1)), -1.0)
end
function JuMP.add_constraint(
    model::JuMP.Model,
    c::DisjointConstraint{
        <:Union{Sets.Ellipsoid,Sets.ShiftedEllipsoid},
        <:Union{Sets.Ellipsoid,Sets.ShiftedEllipsoid},
    },
    name::String = "",
)
    return _quadratic_separation(
        model,
        _as_shifted(c.first),
        _as_shifted(c.second),
        true;
        c.kws...,
    )
end

# A set constraint can produce several independent conic certificates.
function certificate_matrix(crefs::AbstractVector{<:JuMP.ConstraintRef})
    return cat(certificate_matrix.(crefs)...; dims = (1, 2))
end

const _PolytopeVariableRef = SetVariableRef{<:JuMP.AbstractModel,<:Polytope}
const _VertexPolytope = Sets.PolarOf{<:Sets.Piecewise{<:Any,<:Sets.PolarPoint}}
const _FixedEllipsoid = Union{Sets.Ellipsoid{<:Real},Sets.ShiftedEllipsoid{<:Real}}

# Both gauge and vertex coordinates admit conic separation certificates.
# Leave the space undecided so other constraints can select a representation.
set_space(
    space::Space,
    ::InclusionConstraint{<:Sets.ShiftedEllipsoid{<:Real},<:_PolytopeVariableRef},
) = space
set_space(
    space::Space,
    ::InclusionConstraint{<:Sets.Ellipsoid{<:Real},<:_PolytopeVariableRef},
) = space
set_space(
    space::Space,
    ::DisjointConstraint{<:_PolytopeVariableRef,<:Sets.ShiftedEllipsoid{<:Real}},
) = space
set_space(
    space::Space,
    ::DisjointConstraint{<:Sets.ShiftedEllipsoid{<:Real},<:_PolytopeVariableRef},
) = space

function _ellipsoid_geometry(set::Sets.ShiftedEllipsoid{<:Real})
    factor = cholesky(set.Q)
    center = -(factor \ set.b)
    radius = sqrt(dot(set.b, factor \ set.b) - set.β)
    return factor, center, radius
end

function _check_normal_fan(inner::Sets.Piecewise)
    partition = inner.polytope
    isempty(rays(partition)) &&
    isempty(lines(partition)) &&
    isempty(hyperplanes(partition)) &&
    all(h -> h.β > 0, halfspaces(partition)) || throw(
        ArgumentError(
            "the partition must be a bounded polytope with the origin in its interior",
        ),
    )
    return nothing
end

function JuMP.add_constraint(
    model::JuMP.Model,
    c::InclusionConstraint{<:_FixedEllipsoid,<:_VertexPolytope},
    name::String = "",
)
    return _ellipsoid_in_polytope(model, _as_shifted(c.subset), c.supset; c.kws...)
end

function _ellipsoid_in_polytope(
    model,
    ellipsoid,
    polytope;
    margin = 0,
    S_procedure_scaling::Nothing = nothing,
)
    inner = polar(polytope)
    _check_normal_fan(inner)
    factor, center, radius = _ellipsoid_geometry(ellipsoid)
    residuals = JuMP.AffExpr[]
    # The fixed partition supplies the extreme normal directions. The support
    # function of the variable polytope is linear on each partition cone, so
    # these finitely many tests imply containment of the whole ellipsoid.
    for point in points(inner.polytope)
        normal = Polyhedra.coord(point)
        iszero(normal) && continue
        piece = findfirst(cone -> normal in cone, inner.pieces)
        support = dot(center, normal) + radius * sqrt(dot(normal, factor \ normal))
        # The polar family contains the origin. A positive margin also keeps
        # its support function strictly positive in every normal direction.
        bound = max(0, support) + margin * norm(normal)
        push!(residuals, dot(inner.sets[piece].a, normal) - bound)
    end
    return @constraint(model, residuals in MOI.Nonnegatives(length(residuals)))
end

function JuMP.add_constraint(
    model::JuMP.Model,
    c::DisjointConstraint{<:_VertexPolytope,<:Sets.ShiftedEllipsoid{<:Real}},
    name::String = "",
)
    return _polytope_disjoint(model, c.first, c.second; c.kws...)
end
function JuMP.add_constraint(
    model::JuMP.Model,
    c::DisjointConstraint{<:Sets.ShiftedEllipsoid{<:Real},<:_VertexPolytope},
    name::String = "",
)
    return _polytope_disjoint(model, c.second, c.first; c.kws...)
end

function _polytope_disjoint(
    model,
    polytope,
    unsafe;
    margin = 1,
    S_procedure_scaling::Nothing = nothing,
)
    inner = polar(polytope)
    _check_normal_fan(inner)
    # The unsafe quadratic describes an ellipsoid's exterior. Its complement
    # is convex, so testing all variable vertices is exact. Eroding its radius
    # by margin gives strict separation using only second-order cones.
    complement = Sets.ShiftedEllipsoid(Symmetric(-unsafe.Q), -unsafe.b, -unsafe.β)
    factor, center, radius = _ellipsoid_geometry(complement)
    constraints = JuMP.ConstraintRef[
        @constraint(
            model,
            [radius - margin; factor.U * (vertex.a - center)] in SecondOrderCone()
        ) for vertex in inner.sets
    ]
    # A polar set also contains the origin. For a centered ellipsoid this
    # follows from each SOC; otherwise its radius bound must be checked too.
    offset = norm(factor.U * center)
    if !iszero(offset)
        push!(constraints, @constraint(model, margin <= radius - offset))
    end
    return constraints
end

const _LinearPiecewise = Sets.Piecewise{<:Any,<:Sets.PolarPoint}

function JuMP.add_constraint(
    model::JuMP.Model,
    c::InclusionConstraint{<:_FixedEllipsoid,<:_LinearPiecewise},
    name::String = "",
)
    return _ellipsoid_in_gauge(model, _as_shifted(c.subset), c.supset; c.kws...)
end
function _ellipsoid_in_gauge(
    model,
    ellipsoid,
    gauge;
    margin = 0,
    S_procedure_scaling::Nothing = nothing,
)
    _check_normal_fan(gauge)
    factor, center, radius = _ellipsoid_geometry(ellipsoid)
    support_map = radius * inv(Matrix(factor.L))
    # The gauge is the maximum of its linear pieces. Its value is at most
    # 1 - margin throughout the initial ellipsoid exactly when each support
    # bound holds; no S-procedure multiplier is needed.
    return [
        @constraint(
            model,
            [1 - margin - dot(center, piece.a); support_map * piece.a] in SecondOrderCone()
        ) for piece in gauge.sets
    ]
end

function JuMP.add_constraint(
    model::JuMP.Model,
    c::DisjointConstraint{<:_LinearPiecewise,<:Sets.ShiftedEllipsoid{<:Real}},
    name::String = "",
)
    return _gauge_disjoint(model, c.first, c.second; c.kws...)
end
function JuMP.add_constraint(
    model::JuMP.Model,
    c::DisjointConstraint{<:Sets.ShiftedEllipsoid{<:Real},<:_LinearPiecewise},
    name::String = "",
)
    return _gauge_disjoint(model, c.second, c.first; c.kws...)
end
function _gauge_disjoint(
    model,
    gauge,
    unsafe;
    margin = 1,
    S_procedure_scaling::Nothing = nothing,
)
    _check_normal_fan(gauge)
    complement = Sets.ShiftedEllipsoid(Symmetric(-unsafe.Q), -unsafe.b, -unsafe.β)
    _ellipsoid_geometry(complement) # The complement must be a bounded ellipsoid.
    # Every finite homogeneous gauge contains a neighborhood of zero. It
    # cannot lie strictly inside an ellipsoid which does not contain zero.
    if complement.β >= 0
        return @constraint(model, 0 <= -1)
    end
    residuals = JuMP.AffExpr[]
    for point in points(gauge.polytope)
        normal = Polyhedra.coord(point)
        iszero(normal) && continue
        piece = findfirst(cone -> normal in cone, gauge.pieces)
        q = dot(normal, complement.Q * normal)
        b = dot(complement.b, normal)
        root = sqrt(b^2 - complement.β * q)
        # Gauge of the complementary ellipsoid, evaluated in this fixed
        # direction. The two formulas avoid cancellation for either sign.
        bound = b >= 0 ? (root + b) / (-complement.β) : q / (root - b)
        push!(residuals, dot(gauge.sets[piece].a, normal) - bound - margin * norm(normal))
    end
    return @constraint(model, residuals in MOI.Nonnegatives(length(residuals)))
end
