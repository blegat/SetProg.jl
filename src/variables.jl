abstract type AbstractVariable <: JuMP.AbstractVariable end

abstract type HintPoint end

# It will not really be the center and the center for z = 0 is the the same as for <h, x> = 0
struct CenterPoint{T} <: HintPoint
    h::Vector{T}
end

struct InteriorPoint{T} <: HintPoint
    h::Vector{T}
end
_β(model, ::InteriorPoint) = @variable(model)
_b(model, h::InteriorPoint) = @variable(model, [1:length(h.h)])
function polar_perspective_ellipsoid(model, Q::Symmetric{JuMP.VariableRef},
                                     point::CenterPoint, z, x)
    psd_constraint(model, Q)
    ell = Sets.Ellipsoid(Q)
    return polar_perspective_ellipsoid(ell, point, z, x)
end
function polar_perspective_ellipsoid(model, Q::Symmetric{JuMP.VariableRef},
                                     point::InteriorPoint, z, x)
    n = LinearAlgebra.checksquare(Q)
    @assert n == length(point.h)
    β = @variable(model, base_name="β")
    b = @variable(model, [1:length(point.h)], base_name="b")
    psd_constraint(model, Symmetric([β+1 b'; b Q]))
    ell = Sets.ShiftedEllipsoid(Q, b, β)
        return polar_perspective_ellipsoid(ell, point, z, x)
end

### Polytope ###
"""
    Polytope(; symmetric = false, dimension = nothing, piecewise = nothing,
             domain = nothing)

A polyhedral set variable, optionally restricted to a fixed polyhedral domain.
`domain` accepts a Polyhedra representation or polyhedron, `MOI.Nonnegatives(n)`,
or `MOI.Nonnegatives` to infer the orthant dimension from the model.
Domain-restricted variables currently require primal modeling.
"""
struct Polytope{D} <: AbstractVariable
    symmetric::Bool
    dimension::Union{Nothing,Int}
    piecewise::Union{Polyhedra.Rep,Nothing}
    domain::D
end
Polytope(
    symmetric::Bool,
    dimension::Union{Int,Nothing},
    piecewise::Union{Polyhedra.Rep,Nothing},
) = Polytope(symmetric, dimension, piecewise, nothing)

function Polytope(;
    symmetric::Bool = false,
    dimension::Union{Int,Nothing} = nothing,
    piecewise::Union{Polyhedra.Rep,Nothing} = nothing,
    domain::Union{Nothing,Polyhedra.Rep,MOI.Nonnegatives,Type{MOI.Nonnegatives}} = nothing,
)
    function update_dim(object, dim_fun)
        if object !== nothing
            d = dim_fun(object)
            if dimension === nothing
                dimension = d
            elseif dimension != d
                throw(DimensionMismatch())
            end
        end
    end
    update_dim(piecewise, Polyhedra.fulldim)
    if domain isa Polyhedra.Rep
        update_dim(domain, Polyhedra.fulldim)
    elseif domain isa MOI.Nonnegatives
        update_dim(domain, MOI.dimension)
    end
    return Polytope(symmetric, dimension, piecewise, domain)
end
Sets.space_variables(::Polytope) = nothing

_polytope_domain(domain::Polyhedra.Polyhedron, n) = domain
_polytope_domain(domain::Polyhedra.Representation, n) = polyhedron(domain)
function _polytope_domain(::Union{MOI.Nonnegatives,Type{MOI.Nonnegatives}}, n)
    return polyhedron(hrep(-Matrix{Float64}(I, n, n), zeros(n)))
end

function variable_set(model::JuMP.AbstractModel, ell::Polytope, space::Space,
                      space_dimension, space_polyvars)
    ell.domain === nothing || space == PrimalSpace ||
        throw(ArgumentError("domain-restricted polytopes require primal modeling"))
    n = space_dimension
    if ell.symmetric || ell.domain !== nothing
        if ell.piecewise === nothing
            set = Sets.PolarPoint(@variable(model, [1:n], base_name = "a"))
        else
            hashyperplanes(ell.piecewise) && error("hyperplanes not supported for piecewise")
            a = @variable(model, [1:nhalfspaces(ell.piecewise), 1:n], base_name="a")
            sets = [Sets.PolarPoint(a[i, :]) for i in 1:nhalfspaces(ell.piecewise)]
            set = Sets.Piecewise(sets, ell.piecewise)
            for i in eachindex(set.graph)
                for (j, v) in set.graph[i]
                    if i < j # The constraints are the same for (i, j) and (j, i)
                        λ = @variable(model, base_name="λc[$i,$j]")
                        # Ensures continuity
                        @constraint(model, sets[i].a - sets[j].a .== λ * v)
                        # Ensures convexity
                        # We need to ensure that `dot(sets[i].a, x)` is maximal
                        # over when `piece[i]`; see Proposition 5 of [R21].
                        #
                        # [R21] Raković, S. V.
                        # *Control Minkowski–Lyapunov functions*
                        # Automatica, Elsevier BV, 2021, 128, 109598
                        #
                        # In fact, we only need this to hold locally to prove
                        # convexity of each point so we can just check for every
                        # neighbor.
                        # Because <ai, x> = <aj, x> by the previous constraint,
                        # we can check the gradient.
                        @constraint(model, dot(sets[i].a, v) <= dot(sets[j].a, v), base_name="conv[$i,$j]")
                    end
                end
            end
        end
        if ell.domain !== nothing
            set = Sets.Intersection(set, _polytope_domain(ell.domain, n))
        end
        if space == PrimalSpace
            return set
        else
            @assert space == DualSpace
            return Polyhedra.polar(set)
        end
    else
        error("Non-symmetric polytope not supported yet")
    end
end
function JuMP.value(h::Sets.PolarPoint)
    return Sets.PolarPoint(JuMP.value.(h.a))
end



### Ellipsoid ###
struct Ellipsoid <: AbstractVariable
    point::Union{Nothing, HintPoint}
    symmetric::Bool
    dimension::Union{Nothing, Int}
    guaranteed_psd::Bool # Is it already guaranteed that it is PSD ? e.g. by nth_root or log
    superset::Union{Sets.Ellipsoid, Nothing}
    piecewise::Union{Polyhedra.Rep, Nothing}
end
function Ellipsoid(; point::Union{Nothing, HintPoint}=nothing,
                   symmetric::Bool=false,
                   dimension::Union{Int, Nothing}=nothing,
                   superset::Union{Sets.Ellipsoid, Nothing}=nothing,
                   piecewise::Union{Polyhedra.Rep, Nothing}=nothing)
    function update_dim(object, dim_fun)
        if object !== nothing
            d = dim_fun(object)
            if dimension === nothing
                dimension = d
            elseif dimension != d
                throw(DimensionMismatch())
            end
        end
    end
    update_dim(point, point -> length(point.h))
    update_dim(superset, Sets.dimension)
    update_dim(piecewise, Polyhedra.fulldim)
    return Ellipsoid(point, symmetric, dimension, false, superset, piecewise)
end
Sets.space_variables(::Ellipsoid) = nothing

function variable_set(model::JuMP.AbstractModel, ell::Ellipsoid, space::Space,
                      space_dimension, space_polyvars)
    n = space_dimension
    function new_Q()
        # TODO, we should use constraiend variable instead in case direct mode is used.
        Q = @variable(model, [1:n, 1:n], Symmetric, base_name="Q")
        if !ell.guaranteed_psd
            psd_constraint(model, Q)
        end
        return Q
    end
    if ell.symmetric
        function new_piece()
            if space == PrimalSpace
                if ell.superset !== nothing
                    Q = Symmetric(new_Q() + ell.superset.Q)
                else
                    Q = new_Q()
                end
                return Sets.Ellipsoid(Q)
            else
                ell.superset === nothing || error("superset not supported in dual space")
                @assert space == DualSpace
                return Sets.Ellipsoid(new_Q())
            end
        end
        if ell.piecewise === nothing
            set = new_piece()
        else
            hashyperplanes(ell.piecewise) && error("hyperplanes not supported for piecewise")
            sets = [new_piece() for i in 1:nhalfspaces(ell.piecewise)]
            set = Sets.Piecewise(sets, ell.piecewise)
            sos_extension().piecewise_ellipsoid_continuity(model, set, n)
            for i in eachindex(set.graph)
                for (j, v) in set.graph[i]
                    if i < j # The constraints are the same for (i, j) and (j, i)
                        # v corresponds to `-n_ij` in LCSS paper
                        Δ = sets[i].Q * v - sets[j].Q * v
                        inter = set.pieces[i] ∩ set.pieces[j]
                        h = HalfSpace(Δ, zero(eltype(Δ)))
                        @constraint(model, inter ⊆ h)
                    end
                end
            end
        end
        if space == PrimalSpace
            return set
        else
            return Polyhedra.polar(set)
        end
    else
        ell.superset === nothing || error("superset not supported for non-symmetric Ellipsoid")
        ell.piecewise === nothing || error("piecewise not supported for non-symmetric Ellipsoid")
        if space == PrimalSpace
            error("Non-symmetric ellipsoid non implemented yet, use `Ellipsoid(symmetric=true)`.")
        else
            @assert space == DualSpace
            if ell.point === nothing
                throw(ArgumentError("Specify a point for nonsymmetric ellipsoid, e.g. `Ellipsoid(point=InteriorPoint([1.0, 0.0]))"))
            end
            return polar_perspective_ellipsoid(model, new_Q(), ell.point,
                                               data(model).perspective_polyvar,
                                               space_polyvars)
        end
    end
end
function JuMP.value(ell::Sets.Ellipsoid)
    return Sets.Ellipsoid(Symmetric(JuMP.value.(ell.Q)))
end

### PolySet ###
struct PolySet <: AbstractVariable
    point::Union{Nothing, HintPoint}
    symmetric::Bool
    degree::Int
    dimension::Union{Nothing, Int}
    convex::Bool
    variables::Union{Nothing, AbstractVector}
    superset::Union{Sets.PolySet, Nothing}
    basis::Type
end
PolySet(; kws...) = sos_extension().polyset(; kws...)

Sets.space_variables(p::PolySet) = p.variables
_value(convexity_proof::Nothing) = nothing
function JuMP.value(set::Sets.PolySet)
    return Sets.PolySet(set.degree, JuMP.value(set.p))
end
function JuMP.value(set::Sets.ConvexPolySet)
    return Sets.ConvexPolySet(set.degree, JuMP.value(set.p), _value(set.convexity_proof))
end
function JuMP.value(set::Sets.Polar)
    return Polyhedra.polar(JuMP.value(Polyhedra.polar(set)))
end
function JuMP.value(set::Sets.PerspectiveDual)
    return Sets.perspective_dual(JuMP.value(Sets.perspective_dual(set)))
end
function JuMP.value(set::Sets.Householder)
    return Sets.Householder(JuMP.value(set.set), JuMP.value(set.p), set.h,
                            set.z, set.x)
end
function JuMP.value(set::Sets.ShiftedEllipsoid)
    return Sets.ShiftedEllipsoid(Symmetric(JuMP.value.(set.Q)),
                                 JuMP.value.(set.b), JuMP.value(set.β))
end
function JuMP.value(set::Sets.ConvexPolynomialSet)
    return Sets.ConvexPolynomialSet(set.degree, JuMP.value(set.q), set.z, set.x)
end
function JuMP.value(set::Sets.Piecewise)
    return Sets.Piecewise(JuMP.value.(set.sets), set.polytope, set.pieces, set.graph)
end

### SetVariableRef ###
mutable struct SetVariableRef{M <: JuMP.AbstractModel,
                              S <: AbstractVariable} <: JuMP.AbstractVariableRef
    model::M
    set::S
    name::String
    # `variable` is typically `Sets.AbstractSet{JuMP.VariableRef}` but it can
    # also be `Sets.AbstractSet{JuMP.AffExpr}` with `PolySet(superset = ...)`.
    variable::Union{Nothing, Sets.AbstractSet}
    space_index::Union{Nothing, SpaceIndex}
end
JuMP.name(vref::SetVariableRef) = vref.name
JuMP.is_valid(model::JuMP.Model, vref::SetVariableRef) = model === vref.model
function JuMP.build_variable(_error::Function, info::JuMP.VariableInfo, set::AbstractVariable)
    @assert !info.has_lb && !info.has_ub && !info.has_fix && !info.binary && !info.integer && !info.has_start
    return set
end
function JuMP.add_variable(model::JuMP.AbstractModel, set::AbstractVariable, name::String)
    vref = SetVariableRef(model, set, name, nothing, nothing)
    d = data(model)
    @assert d.state == Modeling
    push!(d.variables, vref)
    return vref
end
JuMP.value(vref::SetVariableRef) = JuMP.value(vref.variable)

function clear_spaces(vref::SetVariableRef)
    vref.space_index = nothing
end
Sets.space_variables(::SetVariableRef) = nothing
Sets.perspective_variable(::SetVariableRef) = nothing
function create_spaces(vref::SetVariableRef, spaces::Spaces)
    if vref.space_index === nothing
        if Sets.space_variables(vref.set) === nothing
            if vref.set.dimension === nothing
                vref.space_index = new_space(spaces)
            else
                vref.space_index = new_space(spaces, vref.set.dimension)
            end
        else
            vref.space_index = new_space(spaces, Sets.space_variables(vref.set))
        end
    end
    return vref.space_index
end
space_index(vref::SetVariableRef) = vref.space_index

function load(model::JuMP.AbstractModel, vref::SetVariableRef)
    d = data(model)
    vref.variable = variable_set(model, vref.set, d.space,
                                 space_dimension(d.spaces, vref.space_index),
                                 variable_polyvars(vref.set, d.spaces, vref.space_index))
end

variablify(p::Sets.Projection) = Polyhedra.project(variablify(p.set), p.indices)

variable_polyvars(::Polytope, ::Spaces, ::SpaceIndex) = nothing
function variable_polyvars(ell::Ellipsoid, spaces::Spaces, index::SpaceIndex)
    return ell.symmetric ? nothing : space_polyvars(spaces, index)
end
variable_polyvars(::PolySet, spaces::Spaces, index::SpaceIndex) = space_polyvars(spaces, index)

JuMP.value(set::Sets.Intersection) = Sets.Intersection(JuMP.value(set.set), set.domain)
