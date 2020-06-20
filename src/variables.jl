abstract type AbstractVariable <: JuMP.AbstractVariable end

abstract type HintPoint end

# It will not really be the center and the center for z = 0 is the the same as for <h, x> = 0
struct CenterPoint{T} <: HintPoint
    h::Vector{T}
end

struct InteriorPoint{T} <: HintPoint
    h::Vector{T}
end
_β(model, h::InteriorPoint) = @variable(model)
_b(model, h::InteriorPoint) = @variable(model, [1:length(h.h)])

function polar_perspective_ellipsoid(ell, point::HintPoint, z::SpaceVariable,
                                     x::Vector{SpaceVariable})
    y = [z; x]
    H = Sets._householder(point.h)
    p = y' * H * Sets._HPH(ell) * H * y
    return Sets.perspective_dual(Sets.Householder(ell, p, point.h, z, x))
end
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

function perspective_dual_polyset(set, point::HintPoint, z::SpaceVariable,
                                  x::Vector{SpaceVariable})
    y = [z; x]
    H = Sets._householder(point.h)
    p = Sets.perspective_gauge0(set)(y => H * y)
    return Sets.perspective_dual(Sets.Householder(set, p, point.h, z, x))
end
# For `CenterPoint`, `q` should be non-perspective, need refactoring
function perspective_dual_polyset(degree, q, point::InteriorPoint, z, x)
    set = Sets.ConvexPolynomialSet(degree, q, z, x)
    perspective_dual_polyset(set, point, z, x)
end

### Ellipsoid ###
struct Ellipsoid <: AbstractVariable
    point::Union{Nothing, HintPoint}
    symmetric::Bool
    dimension::Union{Nothing, Int}
    guaranteed_psd::Bool # Is it already guaranteed that it is PSD ? e.g. by nth_root
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
            @polyvar x[1:n]
            q = [quad_form(set.Q, x) for set in sets]
            for i in eachindex(set.graph)
                for (j, v) in set.graph[i]
                    if i < j # The constraints are the same for (i, j) and (j, i)
                        @constraint(model, q[i] == q[j], domain = @set x'v == 0)
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
            return Sets.polar(set)
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
    variables::Union{Nothing, Vector{SpaceVariable}}
    superset::Union{Sets.PolySet, Nothing}
    basis::Type
end
function PolySet(; point::Union{Nothing, HintPoint}=nothing,
                 symmetric::Bool=false,
                 degree::Union{Int, Nothing}=nothing,
                 dimension::Union{Int, Nothing}=nothing,
                 convex::Bool=false,
                 variables::Union{Vector{SpaceVariable}, Nothing}=nothing,
                 superset::Union{Sets.PolySet, Nothing}=nothing,
                 basis::Type=MultivariateBases.MonomialBasis)
    if degree === nothing
        error("Degree of PolySet not specified, use PolySet(degree=..., ...)")
    end
    if isodd(degree)
        throw(ArgumentError("Degree of PolySet not even"))
    end
    if dimension === nothing
        if point !== nothing
            dimension = length(point.h)
        end
    end
    if superset !== nothing
        if dimension === nothing
            dimension = Sets.dimension(superset)
        elseif dimension != Sets.dimension(superset)
            throw(DimensionMismatch())
        end
        if variables === nothing
            variables = Sets.space_variables(superset)
        elseif variables != Sets.space_variables(superset)
            error("Space variables set does not correspond to superset space variables.")
        end
    end
    return PolySet(point, symmetric, degree, dimension, convex, variables, superset, basis)
end
Sets.space_variables(p::PolySet) = p.variables

function constrain_convex(model, p, vars)
    hessian = differentiate(p, vars, 2)
    # We do not just do `@constraint(model, p in SOSConvex())` as we would
    # like to have access to the PSD matrix of variables for the det volume heuristic
    y = [MultivariatePolynomials.similarvariable(eltype(hessian), gensym()) for i in 1:LinearAlgebra.checksquare(hessian)]
    q = dot(y, hessian * y)
    X = SumOfSquares.Certificate.monomials_half_newton_polytope(MultivariatePolynomials.monomials(q), (y,))
    # If `X` is empty, we will need the following bridge
    JuMP.add_bridge(model, SumOfSquares.Bridges.Constraint.EmptyBridge)
    # If `length(X)` is 2, we will need the following bridge
    JuMP.add_bridge(model, SumOfSquares.Bridges.Constraint.PositiveSemidefinite2x2Bridge)
    set = SumOfSquares.matrix_cone(MOI.PositiveSemidefiniteConeTriangle,
                                   length(X))
    Q = @variable(model, [1:MOI.dimension(set)])
    @constraint(model, Q in set)
    s = SumOfSquares.build_gram_matrix(Q, MonomialBasis(X))
    @constraint(model, q == s)
    return MultivariateMoments.getmat(s)
end

function variable_set(model::JuMP.AbstractModel, set::PolySet, space::Space,
                      space_dimension, space_polyvars)
    n = space_dimension
    d = data(model)
    # General all monomials of degree `degree`, we don't want monomials of
    # lower degree as the polynomial is homogeneous
    @assert iseven(set.degree)
    if set.symmetric
        monos = monomials(space_polyvars, div(set.degree, 2))
    else
        monos = monomials(lift_space_variables(d, space_polyvars),
                          div(set.degree, 2))
    end
    basis = MultivariateBases.basis_covering_monomials(set.basis, monos)
    # TODO If `set.convex` and `set.symmetric`, no need for the poly to be SOS, see Lemma 6.33 of [BPT12]
    p = @variable(model, variable_type=SOSPoly(basis))
    if set.convex
        set.superset === nothing || error("superset not supported for convex PolySet")
        if set.symmetric
            convexity_proof = constrain_convex(model, p, space_polyvars)
            if space == PrimalSpace
                return Sets.ConvexPolySet(set.degree, p, convexity_proof)
            else
                @assert space == DualSpace
                return Sets.polar(Sets.ConvexPolySet(set.degree, p, convexity_proof))
            end
        else
            constrain_convex(model, subs(p, d.perspective_polyvar => 1),
                             space_polyvars)
            if space == PrimalSpace
                error("Non-symmetric PolySet in PrimalSpace not implemented yet")
            else
                @assert space == DualSpace
                if set.point === nothing
                    throw(ArgumentError("Specify a point for nonsymmetric polyset, e.g. `PolySet(point=InteriorPoint([1.0, 0.0]))"))
                end
                return perspective_dual_polyset(set.degree, p, set.point, d.perspective_polyvar, space_polyvars)
            end
        end
    else
        if set.symmetric
            if space == PrimalSpace
                if set.superset !== nothing
                    p = SetProg.SumOfSquares.gram_operate(+, set.superset.p, p)
                end
                return Sets.PolySet(set.degree, p)
            else
                error("Non-convex PolySet not supported in $space")
            end
        else
            error("Non-convex nonsymmetric PolySet not implemented yet.")
        end
    end
end
_value(convexity_proof::Nothing) = nothing
function _value(convexity_proof::MultivariateMoments.SymMatrix)
    return MultivariateMoments.SymMatrix(JuMP.value.(convexity_proof.Q),
                                         convexity_proof.n)
end
function JuMP.value(set::Sets.PolySet)
    return Sets.PolySet(set.degree, JuMP.value(set.p))
end
function JuMP.value(set::Sets.ConvexPolySet)
    return Sets.ConvexPolySet(set.degree, JuMP.value(set.p), _value(set.convexity_proof))
end
function JuMP.value(set::Sets.Polar)
    return Sets.polar(JuMP.value(Sets.polar(set)))
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
function JuMP.build_variable(_error, info::JuMP.VariableInfo, set::AbstractVariable)
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
function Sets.perspective_variable(::SetVariableRef) end
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
                                 space_polyvars(d.spaces, vref.space_index))
end
