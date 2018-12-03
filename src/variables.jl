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


### Ellipsoid ###
struct Ellipsoid <: AbstractVariable
    point::Union{Nothing, HintPoint}
    symmetric::Bool
    dimension::Union{Nothing, Int}
    guaranteed_psd::Bool # Is it already guaranteed that it is PSD ? e.g. by nth_root
end
function Ellipsoid(; point::Union{Nothing, HintPoint}=nothing,
                   symmetric::Bool=false,
                   dimension::Union{Int, Nothing}=nothing)
    if dimension === nothing
        if point !== nothing
            dimension = length(point.h)
        end
    end
    return Ellipsoid(point, symmetric, dimension, false)
end
function dual_quad_cone(model, Q::Symmetric{JuMP.VariableRef},
                        point::CenterPoint, y::Vector)
    @constraint(model, Q in PSDCone())
    return Sets.CenterDualQuadCone(Q, y, point.h)
end
function dual_quad_cone(model, Q::Symmetric{JuMP.VariableRef},
                        point::InteriorPoint, y::Vector)
    n = LinearAlgebra.checksquare(Q)
    @assert n == length(point.h)
    β = @variable(model, base_name="β")
    b = @variable(model, [1:length(point.h)], base_name="b")
    @constraint(model, Symmetric([β+1 b'; b Q]) in PSDCone())
    return Sets.InteriorDualQuadCone(Q, b, β, y, point.h)
end
function variable_set(model::JuMP.AbstractModel, ell::Ellipsoid, space::Space,
                      space_dimension, space_polyvars)
    n = space_dimension
    Q = @variable(model, [1:n, 1:n], Symmetric, base_name="Q")
    if !ell.guaranteed_psd
        @constraint(model, Q in PSDCone())
    end
    if ell.symmetric
        if space == PrimalSpace
            return Sets.EllipsoidAtOrigin(Q)
        else
            @assert space == DualSpace
            return Sets.PolarEllipsoidAtOrigin(Q)
        end
    else
        if space == PrimalSpace
            error("TODO")
        else
            @assert space == DualSpace
            if ell.point === nothing
                throw(ArgumentError("Specify a point for nonsymmetric ellipsoid, e.g. `Ellipsoid(point=InteriorPoint([1.0, 0.0]))"))
            end
            return dual_quad_cone(model, Q, ell.point,
                                  [data(model).perspective_polyvar;
                                   space_polyvars])
        end
    end
end
function JuMP.value(ell::Sets.EllipsoidAtOrigin)
    return Sets.EllipsoidAtOrigin(Symmetric(JuMP.value.(ell.Q)))
end
function JuMP.value(ell::Sets.PolarEllipsoidAtOrigin)
    return Sets.PolarEllipsoidAtOrigin(Symmetric(JuMP.value.(ell.Q)))
end

### PolySet ###
struct PolySet <: AbstractVariable
    point::Union{Nothing, HintPoint}
    symmetric::Bool
    degree::Int
    dimension::Union{Nothing, Int}
    convex::Bool
end
function PolySet(; point::Union{Nothing, HintPoint}=nothing,
                 symmetric::Bool=false,
                 degree::Union{Int, Nothing}=nothing,
                 dimension::Union{Int, Nothing}=nothing,
                 convex::Bool=false)
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
    return PolySet(point, symmetric, degree, dimension, convex)
end

function constrain_convex(model, p, vars)
    hessian = differentiate(p, vars, 2)
    return @constraint(model, hessian in SOSMatrixCone())
end

function variable_set(model::JuMP.AbstractModel, set::PolySet, space::Space,
                      space_dimension, space_polyvars)
    n = space_dimension
    d = data(model)
    # General all monomials of degree `degree`, we don't want monomials of
    # lower degree as the polynomial is homogeneous
    @assert iseven(set.degree)
    if set.convex
        if set.symmetric
            monos = monomials(space_polyvars, div(set.degree, 2))
            p = @variable(model, variable_type=SOSPoly(monos))
            cref = constrain_convex(model, p, space_polyvars)
            slack = SumOfSquares.PolyJuMP.getdelegate(cref).slack
            convexity_proof = MultivariateMoments.getmat(slack)
            if space == PrimalSpace
                return Sets.ConvexPolynomialSublevelSetAtOrigin(set.degree, p, convexity_proof)
            else
                @assert space == DualSpace
                return Sets.PolarConvexPolynomialSublevelSetAtOrigin(set.degree, p, convexity_proof)
            end
        else
            monos = monomials(lift_space_variables(d, space_polyvars),
                              div(set.degree, 2))
            p = @variable(model, variable_type=SOSPoly(monos))
            cref = constrain_convex(model, subs(p, d.perspective_polyvar => 1),
                                    space_polyvars)
            if space == PrimalSpace
                error("Non-symmetric PolySet in PrimalSpace not implemented yet")
            else
                @assert space == DualSpace
                if set.point === nothing
                    throw(ArgumentError("Specify a point for nonsymmetric polyset, e.g. `PolySet(point=InteriorPoint([1.0, 0.0]))"))
                end
                return Sets.DualConvexPolynomialCone(set.degree, p, set.point.h,
                                                     d.perspective_polyvar,
                                                     space_polyvars)
            end
        end
    else
        error("Non-convex PolySet not implemented yet")
    end
end
_value(convexity_proof::Nothing) = nothing
function _value(convexity_proof::MultivariateMoments.SymMatrix)
    return MultivariateMoments.SymMatrix(JuMP.value.(convexity_proof.Q),
                                         convexity_proof.n)
end
function JuMP.value(set::Sets.ConvexPolynomialSublevelSetAtOrigin)
    return Sets.ConvexPolynomialSublevelSetAtOrigin(set.degree, JuMP.value(set.p), _value(set.convexity_proof))
end
function JuMP.value(set::Sets.PolarConvexPolynomialSublevelSetAtOrigin)
    return Sets.PolarConvexPolynomialSublevelSetAtOrigin(set.degree, JuMP.value(set.p), _value(set.convexity_proof))
end
function JuMP.value(set::Sets.CenterDualQuadCone)
    return Sets.CenterDualQuadCone(JuMP.value(set.p),
                                   Symmetric(JuMP.value.(set.Q)), set.h, set.H)
end
function JuMP.value(set::Sets.InteriorDualQuadCone)
    return Sets.InteriorDualQuadCone(JuMP.value(set.p),
                                     Symmetric(JuMP.value.(set.Q)),
                                     JuMP.value.(set.b),
                                     JuMP.value(set.β), set.h, set.H)
end
function JuMP.value(set::Sets.DualConvexPolynomialCone)
    return Sets.DualConvexPolynomialCone(set.degree, JuMP.value(set.q),
                                         JuMP.value(set.p), set.h, set.H, set.z,
                                         set.x)
end

### VariableRef ###
mutable struct VariableRef{M <: JuMP.AbstractModel,
                           S <: AbstractVariable} <: JuMP.AbstractVariableRef
    model::M
    set::S
    name::String
    variable::Union{Nothing, Sets.AbstractSet{JuMP.VariableRef}}
    space_index::Union{Nothing, SpaceIndex}
end
JuMP.name(vref::VariableRef) = vref.name
function JuMP.build_variable(_error, info::JuMP.VariableInfo, set::AbstractVariable)
    @assert !info.has_lb && !info.has_ub && !info.has_fix && !info.binary && !info.integer && !info.has_start
    return set
end
function JuMP.add_variable(model::JuMP.AbstractModel, set::AbstractVariable, name::String)
    vref = VariableRef(model, set, name, nothing, nothing)
    d = data(model)
    @assert d.state == Modeling
    push!(d.variables, vref)
    return vref
end
JuMP.value(vref::VariableRef) = JuMP.value(vref.variable)

function clear_spaces(vref::VariableRef)
    vref.space_index = nothing
end
function create_spaces(vref::VariableRef, spaces::Spaces)
    if vref.space_index === nothing
        if vref.set.dimension === nothing
            vref.space_index = new_space(spaces)
        else
            vref.space_index = new_space(spaces, vref.set.dimension)
        end
    end
    return vref.space_index
end
space_index(vref::VariableRef) = vref.space_index

function load(model::JuMP.AbstractModel, vref::VariableRef)
    d = data(model)
    vref.variable = variable_set(model, vref.set, d.space,
                                 space_dimension(d.spaces, vref.space_index),
                                 space_polyvars(d.spaces, vref.space_index))
end
