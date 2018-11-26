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
    dimension::Int
    guaranteed_psd::Bool # Is it already guaranteed that it is PSD ? e.g. by nth_root
end
function Ellipsoid(; point::Union{Nothing, HintPoint}=nothing,
                   symmetric::Bool=false,
                   dimension::Union{Int, Nothing}=nothing)
    if dimension === nothing
        if point !== nothing
            dimension = length(point.h)
        else
            error("Dimension of Ellipsoid not specified, use Ellipsoid(dimension=...)")
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
function dual_quad_cone(model, Q::Symmetric{JuMP.VariableRef}, point::HintPoint)
    n = LinearAlgebra.checksquare(Q)
    d = data(model)
    y = [d.perspective_polyvar; d.polyvars[1:n]]
    return dual_quad_cone(model, Q, point, y)
end
function variable_set(model::JuMP.AbstractModel, ell::Ellipsoid, space::Space)
    n = ell.dimension
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
            return dual_quad_cone(model, Q, ell.point)
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
    dimension::Int
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
        error("Dimension of Ellipsoid not specified, use Ellipsoid(dimension=..., ...)")
    end
    return PolySet(point, symmetric, degree, dimension, convex)
end

function constrain_convex(model, p, vars)
    hessian = differentiate(p, vars, 2)
    return @constraint(model, hessian in SOSMatrixCone())
end

function variable_set(model::JuMP.AbstractModel, set::PolySet, space::Space)
    n = set.dimension
    vars = data(model).polyvars[1:n]
    # General all monomials of degree `degree`, we don't want monomials of
    # lower degree as the polynomial is homogeneous
    @assert iseven(set.degree)
    if set.convex
        monos = monomials(vars, div(set.degree, 2))
        p = @variable(model, variable_type=SOSPoly(monos))
        cref = constrain_convex(model, p, vars)
        slack = SumOfSquares.PolyJuMP.getdelegate(cref).slack
        convexity_proof = MultivariateMoments.getmat(slack)
        if space == PrimalSpace
            return Sets.ConvexPolynomialSublevelSetAtOrigin(set.degree, p, convexity_proof)
        else
            @assert space == DualSpace
            return Sets.PolarConvexPolynomialSublevelSetAtOrigin(set.degree, p, convexity_proof)
        end
    else
        error("TODO")
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

### VariableRef ###
mutable struct VariableRef{M <: JuMP.AbstractModel,
                           S <: AbstractVariable} <: JuMP.AbstractVariableRef
    model::M
    set::S
    name::String
    variable::Union{Nothing, Sets.AbstractSet{JuMP.VariableRef}}
end
JuMP.name(vref::VariableRef) = vref.name
function JuMP.build_variable(_error, info::JuMP.VariableInfo, set::AbstractVariable)
    @assert !info.has_lb && !info.has_ub && !info.has_fix && !info.binary && !info.integer && !info.has_start
    return set
end
function JuMP.add_variable(model::JuMP.AbstractModel, set::AbstractVariable, name::String)
    vref = VariableRef(model, set, name, nothing)
    d = data(model)
    @assert d.state == Modeling
    push!(d.variables, vref)
    return vref
end
function load(model::JuMP.AbstractModel, vref::VariableRef)
    d = data(model)
    variable = variable_set(model, vref.set, d.space)
    vref.variable = variable
end
JuMP.value(vref::VariableRef) = JuMP.value(vref.variable)
