abstract type AbstractVariable <: JuMP.AbstractVariable end
struct Ellipsoid <: AbstractVariable
    dimension::Int
end
function variable_set(model::JuMP.AbstractModel, ell::Ellipsoid, space::Space)
    n = ell.dimension
    Q = @variable(model, [1:n, 1:n], Symmetric)
    if space == PrimalSpace
        return Sets.EllipsoidAtOrigin(Q)
    else
        @assert space == DualSpace
        return Sets.PolarEllipsoidAtOrigin(Q)
    end
end
function JuMP.result_value(ell::Sets.EllipsoidAtOrigin)
    return Sets.EllipsoidAtOrigin(Symmetric(JuMP.result_value.(ell.Q)))
end
function JuMP.result_value(ell::Sets.PolarEllipsoidAtOrigin)
    return Sets.PolarEllipsoidAtOrigin(Symmetric(JuMP.result_value.(ell.Q)))
end

mutable struct VariableRef{M <: JuMP.AbstractModel,
                           S <: AbstractVariable} <: JuMP.AbstractVariableRef
    model::M
    set::S
    name::String
    variable::Union{Nothing, Sets.AbstractSet{JuMP.VariableRef}}
end
JuMP.name(vref::VariableRef) = vref.name
function JuMP.build_variable(_error, info::JuMP.VariableInfo, ell::Ellipsoid)
    @assert !info.has_lb && !info.has_ub && !info.has_fix && !info.binary && !info.integer && !info.has_start
    return ell
end
function JuMP.add_variable(model::JuMP.AbstractModel, ell::Ellipsoid, name::String)
    vref = VariableRef(model, ell, name, nothing)
    push!(data(model).variables, vref)
    return vref
end
function load(model::JuMP.AbstractModel, ell::VariableRef)
    d = data(model)
    ell.variable = variable_set(model, ell.set, d.space)
end
JuMP.result_value(vref::VariableRef) = JuMP.result_value(vref.variable)
