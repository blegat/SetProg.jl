abstract type AbstractVariable <: JuMP.AbstractVariable end
struct Ellipsoid <: AbstractVariable
end
struct VariableRef{M <: JuMP.AbstractModel,
                   S <: AbstractVariable} <: JuMP.AbstractVariableRef
    model::M
    set::S
    name::String
end
JuMP.name(vref::VariableRef) = vref.name
function JuMP.build_variable(_error, info::JuMP.VariableInfo, ell::Ellipsoid)
    @assert !info.has_lb && !info.has_ub && !info.has_fix && !info.binary && !info.integer && !info.has_start
    return ell
end
function JuMP.add_variable(model::JuMP.AbstractModel, ell::Ellipsoid, name::String)
    return VariableRef(model, ell, name)
end
