abstract type SetConstraint  <: JuMP.AbstractConstraint end

## Loading  ##

need_variablify(p::Polyhedra.Rep) = false
variablify(p::Polyhedra.Rep) = p
need_variablify(v::VariableRef) = true
variablify(v::VariableRef) = v.variable

function load(model::JuMP.Model, constraint::SetConstraint)
    @assert need_variablify(constraint) # Avoid infinite recursion
    return JuMP.add_constraint(model, variablify(constraint))
end

## Adding ##

struct SetShape <: JuMP.AbstractShape end
struct ConstraintIndex
    value::Int
end

function JuMP.add_constraint(model::JuMP.Model, constraint::SetConstraint,
                             name::String)
    container = data(model)
    index = ConstraintIndex(container.last_index += 1)
    container.constraints[index] = constraint
    container.names[index] = name
    return JuMP.ConstraintRef(model, index, SetShape())
end

### SetConstraintRef ###

const SetConstraintRef{M} = JuMP.ConstraintRef{M, ConstraintIndex, SetShape}
function JuMP.name(cref::SetConstraintRef)
    return data(cref.model).names[cref.index]
end
function JuMP.constraint_object(cref::SetConstraintRef)
    return data(cref.model).constraints[cref.index]
end

include("membership.jl")
include("inclusion.jl")
