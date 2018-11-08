struct InclusionConstraint{SubSetType, SupSetType} <: JuMP.AbstractConstraint
    subset::SubSetType
    supset::SupSetType
end
JuMP.function_string(print_mode, c::InclusionConstraint) = string(c.subset)
JuMP.in_set_string(print_mode, c::InclusionConstraint) = string("⊆ ", c.supset)

struct PowerSet{S}
    set::S
end
function JuMP.build_constraint(_error::Function, subset,
                               supset_powerset::PowerSet)
    InclusionConstraint(subset, supset_powerset.set)
end

struct SetShape <: JuMP.AbstractShape end
struct ConstraintIndex
    value::Int
end
function JuMP.add_constraint(model::JuMP.Model, constraint::InclusionConstraint,
                             name::String)
    container = data(model)
    index = ConstraintIndex(container.last_index += 1)
    container.constraints[index] = constraint
    container.names[index] = name
    return ConstraintRef(model, index, SetShape())
end
function JuMP.name(cref::JuMP.ConstraintRef{JuMP.Model, ConstraintIndex})
    return data(cref.model).names[cref.index]
end
function JuMP.constraint_object(cref::JuMP.ConstraintRef{JuMP.Model,
                                                         ConstraintIndex})
    return data(cref.model).constraints[cref.index]
end
