abstract type SetConstraint  <: JuMP.AbstractConstraint end

function load(model::JuMP.Model, constraint::SetConstraint)
    load(model, variablify(constraint))
end

variablify(p::Polyhedra.Rep) = p
variablify(v::VariableRef) = v.variable

struct InclusionConstraint{SubSetType, SupSetType} <: SetConstraint
    subset::SubSetType
    supset::SupSetType
end
JuMP.function_string(print_mode, c::InclusionConstraint) = string(c.subset)
JuMP.in_set_string(print_mode, c::InclusionConstraint) = string("⊆ ", c.supset)
function variablify(c::InclusionConstraint)
    return InclusionConstraint(variablify(c.subset), variablify(c.supset))
end

# Primal:
#   set : x^T Q x ≤ 1
#   volume proportional to 1/det(Q)
#   t ≤ det(Q)^(1/n) <=> 1/t^n ≥ 1/det(Q)
#   volume proportional to 1/t^n
#   For t ≤ det(Q)^(1/n) to be tight we need to maximize `t`
#   hence we need to minimize the volume
function set_space(space::Space, ::InclusionConstraint{<:VariableRef, <:Polyhedra.Rep})
    return set_space(space, DualSpace)
end
function set_space(space::Space, ::InclusionConstraint{<:Polyhedra.Rep, <:VariableRef})
    return set_space(space, PrimalSpace)
end
function load(model::JuMP.Model,
              constraint::InclusionConstraint{<:Sets.PolarEllipsoidAtOrigin,
                                              <:Polyhedra.Rep})
    Q = constraint.subset.Q
    p = constraint.supset
    for h in hyperplanes(p)
        @assert iszero(h.β) # Otherwise it is not symmetric around the origin
        @constraint(model, dot(h.a, Q * h.a) == 0.0)
    end
    for h in halfspaces(p)
        @constraint(model, dot(h.a, Q * h.a) ≤ h.β^2)
    end
end

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
const SetConstraintRef{M} = JuMP.ConstraintRef{M, ConstraintIndex, SetShape}

function JuMP.add_constraint(model::JuMP.Model, constraint::InclusionConstraint,
                             name::String)
    container = data(model)
    index = ConstraintIndex(container.last_index += 1)
    container.constraints[index] = constraint
    container.names[index] = name
    return JuMP.ConstraintRef(model, index, SetShape())
end
function JuMP.name(cref::SetConstraintRef)
    return data(cref.model).names[cref.index]
end
function JuMP.constraint_object(cref::SetConstraintRef)
    return data(cref.model).constraints[cref.index]
end

function load(model::JuMP.Model, cref::SetConstraintRef)
    load(model, data(model)[cref.index])
end
