mutable struct Data
    variables::Set{VariableRef}
    constraints::Dict{ConstraintIndex, InclusionConstraint}
    names::Dict{ConstraintIndex, String}
    last_index::Int
    objective::Union{Nothing, AbstractScalarFunction}
    objective_variable::Union{Nothing, JuMP.VariableRef}
    space::Space
end

function set_space(cur::Space, space::Space)
    if space == Undecided || cur == Undecided || cur == space
        return space
    else
        error("Incompatible constraints/objective")
    end
end
function set_space(d::Data, model::JuMP.Model)
    space = Undecided
    for (index, constraint) in d.constraints
        space = set_space(space, constraint)
    end
    if d.objective !== nothing
        space = set_space(space, d.objective, model)
    end
    d.space = space
end

function data(model::JuMP.Model)
    if !haskey(model.ext, :SetProg)
        model.ext[:SetProg] = Data(Set{VariableRef}(),
                                   Dict{ConstraintIndex, InclusionConstraint}(),
                                   Dict{ConstraintIndex, String}(), 0, nothing,
                                   nothing, Undecided)
        model.optimize_hook = optimize_hook
    end
    return model.ext[:SetProg]
end

