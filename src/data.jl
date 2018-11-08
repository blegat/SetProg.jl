mutable struct Data
    constraints::Dict{ConstraintIndex, InclusionConstraint}
    names::Dict{ConstraintIndex, String}
    last_index::Int
    objective::Union{Nothing, AbstractScalarFunction}
end
function data(model::JuMP.Model)
    if !haskey(model.ext, :SetProg)
        model.ext[:SetProg] = Data(Dict{ConstraintIndex, InclusionConstraint}(),
                                   Dict{ConstraintIndex, String}(), 0, nothing)
    end
    return model.ext[:SetProg]
end

