# Avoid variables, constraints and objective being added when loading
@enum State Loading Modeling

mutable struct Data
    variables::Set{VariableRef}
    constraints::Dict{ConstraintIndex, SetConstraint}
    names::Dict{ConstraintIndex, String}
    last_index::Int
    objective_sense::MOI.OptimizationSense
    objective::Union{Nothing, AbstractScalarFunction}
    objective_variable::Union{Nothing, JuMP.VariableRef}
    polyvars::Union{Nothing, Vector{DynamicPolynomials.PolyVar{true}}}
    perspective_polyvar::DynamicPolynomials.PolyVar{true}
    space::Space
    state::State
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
        @polyvar z # perspective variable
        model.ext[:SetProg] = Data(Set{VariableRef}(),
                                   Dict{ConstraintIndex, SetConstraint}(),
                                   Dict{ConstraintIndex, String}(), 0,
                                   MOI.FeasibilitySense, nothing, nothing,
                                   nothing, z, Undecided, Modeling)
        model.optimize_hook = optimize_hook
    end
    return model.ext[:SetProg]
end
