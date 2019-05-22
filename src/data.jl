# Avoid variables, constraints and objective being added when loading
@enum State Loading Modeling

mutable struct Data
    variables::Set{SetVariableRef}
    constraints::Dict{ConstraintIndex, SetConstraint}
    transformed_constraints::Dict{ConstraintIndex, JuMP.ConstraintRef}
    names::Dict{ConstraintIndex, String}
    last_index::Int
    objective_sense::MOI.OptimizationSense
    objective::Union{Nothing, AbstractScalarFunction}
    objective_variable::Union{Nothing, JuMP.VariableRef}
    perspective_polyvar::Union{Nothing, SpaceVariable}
    space::Space
    state::State
    spaces::Union{Nothing, Spaces}
end

function lift_space_variables(data::Data,
                              x::Vector{SpaceVariable})
    @assert data.perspective_polyvar > x[1]
    SpaceVariable[data.perspective_polyvar; x]
end

function set_space(cur::Space, space::Space)
    if space == Undecided || cur == Undecided || cur == space
        return space
    else
        error("Incompatible constraints/objective, some require to do the modeling in the primal space and some in the dual space.")
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
    if space == Undecided
        space = PrimalSpace
    end
    d.space = space
end

# Remove spaces creating in earlier `optimize!`
function clear_spaces(d::Data)
    for vref in d.variables
        clear_spaces(vref)
    end
    for (index, constraint) in d.constraints
        clear_spaces(constraint)
    end
end
# No space index stored in constants
function clear_spaces(::Union{Polyhedra.Rep, Sets.AbstractSet}) end

# No perspective variable
function Sets.perspective_variable(::Polyhedra.Rep) end

synchronize_perspective(::Nothing, ::Nothing) = nothing
synchronize_perspective(::Nothing, z::SpaceVariable) = z
synchronize_perspective(z::SpaceVariable, ::Nothing) = z
function synchronize_perspective(z1::SpaceVariable, z2::SpaceVariable)
    if z1 !== z2
        throw(ArgumentError("Perspective variables do not match"))
    end
    return z1
end
function synchronize_perspective(d::Data)
    z = nothing
    for (index, constraint) in d.constraints
        z = synchronize_perspective(z, Sets.perspective_variable(constraint))
    end
    if z === nothing
        @polyvar z
    end
    d.perspective_polyvar = z
end

function create_spaces(set::Polyhedra.Rep, spaces::Spaces)
    return new_space(spaces, Polyhedra.fulldim(set))
end
function create_spaces(set::Sets.AbstractSet, spaces::Spaces)
    polyvars = Sets.space_variables(set)
    if polyvars === nothing
        return new_space(spaces, Sets.dimension(set))
    else
        return new_space(spaces, polyvars)
    end
end

# Create spaces and identify objects lying the the same space
function create_spaces(d::Data)
    spaces = Spaces()
    for vref in d.variables
        create_spaces(vref, spaces)
    end
    for (index, constraint) in d.constraints
        create_spaces(constraint, spaces)
    end
    d.spaces = spaces
end

function data(model::JuMP.Model)
    if !haskey(model.ext, :SetProg)
        model.ext[:SetProg] = Data(Set{SetVariableRef}(),
                                   Dict{ConstraintIndex, SetConstraint}(),
                                   Dict{ConstraintIndex, JuMP.ConstraintRef}(),
                                   Dict{ConstraintIndex, String}(), 0,
                                   MOI.FEASIBILITY_SENSE, nothing, nothing,
                                   nothing, Undecided, Modeling, nothing)
        model.optimize_hook = optimize_hook
    end
    return model.ext[:SetProg]
end
