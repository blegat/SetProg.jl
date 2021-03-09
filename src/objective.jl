abstract type AbstractScalarFunction <: JuMP.AbstractJuMPScalar end

struct Volume{S} <: AbstractScalarFunction
    set::S
end
Polyhedra.volume(variable::Union{SetVariableRef, Sets.Projection{<:SetVariableRef}, Sets.AbstractSet{<:JuMP.AbstractJuMPScalar}}) = Volume(variable)

function JuMP.set_objective(model::JuMP.Model,
                            sense::MOI.OptimizationSense,
                            func::AbstractScalarFunction)
    d = data(model)
    @assert d.state == Modeling
    d.objective_sense = sense
    d.objective = func
end

function load(model::JuMP.Model, f::AbstractScalarFunction)
    JuMP.set_objective_sense(model, objective_sense(model, f))
    JuMP.set_objective_function(model, objective_function(model, f))
end

include("affine.jl")
include("nth_root.jl")
include("L1_heuristic.jl")
