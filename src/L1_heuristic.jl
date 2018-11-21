struct L1Heuristic{V <: SetProg.VariableRef} <: AbstractScalarFunction
    variable::V
    rectangle_vertex::Vector{Float64}
end
Base.copy(l::L1Heuristic) = l
L1_heuristic(volume::Volume, v::Vector{Float64}) = RootVolume(volume.variable, v)
Base.show(io::IO, l::L1Heuristic) = print(io, "L1-heuristic(", l.variable, ")")

set_space(space::Space, ::L1Heuristic, ::JuMP.Model) = return space

function load(model::JuMP.Model, l::RootVolume)
    t = l1_integral(l.variable.variable, l.rectangle_vertex)
    JuMP.set_objective_sense(model, data(model).objective_sense)
    JuMP.set_objective_function(model, t)
end
