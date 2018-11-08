abstract type AbstractScalarFunction <: JuMP.AbstractJuMPScalar end

struct Volume{V <: SetProg.VariableRef} <: AbstractScalarFunction
    variable::V
end
Polyhedra.volume(variable::SetProg.VariableRef) = Volume(variable)
struct RootVolume{V <: SetProg.VariableRef} <: AbstractScalarFunction
    variable::V
end
Base.copy(rv::RootVolume) = rv
nth_root(volume::Volume) = RootVolume(volume.variable)
Base.show(io::IO, rv::RootVolume) = print(io, "volume^(1/n)(", rv.variable, ")")

function JuMP.set_objective_function(model::JuMP.AbstractModel,
                                     func::AbstractScalarFunction)
    data(model).objective = func
end
