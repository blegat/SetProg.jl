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

function JuMP.set_objective(model::JuMP.Model,
                            sense::MOI.OptimizationSense,
                            func::AbstractScalarFunction)
    data(model).objective_sense = sense
    data(model).objective = func
end

# Primal:
#   set : x^T Q x ≤ 1
#   volume proportional to 1/det(Q)
#   t ≤ det(Q)^(1/n) <=> 1/t^n ≥ 1/det(Q)
#   volume proportional to 1/t^n
#   For t ≤ det(Q)^(1/n) to be tight we need to maximize `t`
#   hence we need to minimize the volume
# Dual:
#   set : x^T Q^{-1} x ≤ 1
#   volume proportional to det(Q)
#   t ≤ det(Q)^(1/n) <=> t^n ≥ det(Q)
#   volume proportional to t^n
#   For t ≤ det(Q)^(1/n) to be tight we need to maximize `t`
#   hence we need to maximize the volume
function set_space(space::Space, ::RootVolume, model::JuMP.Model)
    sense = data(model).objective_sense
    if sense == MOI.MinSense
        return set_space(space, PrimalSpace)
    else
        # The sense cannot be FeasibilitySense since the objective function is
        # not nothing
        @assert sense == MOI.MaxSense
        return set_space(space, DualSpace)
    end
end
function root_volume(model::JuMP.Model, ell::Union{Sets.EllipsoidAtOrigin,
                                                   Sets.PolarEllipsoidAtOrigin})
    Q = ell.Q
    n = Sets.dimension(ell)
    t = @variable(model)
    upper_tri = [Q[i, j] for j in 1:n for i in 1:j]
    @constraint(model, [t; upper_tri] in MOI.RootDetConeTriangle(n))
    return t
end
function load(model::JuMP.Model, rv::RootVolume)
    t = root_volume(model, rv.variable.variable)
    JuMP.set_objective_sense(model, MOI.MaxSense)
    JuMP.set_objective_function(model, t)
end
