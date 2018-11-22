struct RootVolume{V <: SetProg.VariableRef} <: AbstractScalarFunction
    variable::V
end
Base.copy(rv::RootVolume) = rv
nth_root(volume::Volume) = RootVolume(volume.variable)
Base.show(io::IO, rv::RootVolume) = print(io, "volume^(1/n)(", rv.variable, ")")

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

function ellipsoid_root_volume(model::JuMP.Model, Q::AbstractMatrix)
    n = LinearAlgebra.checksquare(Q)
    t = @variable(model)
    upper_tri = [Q[i, j] for j in 1:n for i in 1:j]
    @constraint(model, [t; upper_tri] in MOI.RootDetConeTriangle(n))
    return t
end

function root_volume(model::JuMP.Model, ell::Union{Sets.EllipsoidAtOrigin,
                                                   Sets.PolarEllipsoidAtOrigin,
                                                   Sets.DualQuadCone})
    return ellipsoid_root_volume(model, ell.Q)
end

"""
    root_volume(model::JuMP.Model, set::Union{Sets.PolynomialSublevelSetAtOrigin,
                                              Sets.PolarPolynomialSublevelSetAtOrigin})

Section IV.A of [MLB05].

[MLB05] A. Magnani, S. Lall and S. Boyd.
*Tractable fitting with convex polynomials via sum-of-squares*.
Proceedings of the 44th IEEE Conference on Decision and Control, and European Control Conference 2005,
**2005**.
"""
function root_volume(model::JuMP.Model, set::Union{Sets.PolynomialSublevelSetAtOrigin,
                                                   Sets.PolarPolynomialSublevelSetAtOrigin})
    if set.convexity_proof === nothing
        error("Cannot optimize volume of non-convex polynomial sublevel set. Use PolySet(convex=true, ...)")
    end
    return ellipsoid_root_volume(model, set.convexity_proof)
end

objective_sense(::JuMP.Model, ::RootVolume) = MOI.MaxSense

function load(model::JuMP.Model, rv::RootVolume)
    t = root_volume(model, rv.variable.variable)
    JuMP.set_objective_sense(model, MOI.MaxSense)
    JuMP.set_objective_function(model, t)
end
