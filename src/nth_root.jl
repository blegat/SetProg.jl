struct RootVolume{V <: SetVariableRef} <: AbstractScalarFunction
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
function set_space(space::Space, rv::RootVolume, model::JuMP.Model)
    if rv.variable isa Ellipsoid
        rv.variable.guaranteed_psd = true
    end
    sense = data(model).objective_sense
    if sense == MOI.MIN_SENSE
        return set_space(space, PrimalSpace)
    else
        # The sense cannot be FEASIBILITY_SENSE since the objective function is
        # not nothing
        @assert sense == MOI.MAX_SENSE
        return set_space(space, DualSpace)
    end
end

function ellipsoid_root_volume(model::JuMP.Model, Q::AbstractMatrix)
    n = LinearAlgebra.checksquare(Q)
    t = @variable(model, base_name="t")
    upper_tri = [Q[i, j] for j in 1:n for i in 1:j]
    @constraint(model, [t; upper_tri] in MOI.RootDetConeTriangle(n))
    return t
end

function root_volume(model::JuMP.Model, ell::Union{Sets.PolarOrNot{<:Sets.EllipsoidAtOrigin},
                                                   Sets.HouseDualOf{<:Sets.AbstractEllipsoid}})
    return ellipsoid_root_volume(model, Sets.convexity_proof(ell))
end

"""
    root_volume(model::JuMP.Model,
                set::Sets.PolarOrNot{<:Sets.ConvexPolynomialSublevelSetAtOrigin})

Section IV.A of [MLB05].

[MLB05] A. Magnani, S. Lall and S. Boyd.
*Tractable fitting with convex polynomials via sum-of-squares*.
Proceedings of the 44th IEEE Conference on Decision and Control, and European Control Conference 2005,
**2005**.
"""
function root_volume(model::JuMP.Model,
                     set::Sets.PolarOrNot{<:Sets.ConvexPolynomialSublevelSetAtOrigin})
    if Sets.convexity_proof(set) === nothing
        error("Cannot optimize volume of non-convex polynomial sublevel set.",
              " Use PolySet(convex=true, ...)")
    end
    return ellipsoid_root_volume(model, Sets.convexity_proof(set))
end

objective_sense(::JuMP.Model, ::RootVolume) = MOI.MAX_SENSE
function objective_function(model::JuMP.Model, rv::RootVolume)
    return root_volume(model, rv.variable.variable)
end
