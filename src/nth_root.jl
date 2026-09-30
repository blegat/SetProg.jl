struct RootVolume{S} <: AbstractScalarFunction
    set::S
end
Base.copy(rv::RootVolume) = rv
nth_root(volume::Volume) = RootVolume(volume.set)
Base.show(io::IO, rv::RootVolume) = print(io, "volume^(1/n)(", rv.set, ")")

struct LogVolume{S} <: AbstractScalarFunction
    set::S
end
Base.copy(lv::LogVolume) = lv
Base.log(volume::Volume) = LogVolume(volume.set)
Base.show(io::IO, lv::LogVolume) = print(io, "log(volume(", lv.set, "))")

const DetVolume = Union{RootVolume, LogVolume}

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
# The same reasoning holds for `t ≤ log(det(Q))`.
function set_space(space::Space, rv::DetVolume, model::JuMP.Model)
    if rv.set isa Ellipsoid
        rv.set.guaranteed_psd = true
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

_upper_tri(Q) = [Q[i, j] for j in 1:size(Q, 2) for i in 1:j]

function ellipsoid_root_volume(model::JuMP.Model, Q::AbstractMatrix)
    n = LinearAlgebra.checksquare(Q)
    t = @variable(model, base_name="t")
    @constraint(model, [t; _upper_tri(Q)] in MOI.RootDetConeTriangle(n))
    return t
end

# `(t, u, Q)` is in the `LogDetConeTriangle` iff `t ≤ u * log(det(Q / u))`
# so fixing `u = 1` gives `t ≤ log(det(Q))`.
function ellipsoid_log_volume(model::JuMP.Model, Q::AbstractMatrix)
    n = LinearAlgebra.checksquare(Q)
    t = @variable(model, base_name="t")
    @constraint(model, [t; 1; _upper_tri(Q)] in MOI.LogDetConeTriangle(n))
    return t
end

function volume_matrix(ell::Union{Sets.PolarOrNot{<:Sets.Ellipsoid},
                                  Sets.HouseDualOf{<:Sets.AbstractEllipsoid}})
    return Sets.convexity_proof(ell)
end

"""
    volume_matrix(set::Sets.PolarOrNot{<:Sets.ConvexPolySet})

Return the matrix whose determinant is used as volume heuristic, see
Section IV.A of [MLB05].

[MLB05] A. Magnani, S. Lall and S. Boyd.
*Tractable fitting with convex polynomials via sum-of-squares*.
Proceedings of the 44th IEEE Conference on Decision and Control, and European Control Conference 2005,
**2005**.
"""
function volume_matrix(set::Sets.PolarOrNot{<:Sets.ConvexPolySet})
    if Sets.convexity_proof(set) === nothing
        error("Cannot optimize volume of non-convex polynomial sublevel set.",
              " Use PolySet(convex=true, ...)")
    end
    return Sets.convexity_proof(set)
end

function root_volume(model::JuMP.Model, set)
    return ellipsoid_root_volume(model, volume_matrix(set))
end

function log_volume(model::JuMP.Model, set)
    return ellipsoid_log_volume(model, volume_matrix(set))
end

objective_sense(::JuMP.Model, ::DetVolume) = MOI.MAX_SENSE
function objective_function(model::JuMP.Model, rv::RootVolume)
    return root_volume(model, variablify(rv.set))
end
function objective_function(model::JuMP.Model, lv::LogVolume)
    return log_volume(model, variablify(lv.set))
end
