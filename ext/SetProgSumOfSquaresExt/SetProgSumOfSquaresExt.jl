module SetProgSumOfSquaresExt

using LinearAlgebra
using JuMP
using Polyhedra
using SumOfSquares
using MultivariateMoments
import MultivariatePolynomials as MP
import MultivariateBases
import MultivariateBases as MB
import DynamicPolynomials
import SetProg
import SetProg:
    _unstable_constantify,
    data,
    lift_space_variables,
    Data,
    DualSpace,
    HintPoint,
    InclusionConstraint,
    InteriorPoint,
    MA,
    MOI,
    MembershipConstraint,
    PolySet,
    PowerSet,
    PrimalSpace,
    SetConstraintRef,
    Sets,
    Space,
    SpaceIndex,
    Spaces,
    _moment_matrix,
    _value,
    all_exponents,
    apply_map,
    apply_matrix,
    apply_transformation,
    constrain_convex,
    decompose,
    evaluate_monomials,
    integrate,
    integrate_decompositions,
    integrate_monomials,
    l1_integral,
    linear_algebraic_surface,
    perspective_dual_polyset,
    polar_perspective_ellipsoid,
    poly_eval,
    power_integrate,
    quad_form,
    rectangle_integrate,
    s_procedure,
    sublevel_eval,
    transformation,
    transformed_constraint,
    variable_set

function new_polyvars(n::Int)
    DynamicPolynomials.@polyvar x[1:n]
    return x
end
function new_perspective_variable()
    DynamicPolynomials.@polyvar z
    return z
end

function piecewise_ellipsoid_continuity(model, set, n)
    x = new_polyvars(n)
    q = [quad_form(piece.Q, x) for piece in set.sets]
    for i in eachindex(set.graph)
        for (j, v) in set.graph[i]
            if i < j
                @constraint(model, q[i] == q[j], domain = @set x'v == 0)
            end
        end
    end
end

module PolynomialSets
using LinearAlgebra
using Polyhedra
using SumOfSquares
import MultivariatePolynomials as MP
import MultivariateBases as MB
import DynamicPolynomials
import SetProg.Sets:
    AbstractSet,
    ConvexPolySet,
    ConvexPolynomialSet,
    Householder,
    PolySet,
    UnknownSet,
    _print_gauge_function,
    dimension,
    gauge1,
    perspective_variable,
    scaling_function,
    space_variables,
    zero_eliminate

function PolySet(degree::Int, p::GramMatrix{T,B,U}) where {T,B,U}
    return PolySet{T,B,U,typeof(p)}(degree, p)
end
function ConvexPolySet(degree::Int, p::GramMatrix{T,B,U}, ::Nothing) where {T,B,U}
    return ConvexPolySet{T,B,U,typeof(p),Nothing}(degree, p, nothing)
end
function ConvexPolynomialSet(
    degree::Int,
    q::GramMatrix{T,B,U},
    z::V,
    x::Vector{V},
) where {T,B,U,V}
    return ConvexPolynomialSet{T,B,U,typeof(q),V}(degree, q, z, x)
end
function Householder(
    set::AbstractSet{T},
    p::DynamicPolynomials.Polynomial{O,M,U},
    h::Vector{Float64},
    z::V,
    x::Vector{V},
) where {T,O,M,U,V}
    return Householder{T,typeof(set),U,typeof(p),V}(set, p, h, z, x)
end
include("Sets/polynomials.jl")
include("Sets/transformations.jl")
end

include("L1_heuristic.jl")
include("constraints.jl")
include("eval.jl")
include("inclusion.jl")
include("separation.jl")
include("map.jl")
include("membership.jl")
include("tangent.jl")
include("utilities.jl")
include("variables.jl")

end
