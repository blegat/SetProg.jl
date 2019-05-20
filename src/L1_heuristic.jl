struct L1Heuristic{V <: SetVariableRef} <: AbstractScalarFunction
    variable::V
    rectangle_vertex::Vector{Float64}
end
Base.copy(l::L1Heuristic) = l
L1_heuristic(volume::Volume, v::Vector{Float64}) = L1Heuristic(volume.variable, v)
Base.show(io::IO, l::L1Heuristic) = print(io, "L1-heuristic(", l.variable, ")")

set_space(space::Space, ::L1Heuristic, ::JuMP.Model) = return space

function power_integrate(exponent, bound)
    exp = exponent + 1
    @assert isodd(exp)
    return (2/exp) * bound^exp
end
function rectangle_integrate(m::MultivariatePolynomials.AbstractMonomialLike,
                             vertex)
    exp = exponents(m)
    if any(isodd, exp)
        return zero(power_integrate(2, vertex[1]))
    else
        @assert length(exp) == length(vertex)
        return prod(power_integrate.(exp, vertex))
    end
end
function rectangle_integrate(t::MultivariatePolynomials.AbstractTermLike,
                             vertex)
    return coefficient(t) * rectangle_integrate(monomial(t), vertex)
end
function rectangle_integrate(p::MultivariatePolynomials.AbstractPolynomialLike,
                            vertex)
    return sum(rectangle_integrate(t, vertex) for t in terms(p))
end

"""
As studied in [DPW96].

[DPW96] C. Durieu, B. T. Polyak and E. Walter.
*Trace versus determinant in ellipsoidal outer-bounding, with application to
state estimation*.
IFAC Proceedings Volumes, **1996**.
"""
function l1_integral(ell::Sets.EllipsoidAtOrigin, vertex)
    n = Sets.dimension(ell)
    @assert n == length(vertex)
    # The integral of off-diagonal entries xy are zero between the rectangle
    # is symmetric
    # the integral of x^2 is (x^3 - (-x^3)) / 3 = 2/3 * x^3
    c = 2/3
	return sum(i -> ell.Q[i, i] * c * vertex[i]^2, 1:n)
end

"""
As developed in [HL12].

[HL12] D. Henrion and C. Louembet.
*Convex inner approximations of nonconvex semialgebraic sets applied to
fixed-order controller design*.
International Journal of Control, **2012**.
"""
function l1_integral(set::Union{Sets.PolynomialSublevelSetAtOrigin,
                                Sets.ConvexPolynomialSublevelSetAtOrigin},
                     vertex)
    return rectangle_integrate(set.p, vertex)
end


# The polar of the rectangle with vertices (-v, v) is not a rectangle but the
# smaller rectangle contained in it has vertices (-1 ./ v, 1 ./ v).
# the set is not necessarily inside this rectangle, we just know that it is
# outside the polar of the rectangle with vertices (-v, v). However, since it
# is homogeneous, only the ratios between the dimensions is important
function l1_integral(set::Sets.PolarOf{<:Union{Sets.EllipsoidAtOrigin,
                                               Sets.ConvexPolynomialSublevelSetAtOrigin}},
                     vertex)
    return l1_integral(Sets.polar(set), 1 ./ vertex)
end
function l1_integral(set::Sets.HouseDualOf{<:Sets.AbstractEllipsoid},
                     vertex)
    return l1_integral(Sets.polar(Sets.EllipsoidAtOrigin(set.set.set.Q)),
                       vertex)
end
function l1_integral(set::Sets.HouseDualOf{<:Sets.ConvexPolynomialSet},
                     vertex)
    return rectangle_integrate(subs(set.set.set.q, set.set.set.z => 0),
                               1 ./ vertex)
end

function invert_objective_sense(::Union{Sets.Polar,
                                        Sets.PerspectiveDual})
    return false
end
function invert_objective_sense(::Union{Sets.EllipsoidAtOrigin,
                                        Sets.PolynomialSublevelSetAtOrigin,
                                        Sets.ConvexPolynomialSublevelSetAtOrigin})
    return true
end
function objective_sense(model::JuMP.Model, l::L1Heuristic)
    sense = data(model).objective_sense
    if invert_objective_sense(l.variable.variable)
        if sense == MOI.MAX_SENSE
            return MOI.MIN_SENSE
        elseif sense == MOI.MIN_SENSE
            return MOI.MAX_SENSE
        else
            error("Unsupported objective sense $sense for `L1_heuristic`.")
        end
    else
        return sense
    end
end
function objective_function(::JuMP.Model, l::L1Heuristic)
    return l1_integral(l.variable.variable, l.rectangle_vertex)
end
