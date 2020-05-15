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

struct PowerOfLinearForm
    coefficients::Vector{Int}
    power::Int
end
struct Decomposition
    coefficients::Vector{Int}
    forms::Vector{PowerOfLinearForm}
    deno::Int
end

using Combinatorics
using Polyhedra
using LinearAlgebra
function _volume_simplex(s)
    return abs(det([hcat(collect(points(s))...); ones(npoints(s))']))
end
function volume_simplex(s)
    return _volume_simplex(s) / factorial(fulldim(s))
end
# See (6) of [BBDKV11]
function _linear_simplex(l::PowerOfLinearForm, s)
    ls = [l.coefficients'si for si in points(s)]
    sum(prod(i -> ls[i]^k[i], eachindex(ls)) for k in multiexponents(npoints(s), l.power))
end
function _integrate(l::PowerOfLinearForm, s)
    frac = (factorial(l.power) / factorial(l.power + fulldim(s))) # /!\ TODO: overflow
    return frac * _volume_simplex(s)  * _linear_simplex(l, s)
end
function integrate(l::PowerOfLinearForm, s, cache::Dict{PowerOfLinearForm, Float64})
    if !haskey(cache, l)
        cache[l] = _integrate(l, s)
    end
    return cache[l]
end

function integrate_simplex!(ints::Vector{Float64}, decs::Vector{Decomposition}, simplex,
                            cache = Dict{PowerOfLinearForm, Float64}())
    empty!(cache)
    for (i, dec) in enumerate(decs)
        ints[i] += sum(dec.coefficients[j] * integrate(dec.forms[j], simplex, cache)
                       for j in eachindex(dec.forms)) / dec.deno
    end
end
function integrate_decompositions(decs::Vector{Decomposition}, polytope::Polyhedron,
                                  cache = Dict{PowerOfLinearForm, Float64}())
    ints = zeros(length(decs))
    for simplex in simplices(polytope)
        integrate_simplex!(ints, decs, simplex, cache)
    end
    return ints
end
function integrate_monomials(monos::AbstractVector{<:AbstractMonomial}, polytope::Polyhedron)
    return integrate_decompositions(decompose.(monos), polytope)
end
function integrate(p::AbstractPolynomial, polytope::Polyhedron)
    return coefficients(p)'integrate_monomials(monomials(p), polytope)
end

"""
As studied in [DPW96].

[DPW96] C. Durieu, B. T. Polyak and E. Walter.
*Trace versus determinant in ellipsoidal outer-bounding, with application to
state estimation*.
IFAC Proceedings Volumes, **1996**.
"""
function l1_integral(ell::Sets.Ellipsoid, vertex)
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
function l1_integral(set::Union{Sets.PolySet,
                                Sets.ConvexPolySet},
                     vertex)
    return rectangle_integrate(set.p, vertex)
end

# See (13) of [BBDKV11]
#[BBDKV11] Baldoni, V., Berline, N., De Loera, J., Köppe, M., & Vergne, M. (2011).
# How to integrate a polynomial over a simplex. Mathematics of Computation, 80(273), 297-325.
function decompose(exps::Vector{Int})
    f(i) = 0:i
    d = sum(exps)
    deno = factorial(d)
    coefficients = Int[]
    forms = PowerOfLinearForm[]
    for p in Iterators.product(f.(exps)...)
        dp = sum(p)
        iszero(dp) && continue
        coef = prod(i -> binomial(exps[i], p[i]), eachindex(p))
        if isodd(d - dp)
            coef = -coef
        end
        push!(coefficients, coef)
        push!(forms, PowerOfLinearForm(collect(p), d))
    end
    return Decomposition(coefficients, forms, deno)
end
decompose(mono::AbstractMonomial) = decompose(exponents(mono))

using QHull
function simplices(v::VRepresentation)
    if hasrays(v)
        error("Not a polytope")
    end
    ch = QHull.chull(MixedMatVRep(v).V)
    return [vrep(ch.points[simplex, :]) for simplex in ch.simplices]
end
simplices(polytope::Polyhedron) = simplices(vrep(polytope))
function all_exponents(set::Sets.Piecewise{<:Any, <:Sets.PolySet})
    # TODO check that all polysets have same degree
    s = set.sets[1]
    return exponents.(monomials(Sets.space_variables(s), s.degree))
end
function ell_exponents(i, j, n)
    exps = zeros(Int, n)
    exps[i] += 1
    exps[j] += 1
    return exps
end
function all_exponents(set::Sets.Piecewise{<:Any, <:Sets.Ellipsoid})
    return [ell_exponents(i, j, Sets.dimension(set)) for j in 1:Sets.dimension(set) for i in 1:j]
end
function add_integrals!(total, integral::Function, set::Sets.PolySet)
    p = polynomial(set.p)
    for t in terms(p)
        total = MA.add_mul!(total, integral(exponents(monomial(t))), coefficient(t))
    end
    return total
end
function add_integrals!(total, integral::Function, set::Sets.Ellipsoid)
    for j in 1:Sets.dimension(set)
        for i in 1:j
            total = MA.add_mul!(total, integral(ell_exponents(i, j, Sets.dimension(set))), set.Q[i, j])
        end
    end
    return total
end
function l1_integral(set::Sets.Piecewise{T, <:Union{Sets.Ellipsoid{T}, Sets.PolySet{T}}},
                     vertex) where T
    decs = Decomposition[]
    val = Dict(exps => length(push!(decs, decompose(exps)))
               for exps in all_exponents(set))
    cache = Dict{PowerOfLinearForm, Float64}()
    ints = [zeros(length(decs)) for i in 1:length(set.sets)]
    for simplex in simplices(set.polytope)
        piece = findfirst(h -> simplex ⊆ Polyhedra.hyperplane(h), collect(halfspaces(set.polytope)))
        integrate_simplex!(ints[piece], decs, convexhull(simplex, Polyhedra.origin(Polyhedra.pointtype(set.polytope), Sets.dimension(set))), cache)
    end
    U = MA.promote_operation(*, Float64, T)
    total = zero(MA.promote_operation(+, U, U))
    for (integrals, set) in zip(ints, set.sets)
        total = add_integrals!(total, exps -> integrals[val[exps]], set)
    end
    return total
end

# The polar of the rectangle with vertices (-v, v) is not a rectangle but the
# smaller rectangle contained in it has vertices (-1 ./ v, 1 ./ v).
# the set is not necessarily inside this rectangle, we just know that it is
# outside the polar of the rectangle with vertices (-v, v). However, since it
# is homogeneous, only the ratios between the dimensions is important
function l1_integral(set::Sets.PolarOf{<:Union{Sets.Ellipsoid,
                                               Sets.Piecewise, # vertex ignored by Piecewise anyway ^^
                                               Sets.ConvexPolySet}},
                     vertex)
    return l1_integral(Sets.polar(set), 1 ./ vertex)
end
function l1_integral(set::Sets.HouseDualOf{<:Sets.AbstractEllipsoid},
                     vertex)
    return l1_integral(Sets.polar(Sets.Ellipsoid(set.set.set.Q)),
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
function invert_objective_sense(::Union{Sets.Ellipsoid,
                                        Sets.PolySet,
                                        Sets.ConvexPolySet,
                                        Sets.Piecewise})
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
