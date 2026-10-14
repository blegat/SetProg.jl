
function polar_perspective_ellipsoid(ell, point::HintPoint, z, x::AbstractVector)
    y = [z; x]
    H = Sets._householder(point.h)
    p = y' * H * Sets._perspective_cat(ell) * H * y
    return Sets.perspective_dual(Sets.Householder(ell, p, point.h, z, x))
end

function perspective_dual_polyset(set, point::HintPoint, z, x::AbstractVector)
    y = [z; x]
    H = Sets._householder(point.h)
    p = Sets.perspective_gauge0(set)(y => H * y)
    return Sets.perspective_dual(Sets.Householder(set, p, point.h, z, x))
end
# For `CenterPoint`, `q` should be non-perspective, need refactoring
function perspective_dual_polyset(degree, q, point::InteriorPoint, z, x)
    set = Sets.ConvexPolynomialSet(degree, q, z, x)
    perspective_dual_polyset(set, point, z, x)
end
function polyset(;
    point::Union{Nothing,HintPoint} = nothing,
    symmetric::Bool = false,
    degree::Union{Int,Nothing} = nothing,
    dimension::Union{Int,Nothing} = nothing,
    convex::Bool = false,
    variables::Union{AbstractVector,Nothing} = nothing,
    superset::Union{Sets.PolySet,Nothing} = nothing,
    basis::Type = MultivariateBases.Monomial,
)
    if degree === nothing
        error("Degree of PolySet not specified, use PolySet(degree=..., ...)")
    end
    if isodd(degree)
        throw(ArgumentError("Degree of PolySet not even"))
    end
    if dimension === nothing
        if point !== nothing
            dimension = length(point.h)
        end
    end
    if superset !== nothing
        if dimension === nothing
            dimension = Sets.dimension(superset)
        elseif dimension != Sets.dimension(superset)
            throw(DimensionMismatch())
        end
        if variables === nothing
            variables = Sets.space_variables(superset)
        elseif variables != Sets.space_variables(superset)
            error("Space variables set does not correspond to superset space variables.")
        end
    end
    return PolySet(point, symmetric, degree, dimension, convex, variables, superset, basis)
end

function constrain_convex(model, p, vars)
    hessian = MP.differentiate(p, vars, 2)
    # We do not just do `@constraint(model, p in SOSConvex())` as we would
    # like to have access to the PSD matrix of variables for the det volume heuristic
    y = [
        MP.similar_variable(eltype(hessian), gensym()) for
        i = 1:LinearAlgebra.checksquare(hessian)
    ]
    q = dot(y, hessian * y)
    multipartite = SumOfSquares.Certificate.NewtonDegreeBounds(tuple(y))
    X = SumOfSquares.Certificate.monomials_half_newton_polytope(
        MP.monomials(q),
        multipartite,
    )
    # If `X` is empty, we will need the following bridge
    JuMP.add_bridge(model, SumOfSquares.Bridges.Constraint.EmptyBridge)
    # If `length(X)` is 2, we will need the following bridge
    JuMP.add_bridge(model, SumOfSquares.Bridges.Constraint.PositiveSemidefinite2x2Bridge)
    set = SumOfSquares.matrix_cone(MOI.PositiveSemidefiniteConeTriangle, length(X))
    Q = @variable(model, [1:MOI.dimension(set)])
    @constraint(model, Q in set)
    s = SumOfSquares.build_gram_matrix(
        Q,
        MultivariateBases.SubBasis{MultivariateBases.Monomial}(X),
        MOI.PositiveSemidefiniteConeTriangle,
        Float64,
    )
    @constraint(model, q == s)
    return MultivariateMoments.value_matrix(s)
end

function variable_set(
    model::JuMP.AbstractModel,
    set::PolySet,
    space::Space,
    space_dimension,
    space_polyvars,
)
    n = space_dimension
    d = data(model)
    # General all monomials of degree `degree`, we don't want monomials of
    # lower degree as the polynomial is homogeneous
    @assert iseven(set.degree)
    if set.symmetric
        monos = MP.monomials(space_polyvars, div(set.degree, 2))
    else
        monos = MP.monomials(lift_space_variables(d, space_polyvars), div(set.degree, 2))
    end
    basis = MultivariateBases.explicit_basis_covering(
        MultivariateBases.FullBasis{set.basis}(MP.variables(monos)),
        MultivariateBases.SubBasis{MultivariateBases.Monomial}(monos),
    )
    # TODO If `set.convex` and `set.symmetric`, no need for the poly to be SOS,
    # see Lemma 6.33 of [BPT12].
    p = @variable(model, variable_type=SOSPoly(basis))
    if set.convex
        set.superset === nothing || error("superset not supported for convex PolySet")
        if set.symmetric
            convexity_proof = constrain_convex(model, p, space_polyvars)
            if space == PrimalSpace
                return Sets.ConvexPolySet(set.degree, p, convexity_proof)
            else
                @assert space == DualSpace
                return Polyhedra.polar(Sets.ConvexPolySet(set.degree, p, convexity_proof))
            end
        else
            constrain_convex(model, MP.subs(p, d.perspective_polyvar => 1), space_polyvars)
            if space == PrimalSpace
                error("Non-symmetric PolySet in PrimalSpace not implemented yet")
            else
                @assert space == DualSpace
                if set.point === nothing
                    throw(
                        ArgumentError(
                            "Specify a point for nonsymmetric polyset, e.g. `PolySet(point=InteriorPoint([1.0, 0.0]))",
                        ),
                    )
                end
                return perspective_dual_polyset(
                    set.degree,
                    p,
                    set.point,
                    d.perspective_polyvar,
                    space_polyvars,
                )
            end
        end
    else
        if set.symmetric
            if space == PrimalSpace
                if set.superset !== nothing
                    p = SumOfSquares.gram_operate(+, set.superset.p, p)
                end
                return Sets.PolySet(set.degree, p)
            else
                error("Non-convex PolySet not supported in $space")
            end
        else
            error("Non-convex nonsymmetric PolySet not implemented yet.")
        end
    end
end
function _value(convexity_proof::MultivariateMoments.SymMatrix)
    return MultivariateMoments.SymMatrix(JuMP.value.(convexity_proof.Q), convexity_proof.n)
end
