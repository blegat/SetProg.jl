
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
    degree::Union{Sets.PolynomialDegrees,Nothing} = nothing,
    dimension::Union{Int,Nothing} = nothing,
    convex::Bool = false,
    psd::Bool = true,
    variables::Union{AbstractVector,Nothing} = nothing,
    superset::Union{Sets.PolySet,Nothing} = nothing,
    basis::Type = MultivariateBases.Monomial,
)
    if degree === nothing
        error("Degree of PolySet not specified, use PolySet(degree=..., ...)")
    end
    isempty(degree) && throw(ArgumentError("the degree selection must not be empty"))
    minimum(degree) >= 0 || throw(ArgumentError("polynomial degrees must be nonnegative"))
    if isodd(maximum(degree))
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
    return PolySet(
        point,
        symmetric,
        degree,
        dimension,
        convex,
        psd,
        variables,
        superset,
        basis,
    )
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
    explicit_degrees = !(set.degree isa Int)
    explicit_degrees &&
        space != PrimalSpace &&
        throw(
            ArgumentError("explicit polynomial degree selections require primal modeling"),
        )
    d = data(model)
    if explicit_degrees
        degrees = div(minimum(set.degree), 2):div(maximum(set.degree), 2)
        monos = MP.monomials(space_polyvars, degrees)
    elseif set.symmetric
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
    p = if set.psd || set.superset !== nothing
        @variable(model, variable_type=SOSPoly(basis))
    else
        Q = @variable(model, [1:length(basis), 1:length(basis)], Symmetric)
        SumOfSquares.GramMatrix(Q, basis)
    end
    if !explicit_degrees && set.basis !== MultivariateBases.Monomial
        # A nonhomogeneous basis, such as Chebyshev, introduces lower-degree
        # terms even when covering only homogeneous monomials.
        for term in MP.terms(MP.polynomial(p))
            if MP.degree(term) != set.degree
                @constraint(model, MP.coefficient(term) == 0)
            end
        end
    end
    if explicit_degrees
        set.point === nothing || throw(
            ArgumentError("point hints are not supported with explicit degree ranges"),
        )
        if set.superset !== nothing
            set.convex && error("superset not supported for convex PolySet")
            p = SumOfSquares.gram_operate(+, set.superset.p, p)
        end
        # A Gram product may introduce degrees absent from the selection.
        # Symmetric sets also exclude odd-degree terms.
        for term in MP.terms(MP.polynomial(p))
            degree = MP.degree(term)
            if degree ∉ set.degree || (set.symmetric && isodd(degree))
                @constraint(model, MP.coefficient(term) == 0)
            end
        end
        if set.convex
            proof = constrain_convex(model, p, space_polyvars)
            return Sets.ConvexPolySet(set.degree, p, proof)
        end
        return Sets.PolySet(set.degree, p)
    end
    if set.convex
        set.superset === nothing || error("superset not supported for convex PolySet")
        if set.symmetric
            # Exclude lower-degree expressions already constrained to zero,
            # so the Hessian certificate does not acquire forced zero rows.
            convex_polynomial =
                set.basis === MultivariateBases.Monomial ? p :
                MP.filter_terms(t -> MP.degree(t) == set.degree, p)
            convexity_proof = constrain_convex(model, convex_polynomial, space_polyvars)
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
