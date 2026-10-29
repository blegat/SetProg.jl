# Set-variable operands are resolved by the same loading machinery as maps.
Sets._lp_scalar_type(::SetVariableRef) = JuMP.AffExpr
Sets._lp_dimension(set::SetVariableRef) = set.set.dimension

need_variablify(set::Sets.AbstractLpSum) =
    need_variablify(set.first) || need_variablify(set.second)
function variablify(set::LpSum)
    return LpSum(variablify(set.first), variablify(set.second), set.p)
end
function variablify(set::PolarLpSum)
    return PolarLpSum(variablify(set.first), variablify(set.second), set.p)
end
function JuMP.value(set::LpSum)
    return LpSum(JuMP.value(set.first), JuMP.value(set.second), set.p)
end
function JuMP.value(set::PolarLpSum)
    return PolarLpSum(JuMP.value(set.first), JuMP.value(set.second), set.p)
end
function clear_spaces(set::Sets.AbstractLpSum)
    clear_spaces(set.first)
    clear_spaces(set.second)
end
function create_spaces(set::Sets.AbstractLpSum, spaces::Spaces)
    first = create_spaces(set.first, spaces)
    second = create_spaces(set.second, spaces)
    return merge_spaces(spaces, first, second)
end
needs_polyvars(set::Sets.AbstractLpSum) =
    needs_polyvars(set.first) || needs_polyvars(set.second)

const _PolarLpExpression = Sets.PolarOf{<:Sets.AbstractLpSum}
need_variablify(set::_PolarLpExpression) = need_variablify(set.set)
variablify(set::_PolarLpExpression) = polar(variablify(set.set))
clear_spaces(set::_PolarLpExpression) = clear_spaces(set.set)
create_spaces(set::_PolarLpExpression, spaces::Spaces) = create_spaces(set.set, spaces)
needs_polyvars(set::SetVariableRef) = needs_polyvars(set.set)
needs_polyvars(::AbstractVariable) = true
needs_polyvars(set::Ellipsoid) = !set.symmetric
needs_polyvars(::Polytope) = false
needs_polyvars(set::Domain) = needs_polyvars(set.set)

# A squared-gauge sum can be written as a Schur-complement certificate when
# each term is a fixed ellipsoid or the polar of an affine matrix ellipsoid.
# Keeping inverse matrices in the diagonal blocks avoids products of variables.
function _quadratic_lp_terms(set::PolarLpSum)
    return _quadratic_lp_terms(Val(set.p == 2 ? 2 : set.p), set)
end
function _quadratic_lp_terms(::Val{2}, set::PolarLpSum)
    return (_quadratic_lp_terms(set.first)..., _quadratic_lp_terms(set.second)...)
end
function _quadratic_lp_terms(set::Sets.LinearPreImage)
    return _quadratic_lp_terms(set.set, set.A)
end
function _quadratic_lp_terms(set::Sets.Ellipsoid{<:Number}, A)
    return ((A, inv(Matrix(set.Q))),)
end
function _quadratic_lp_terms(set::Sets.PolarOf{<:Sets.Ellipsoid}, A)
    return ((A, polar(set).Q),)
end
function _quadratic_lp_terms(
    set::Union{Sets.Ellipsoid{<:Number},Sets.PolarOf{<:Sets.Ellipsoid}},
)
    return _quadratic_lp_terms(
        set,
        Matrix{Bool}(I, Sets.dimension(set), Sets.dimension(set)),
    )
end

function JuMP.build_constraint(
    _error::Function,
    subset::Sets.Ellipsoid,
    superset::PowerSet{<:PolarLpSum};
    psd_margin::Real = 0,
)
    psd_margin >= 0 || throw(ArgumentError("psd_margin must be nonnegative"))
    terms = _quadratic_lp_terms(superset.set)
    n = Sets.dimension(subset)
    for (A, D) in terms
        size(A, 2) == n && size(A, 1) == LinearAlgebra.checksquare(D) ||
            throw(DimensionMismatch())
    end
    T = promote_type(
        JuMP.AffExpr,
        eltype(subset.Q),
        (eltype(A) for (A, _) in terms)...,
        (eltype(D) for (_, D) in terms)...,
        typeof(psd_margin),
    )
    dimension = n + sum(size(A, 1) for (A, _) in terms)
    matrix = zeros(T, dimension, dimension)
    matrix[1:n, 1:n] .= subset.Q
    offset = n
    for (A, D) in terms
        indices = (offset+1):(offset+size(A, 1))
        matrix[indices, 1:n] .= A
        matrix[1:n, indices] .= A'
        matrix[indices, indices] .= D
        offset += size(A, 1)
    end
    return JuMP.build_constraint(_error, matrix - psd_margin * I, JuMP.PSDCone())
end

function JuMP.value(set::Sets.LinearPreImage)
    return Sets.LinearPreImage(JuMP.value(set.set), JuMP.value.(set.A))
end
