const _SignedSet = Union{Sets.ShiftedEllipsoid,Sets.PolySet}
SetProg._fixed_set(set::Sets.PolySet) = MP.coefficient_type(MP.polynomial(set.p)) <: Number

function polynomial_multiplier(model, variables, degree, name)
    basis = MP.monomials(variables, 0:div(degree, 2))
    return @variable(model, variable_type = SOSPoly(basis), base_name = name)
end
function _zero_polynomial(set::Sets.ShiftedEllipsoid, variables)
    return dot(variables, set.Q * variables) + 2dot(set.b, variables) + set.β
end
_zero_polynomial(set::Sets.PolySet, variables) = MP.polynomial(set.p) - 1
_scaling_polynomial(λ::Union{Real,JuMP.AbstractJuMPScalar}) = λ
_scaling_polynomial(λ) = MP.polynomial(λ)

function _quadratic_certificate_matrix(p, variables)
    n = length(variables)
    matrix = zeros(JuMP.AffExpr, n + 1, n + 1)
    matrix[1, 1] = MP.coefficient(p, one(MP.monomial(first(variables))))
    for i = 1:n
        matrix[1, i+1] = MP.coefficient(p, MP.monomial(variables[i])) / 2
        matrix[i+1, 1] = matrix[1, i+1]
        matrix[i+1, i+1] = MP.coefficient(p, MP.monomial(variables[i]^2))
        for j = (i+1):n
            matrix[i+1, j+1] =
                MP.coefficient(p, MP.monomial(variables[i] * variables[j])) / 2
            matrix[j+1, i+1] = matrix[i+1, j+1]
        end
    end
    return matrix
end
function _polynomial_separation(
    model,
    first,
    second,
    disjoint;
    S_procedure_scaling = nothing,
    margin = disjoint ? 1 : 0,
)
    variables = if first isa Sets.PolySet
        Sets.space_variables(first)
    elseif second isa Sets.PolySet
        Sets.space_variables(second)
    else
        new_polyvars(Sets.dimension(first))
    end
    λ, fixed_first = SetProg._scaling(model, first, second, S_procedure_scaling, variables)
    λ = _scaling_polynomial(λ)
    a, b = _zero_polynomial(first, variables), _zero_polynomial(second, variables)
    p = if disjoint
        fixed_first ? λ * a + b : a + λ * b
    else
        fixed_first ? λ * a - b : a - λ * b
    end
    degree = max(2, MP.maxdegree(p))
    p -= margin * (1 + sum(variables .^ 2))^cld(degree, 2)
    if degree <= 2
        return SetProg._matrix_certificate(
            model,
            _quadratic_certificate_matrix(p, variables),
        )
    end
    return @constraint(model, p in SOSCone())
end
function JuMP.add_constraint(
    model::JuMP.Model,
    c::InclusionConstraint{<:_SignedSet,<:_SignedSet},
    name::String = "",
)
    return _polynomial_separation(model, c.subset, c.supset, false; c.kws...)
end
function JuMP.add_constraint(
    model::JuMP.Model,
    c::SetProg.DisjointConstraint{<:_SignedSet,<:_SignedSet},
    name::String = "",
)
    return _polynomial_separation(model, c.first, c.second, true; c.kws...)
end
function polynomial_certificate_matrix(cref)
    object = JuMP.constraint_object(cref)
    if object.set isa SumOfSquares.PositiveSemidefinite2x2ConeTriangle
        values = JuMP.value.(object.func)
        return Symmetric([values[1] values[2]; values[2] values[3]])
    end
    return Matrix(SumOfSquares.gram_matrix(cref).Q)
end

SetProg.set_space(
    space::Space,
    ::InclusionConstraint{<:Sets.PolySet,<:SetProg.SetVariableRef},
) = SetProg.set_space(space, PrimalSpace)
