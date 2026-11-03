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
