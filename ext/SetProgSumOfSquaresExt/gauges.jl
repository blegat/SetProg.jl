function SetProg.gauge_value(set::Union{Sets.PolySet,Sets.ConvexPolySet}, x::AbstractVector)
    return MP.polynomial(set.p)(Sets.space_variables(set) => collect(x))
end
