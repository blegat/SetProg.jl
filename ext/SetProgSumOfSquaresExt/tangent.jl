function linear_algebraic_surface(set::Sets.PolarOf{<:Sets.ConvexPolySet}, A, E)
    v = MP.variables(set.set.p)
    DynamicPolynomials.@polyvar z[1:size(E, 1)]
    Ez = E' * z
    return JuMP.build_constraint(
        error,
        -dot(A' * z, [d(v => Ez) for d in MP.differentiate(set.set.p, v)]),
        SOSCone(),
    )
end
