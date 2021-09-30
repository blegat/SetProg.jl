using Test     #src
# # Continuous-time Controlled Invariant Set with State-Dependent Switching
#
#md # [![Binder](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Continuous-time Controlled Invariant Set with State-Dependent Switching.ipynb)
#md # [![nbviewer](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Continuous-time Controlled Invariant Set with State-Dependent Switching.ipynb)
#
# ## Introduction
#

function buck(R, L, Ro, Co)
    A1 = [-R/L  -1/L
           1/Co -1/(Ro * Co)]
    B1 = [1/L, 0.0]
    B2 = zeros(2)
    return (A1, A1), (B1, B2)
end

function boost(R, L, Ro, Co)
    A1 = [-R/L 0
           0   -1/(Ro * Co)]
    A2 = [-R/L  -1/L
           1/Co -1/(Ro * Co)]
    B1 = [1/L, 0.0]
    return (A1, A2), (B1, B1)
end

function buck_boost(R, L, Ro, Co)
    A1 = [-R/L 0
           0   -1/(Ro * Co)]
    A2 = [-R/L  -1/L
           1/Co -1/(Ro * Co)]
    B1 = [1/L, 0.0]
    B2 = zeros(2)
    return (A1, A2), (B1, B2)
end

using SetProg
function maximal_invariant(family, γ = nothing; dirs=dirs)
    model = Model()
    @variable(model, S, family)
    x = boundary_point(S, :x)
    @constraint(model, C * x in E * tangent_cone(S, x))
    return 0, 0
end

sol_ell, γ_ell = maximal_invariant(Ellipsoid(symmetric=true))
