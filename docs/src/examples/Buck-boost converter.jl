using Test     #src
# # Continuous-time Controlled Invariant Set with State-Dependent Switching
#
#md # [![Binder](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/Continuous-time Controlled Invariant Set with State-Dependent Switching.ipynb)
#md # [![nbviewer](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/Continuous-time Controlled Invariant Set with State-Dependent Switching.ipynb)
#
# ## Introduction

# Taken from [SwitchOnSafety](https://github.com/blegat/SwitchOnSafety.jl).
# We threat zero rows separately to avoid having approximate identity in `E` as it can lead to numerical difficuly for the SDP solver down the road.

using LinearAlgebra
function algebraiclift(A, B)
    n = LinearAlgebra.checksquare(A)
    z = findall(i -> all(j -> iszero(B[i, j]), 1:size(B, 2)), 1:n)
    nz = setdiff(1:n, z)
    Ez = Matrix(one(eltype(A)) * I, n, n)[z, :]
    null = transpose(nullspace(transpose(B[nz, :])))
    Enz = zeros(eltype(null), size(null, 1), n)
    Enz[:, nz] = null
    return [A[z, :]; Enz * A], [Ez; Enz]
end

function unwind(AB::AbstractArray{Tuple{S,T}}) where {S,T}
    A = similar(AB, S)
    B = similar(AB, T)
    for i in eachindex(AB)
        A[i] = AB[i][1]
        B[i] = AB[i][2]
    end
    return A, B
end

using SetProg
function system(A, B, τ)
    n = length(A)
    R = JuMP.Containers.@container(
        [i = 1:n, j = 1:n; i != j],
        exp(A[j] * τ),
    )
    S = JuMP.Containers.@container(
        [i = 1:n, j = 1:n; i != j],
        A[j] \ (exp(A[j] * τ) - I) * B[j],
    )
    return unwind(algebraiclift.(A, B))..., unwind(algebraiclift.(R, S))...
end

# We use `ms`, `mH` and `mF`

R  = 2      # [Ω]
L  = 500e-3 # [mH]
Co = 470e-3 # [mF]
Ro = 50     # [Ω]
τ  = 1      # [ms]

function buck(R, L, Ro, Co)
    A1 = [-R/L  -1/L
           1/Co -1/(Ro * Co)]
    B1 = [1/L, 0.0]
    B2 = zeros(2)
    return [A1, A1], [B1, B2]
end

function boost(R, L, Ro, Co)
    A1 = [-R/L 0
           0   -1/(Ro * Co)]
    A2 = [-R/L  -1/L
           1/Co -1/(Ro * Co)]
    B1 = B2 = [1/L, 0.0]
    return [A1, A1], [B1, B2]
end

function buck_boost(R, L, Ro, Co)
    A1 = [-R/L 0
           0   -1/(Ro * Co)]
    A2 = [-R/L  -1/L
           1/Co -1/(Ro * Co)]
    B1 = [1/L, 0.0]
    B2 = zeros(2)
    return [A1, A1], [B1, B2]
end

import GLPK
lp_solver = optimizer_with_attributes(GLPK.Optimizer, MOI.Silent() => true, "presolve" => GLPK.GLP_ON)
import CSDP
sdp_solver = optimizer_with_attributes(CSDP.Optimizer) #, MOI.Silent() => true)
using MosekTools
sdp_solver = optimizer_with_attributes(Mosek.Optimizer, MOI.Silent() => true)
using Polyhedra
interval = HalfSpace([1.0], 1.0) ∩ HalfSpace([-1.0], 1.0)
lib = Polyhedra.DefaultLibrary{Float64}(lp_solver)
□ = polyhedron(interval * interval, lib)
dirs = [[-1 + √3, -1 + √3], [-1, 1]]
all_dirs = [dirs; (-).(dirs)]

using SetProg
function maximal_invariant(family, C, E, T, U, γ = nothing; dirs=dirs)
    model = Model(sdp_solver)
    n = length(C)
    @assert length(E) == n
    @variable(model, S[1:n], family)
    @constraint(model, [i in 1:n], S[i] ⊆ □)
    for i in 1:n
        x = boundary_point(S[i], :x)
        @constraint(model, C[i] * x in E[i] * tangent_cone(S[i], x))
    end
    @constraint(model, [i=1:n, j=1:n; i != j], U[i, j] * S[i] ⊆ T[i, j] * S[j])
    if γ === nothing
        @variable(model, γ)
    end
    for point in dirs
        @constraint(model, [point in dirs, i in 1:n], γ * point in S[i])
    end
    @objective(model, Max, γ)
    JuMP.optimize!(model)
    println(model)
    @show solve_time(model)
    @show JuMP.termination_status(model)
    @show JuMP.objective_value(model)
    if JuMP.termination_status(model) == MOI.OPTIMAL
        return JuMP.value(S), JuMP.objective_value(model)
    else
        return
    end
end

maximal_invariant(Ellipsoid(symmetric=true), system(buck(R, L, Ro, Co)..., τ)...)
maximal_invariant(Ellipsoid(symmetric=true), system(boost(R, L, Ro, Co)..., τ)...)
maximal_invariant(Ellipsoid(symmetric=true), system(buck_boost(R, L, Ro, Co)..., τ)...)
#sol_ell, γ_ell = maximal_invariant(Ellipsoid(symmetric=true), system(buck(R, L, Ro, Co)..., τ)...)
