using LinearAlgebra
using DynamicPolynomials
@polyvar x y
using SetProg
gram = SetProg.SumOfSquares.GramMatrix(Matrix{Float64}(2I, 2, 2), [x^2, y^2])

set = SetProg.Sets.ConvexPolySet(4, gram, nothing)

using Plots
plot(ratio=:equal)
plot!(SetProg.Sets.polar(set), npoints=128, alpha=0.6)
plot!(set, npoints=128, alpha=0.95)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

