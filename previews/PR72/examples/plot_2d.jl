using Test     #src
# # 2D plotting
#
#md # [![Binder](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/generated/plot_2d.ipynb)
#md # [![nbviewer](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/generated/plot_2d.ipynb)
#
# In this notebook we show how to plot the sublevel set of $2x^4 + 2y^4$ and
# its polar. We first build the corresponding Gram matrix: the polynomial is
# $\langle z, Q z \rangle$ where $z = (x^2, y^2)$ and $Q = 2I$.

using LinearAlgebra
using DynamicPolynomials
@polyvar x y
using SetProg
gram = SetProg.SumOfSquares.GramMatrix(Matrix{Float64}(2I, 2, 2), [x^2, y^2])

# Now we create the set from this Gram matrix representation of the polynomial.
# We claim that it is convex (as we want to plot the polar) by creating a
# `ConvexPolySet` even if we give `nothing` in place of the convexity
# certificate.

set = SetProg.Sets.ConvexPolySet(4, gram, nothing)

# We can now plot the set together with its polar:

using Plots
plot(ratio=:equal)
plot!(SetProg.Sets.polar(set), npoints=128, alpha=0.6)
plot!(set, npoints=128, alpha=0.95)
