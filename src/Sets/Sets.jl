module Sets
struct Ellipsoid{T}
    Q::Matrix{T}
    c::Vector{T}
end
using LinearAlgebra
using RecipesBase
@recipe function f(ell::Ellipsoid)
    @assert LinearAlgebra.checksquare(ell.Q) == 2
    αs = range(0, stop=2π, length=1024)
    ps = [[cos(α), sin(α)] for α in αs]
    r = [sqrt(dot(p, ell.Q * p)) for p in ps]
    seriestype --> :shape
    legend --> false
    ell.c[1] .+ cos.(αs) ./ r, ell.c[2] .+ sin.(αs) ./ r
end
end
