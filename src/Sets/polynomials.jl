using MultivariatePolynomials

abstract type AbstractPolySet{T} <: AbstractSet{T} end

struct PolynomialSublevelSet{T, P<:AbstractPolynomial{T}}
    p::P
end
