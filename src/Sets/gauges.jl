"""
    Intersection(set, domain)

The intersection of a SetProg set and a fixed polyhedral domain. This also
represents blocking polars, whose domain remains the nonnegative orthant.
"""
struct Intersection{T,S<:AbstractSet{T},D<:Polyhedra.Polyhedron} <: AbstractSet{T}
    set::S
    domain::D
end
space_variables(set::Intersection) = space_variables(set.set)
perspective_variable(set::Intersection) = perspective_variable(set.set)
dimension(set::Intersection) = dimension(set.set)
