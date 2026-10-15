"""
    Intersection(set, domain)

The intersection of a SetProg set and a fixed polyhedral domain.
"""
struct Intersection{T,S<:AbstractSet{T},D<:Polyhedra.Polyhedron} <: AbstractSet{T}
    set::S
    domain::D
end
space_variables(set::Intersection) = space_variables(set.set)
perspective_variable(set::Intersection) = perspective_variable(set.set)
dimension(set::Intersection) = dimension(set.set)
