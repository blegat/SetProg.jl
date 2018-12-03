using DataStructures

struct SpaceIndex
    value::Int
end

struct Spaces
    indices::IntDisjointSets # Union-Find datastructure
    dimensions::Dict{Int, Int}
    polyvars::Dict{Int, Vector{SpaceVariable}}
end
Spaces() = Spaces(IntDisjointSets(0), Dict{Int, Int}(),
                  Dict{Int, Vector{SpaceVariable}}())
function new_space(spaces::Spaces)
    return SpaceIndex(push!(spaces.indices))
end
function new_space(spaces::Spaces,
                   dim::Int)
    space_index = new_space(spaces)
    spaces.dimensions[space_index.value] = dim
    return space_index
end
function new_space(spaces::Spaces,
                   polyvars::Vector{SpaceVariable})
    space_index = new_space(spaces)
    spaces.dimensions[space_index.value] = length(polyvars)
    spaces.polyvars[space_index.value] = polyvars
    return space_index
end
function merge_property(d::Dict, root, key1, key2, name)
    value = nothing
    if key1 in keys(d)
        v1 = d[key1]
        if key2 in keys(d) && v1 != d[key2]
            error("Sets lie on the same spaces but have different $name")
        end
        value = v1
    elseif key2 in keys(d)
        value = d[key2]
    end
    if value !== nothing
        d[root] = value
    end
end
function merge_spaces(spaces::Spaces, a::SpaceIndex, b::SpaceIndex)
    root_a = find_root(spaces.indices, a.value)
    root_b = find_root(spaces.indices, b.value)
    root = root_union!(spaces.indices, root_a, root_b)
    # the properties my not be set at a.value when it is set at root_a.value
    merge_property(spaces.dimensions, root, root_a, root_b, "dimension")
    merge_property(spaces.polyvars, root, root_a, root_b,
                   "polynomial variables")
    return SpaceIndex(root)
end
function space_dimension(spaces::Spaces, si::SpaceIndex)
    idx = find_root(spaces.indices, si.value)
    if !haskey(spaces.dimensions, idx)
        error("Missing dimension information, use Ellipsoid(dimension=...) or PolySet(dimension=...)")
    end
    return spaces.dimensions[idx]
end
function space_polyvars(spaces::Spaces, si::SpaceIndex)
    idx = find_root(spaces.indices, si.value)
    if !haskey(spaces.polyvars, idx)
        dim = space_dimension(spaces, si)
        @polyvar x[1:dim]
        spaces.polyvars[idx] = x
    end
    return spaces.polyvars[idx]
end
