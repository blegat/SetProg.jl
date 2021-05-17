struct PolarPoint{T} <: AbstractSet{T}
    a::Vector{T}
end
dimension(h::PolarPoint) = length(h.a)

function scaling_function(h::PolarPoint)
    @assert dimension(h) == 2
    return (x, y) -> begin
        return h.a[1] * x + h.a[2] * y
    end
end

function _print_gauge_function(h::PolarPoint; digits=6)
    DynamicPolynomials.@polyvar x[1:2]
    print(" ")
    a = h.a
    if digits !== nothing
        a = round.(a, digits=6)
    end
    println(x' * a)
end
