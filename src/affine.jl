struct AffineTerm{T, F}
    coefficient::T
    func::F
end

struct AffineExpression{T, F <: AbstractScalarFunction} <: AbstractScalarFunction
    terms::Vector{AffineTerm{T, F}}
    constant::T
end

function Base.:(+)(f::F, g::F) where F <: AbstractScalarFunction
    return AffineExpression([AffineTerm(1.0, f), AffineTerm(1.0, g)], 0.0)
end
function Base.:(+)(a::AffineExpression{T, F}, f::F) where {T, F <: AbstractScalarFunction}
    return AffineExpression([a.terms; AffineTerm(1.0, f)], a.constant)
end
function Base.:(+)(f::F, a::AffineExpression{T, F}) where {T, F <: AbstractScalarFunction}
    return a + f
end

function objective_sense(model::JuMP.Model, f::AffineExpression)
    sense = objective_sense(model, f.terms[1])
    for i in 2:length(f.terms)
        @assert sense == objective_sense(model, f.terms[i])
    end
    return sense
end

function objective_function(model::JuMP.Model, aff::AffineExpression)
    obj = convert(JuMP.AffExpr, aff.constant)
    for t in aff.terms
        obj = JuMP.destructive_add!(obj, t.coefficient,
                                    objective_function(t.func))
    end
    return obj
end
