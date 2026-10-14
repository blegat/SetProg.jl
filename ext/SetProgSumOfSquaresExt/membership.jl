
# TODO simplify when https://github.com/JuliaOpt/SumOfSquares.jl/issues/3 is done
function _moment_matrix(f, monos)
    # Takes into account that some monomials are the same so they should
    # have the same value
    d = Dict{typeof(exponents(first(monos))),typeof(f(1, 1))}()
    function element(i, j)
        mono = monos[i] * monos[j]
        exps = exponents(mono)
        if !haskey(d, exps)
            d[exps] = f(i, j)
        end
        return d[exps]
    end
    return MomentMatrix(element, monos)
end

function JuMP.add_constraint(
    model::JuMP.Model,
    constraint::SetProg.MembershipConstraint{Vector{T},SetProg.Sets.ConvexPolySet{Float64}},
    name::String = "",
) where {T}
    set = constraint.set
    @assert iseven(set.degree)
    monos = monomials(SetProg.Sets.space_variables(set), 0:div(set.degree, 2))
    U = promote_type(T, Int, JuMP.VariableRef)
    function variable(i, j)
        if i == length(monos) && j == length(monos)
            return one(U)
        elseif i == length(monos) || j == length(monos)
            return convert(U, constraint.member[length(monos)-min(i, j)])
        else
            return convert(U, JuMP.VariableRef(model))
        end
    end
    ν = _moment_matrix(variable, monos)
    @constraint(model, ν.Q.Q in MOI.PositiveSemidefiniteConeTriangle(length(monos)))
    p = SetProg.Sets.gauge1(set)
    scalar_product = dot(measure(ν), polynomial(p))
    @constraint(model, scalar_product in MOI.LessThan(1.0))
end
