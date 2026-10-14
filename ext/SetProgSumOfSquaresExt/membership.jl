
# TODO simplify when https://github.com/JuliaOpt/SumOfSquares.jl/issues/3 is done
function _moment_matrix(f, monos)
    # Takes into account that some monomials are the same so they should
    # have the same value
    first_exponents = MP.exponents(first(monos)^2)
    first_value = f(1, 1)
    d = Dict(first_exponents => first_value)
    function element(i, j)
        mono = monos[i] * monos[j]
        exps = MP.exponents(mono)
        if !haskey(d, exps)
            d[exps] = f(i, j)
        end
        return d[exps]
    end
    return MomentMatrix{typeof(first_value)}(element, monos)
end

function JuMP.add_constraint(
    model::JuMP.Model,
    constraint::SetProg.MembershipConstraint{
        Vector{T},
        <:SetProg.Sets.ConvexPolySet{Float64},
    },
    name::String = "",
) where {T}
    set = constraint.set
    @assert iseven(set.degree)
    vars = SetProg.Sets.space_variables(set)
    monos = MP.monomials(vars, 0:div(set.degree, 2))
    U = promote_type(T, Int, JuMP.VariableRef)
    function variable(i, j)
        mono = monos[i] * monos[j]
        if isone(mono)
            return one(U)
        elseif MP.degree(mono) == 1
            index = findfirst(==(mono), vars)
            return convert(U, constraint.member[index])
        else
            return convert(U, @variable(model))
        end
    end
    ν = _moment_matrix(variable, monos)
    @constraint(model, ν.Q.Q in MOI.PositiveSemidefiniteConeTriangle(length(monos)))
    p = SetProg.Sets.gauge1(set)
    # The cache shares symbolic entries for repeated monomials. Build the
    # moment vector from unique monomials, without numeric equality checks.
    moments = Dict(zip(MultivariateMoments.vectorized_basis(ν), ν.Q.Q))
    basis = MB.SubBasis{MB.Monomial}(MP.monomial_vector(collect(keys(moments))))
    μ = moment_vector([moments[m] for m in MB.keys_as_monomials(basis)], basis)
    scalar_product = expectation(μ, MP.polynomial(p))
    return @constraint(model, scalar_product in MOI.LessThan(1.0), base_name = name)
end
