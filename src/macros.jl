function JuMP.parse_one_operator_constraint(_error::Function, vectorized::Bool,
                                            ::Val{:⊂}, lhs, rhs)
    _error("Unrecognized symbol ⊂ you mean ⊆ ?")
end
function JuMP.parse_one_operator_constraint(_error::Function, vectorized::Bool,
                                            ::Val{:⊆}, lhs, rhs)
    parse_code = :()
    if vectorized
        build_call = :(JuMP.build_constraint.($_error, $(esc(lhs)), $(esc(:(SetProg.PowerSet.($rhs))))))
    else
        build_call = :(JuMP.build_constraint($_error, $(esc(lhs)), $(esc(:(SetProg.PowerSet($rhs))))))
    end
    return parse_code, build_call
end
function JuMP.parse_one_operator_constraint(_error::Function, vectorized::Bool,
                                            ::Val{:⊃}, lhs, rhs)
    _error("Unrecognized symbol ⊃, did you mean ⊇ ?")
end
function JuMP.parse_one_operator_constraint(_error::Function, vectorized::Bool,
                                            ::Val{:⊇}, lhs, rhs)
    parse_one_operator_constraint(_error, vectorized, Val(:⊆), rhs, lhs)
end
