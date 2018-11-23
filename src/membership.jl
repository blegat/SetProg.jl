### MembershipConstraint ###
struct MembershipConstraint{P, S} <: SetConstraint
    member::P
    set::S
end
JuMP.function_string(print_mode, c::MembershipConstraint) = string(c.member)
JuMP.in_set_string(print_mode, c::MembershipConstraint) = string(JuMP.math_symbol(print_mode, :in), c.set)
function JuMP.build_constraint(_error::Function, member,
                               set::Sets.AbstractSet)
    MembershipConstraint(member, set)
end

function JuMP.add_constraint(model::JuMP.Model,
                             constraint::MembershipConstraint{<:AbstractVector,
                                                              <:Union{Sets.EllipsoidAtOrigin,
                                                                      Sets.PolynomialSublevelSetAtOrigin}},
                             name::String = "")
    @constraint(model, sublevel_eval(constraint.set, constraint.member) in MOI.LessThan(1.0))
end
function JuMP.add_constraint(model::JuMP.Model,
                             constraint::MembershipConstraint{<:Polyhedra.Line,
                                                              <:Union{Sets.EllipsoidAtOrigin,
                                                                      Sets.PolynomialSublevelSetAtOrigin}},
                             name::String = "")
    # We must have (λl)^T Q (λl) ≤ 1 for all λ hence we must have l^T Q l ≤ 0
    # As Q is positive definite, it means l^T Q l = 0
    l = Polyhedra.coord(constraint.member)
    @constraint(model, sublevel_eval(constraint.set, l) in MOI.EqualTo(0.0))
end
function JuMP.add_constraint(model::JuMP.Model,
                             constraint::MembershipConstraint{<:Polyhedra.Ray,
                                                              <:Union{Sets.EllipsoidAtOrigin,
                                                                      Sets.PolynomialSublevelSetAtOrigin}},
                             name::String = "")
    # We must have (λl)^T Q (λl) ≤ 1 for all λ > 0 hence we must have l^T Q l ≤ 0
    # As Q is positive definite, it means l^T Q l = 0
    r = Polyhedra.coord(constraint.member)
    @constraint(model, sublevel_eval(constraint.set, r) in MOI.EqualTo(0.0))
end

