export boundary_point, tangent_cone

struct BoundaryPoint{S} <: SymbolicVariable
    set::S
    symbol::Symbol
end
struct TangentCone{S} <: SymbolicVariable
    set::S
    symbol::Symbol
end
boundary_point(S, s::Symbol) = BoundaryPoint(S, s)
tangent_cone(S, s::Symbol) = TangentCone(S, s)
function tangent_cone(S, b::BoundaryPoint)
    @assert S === b.set
    tangent_cone(S, b.symbol)
end
variablify(t::BoundaryPoint) = BoundaryPoint(variablify(t.set), t.symbol)
need_variablify(t::TangentCone) = need_variablify(t.set)
variablify(t::TangentCone) = TangentCone(variablify(t.set), t.symbol)
create_spaces(b::BoundaryPoint, spaces::Spaces) = create_spaces(b.set, spaces)
create_spaces(t::TangentCone, spaces::Spaces) = create_spaces(t.set, spaces)
function create_spaces(c::MembershipConstraint, spaces::Spaces)
    sub = create_spaces(c.member, spaces)
    sup = create_spaces(c.set, spaces)
    return merge_spaces(spaces, sub, sup)
end
clear_spaces(b::BoundaryPoint) = clear_spaces(b.set)
clear_spaces(t::TangentCone) = clear_spaces(t.set)
function clear_spaces(c::MembershipConstraint)
    clear_spaces(c.member)
    clear_spaces(c.set)
end
Sets.perspective_variable(b::BoundaryPoint) = Sets.perspective_variable(b.set)
Sets.perspective_variable(t::TangentCone) = Sets.perspective_variable(t.set)
function Sets.perspective_variable(c::MembershipConstraint)
    return synchronize_perspective(
        Sets.perspective_variable(c.member),
        Sets.perspective_variable(c.set)
    )
end
function set_space(
    space::Space,
    ::MembershipConstraint{
        <:Sets.LinearImage{<:BoundaryPoint},
        <:Sets.LinearImage{<:TangentCone}}
)
    return set_space(space, DualSpace)
end
using LinearAlgebra
function linear_algebraic_surface(set::Sets.PolarOf{<:Sets.Ellipsoid}, A, E)
    Q = set.set.Q
    return psd_constraint(Symmetric(-A * Q * E' - E * Q * A'))
end
function linear_algebraic_surface(set::Sets.PolarOf{<:Sets.ConvexPolySet}, A, E)
    v = variables(set.set.p)
    DynamicPolynomials.@polyvar z[1:size(E, 1)]
    Ez = E' * z
    return JuMP.build_constraint(
        error,
        -dot(A' * z, [d(v => Ez) for d in differentiate(set.set.p, v)]),
        SOSCone(),
    )
end
function JuMP.build_constraint(
    _error::Function,
    member::Sets.LinearImage{BoundaryPoint{S}},
    set::Sets.LinearImage{TangentCone{S}}
   ) where S <: Sets.PolarOf{<:Union{Sets.Ellipsoid, Sets.ConvexPolySet}}
    _set = member.set.set
    @assert set.set.set === _set
    return linear_algebraic_surface(_set, member.A, set.A)
end
function add_linear_algebraic_surface_domain(
    model, set::Sets.Ellipsoid, domain, A, E)
    Q = set.Q
    P = Symmetric(-A * Q * E' - E * Q * A')
    return _add_constraint_or_not(model, psd_in_domain(model, P, E' \ domain))
end
function add_linear_algebraic_surface(model, set::Sets.PolarOf{<:Sets.Piecewise}, A, E)
    for (i, si) in enumerate(set.set.sets)
        add_linear_algebraic_surface_domain(
            model, si, set.set.pieces[i], A, E
        )
    end
end
function JuMP.add_constraint(
    model::JuMP.Model,
    constraint::MembershipConstraint{
        <:Sets.LinearImage{<:BoundaryPoint},
        <:Sets.LinearImage{<:TangentCone}
    }
)
    _set = constraint.member.set.set
    @assert constraint.set.set.set === _set
    return add_linear_algebraic_surface(model, _set, constraint.member.A, constraint.set.A)
end
