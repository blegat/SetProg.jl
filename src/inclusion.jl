### InclusionConstraint ###
struct InclusionConstraint{SubSetType, SupSetType} <: SetConstraint
    subset::SubSetType
    supset::SupSetType
end
function need_variablify(c::InclusionConstraint)
    return need_variablify(c.subset) || need_variablify(c.supset)
end
function variablify(c::InclusionConstraint)
    return InclusionConstraint(variablify(c.subset), variablify(c.supset))
end
function clear_spaces(c::InclusionConstraint)
    clear_spaces(c.subset)
    clear_spaces(c.supset)
end
function Sets.perspective_variable(c::InclusionConstraint)
    return synchronize_perspective(Sets.perspective_variable(c.subset),
                                   Sets.perspective_variable(c.supset))
end
function create_spaces(c::InclusionConstraint, spaces::Spaces)
    sub = create_spaces(c.subset, spaces)
    sup = create_spaces(c.supset, spaces)
    merge_spaces(spaces, sub, sup)
end

JuMP.function_string(print_mode, c::InclusionConstraint) = string(c.subset)
function JuMP.in_set_string(print_mode, c::InclusionConstraint)
    string(print_mode == JuMP.IJuliaMode ? "\\subseteq" : "⊆", " ", c.supset)
end
struct PowerSet{S}
    set::S
end
function JuMP.build_constraint(_error::Function, subset,
                               supset_powerset::PowerSet)
    InclusionConstraint(subset, supset_powerset.set)
end

### InclusionConstraint for sets ###

## Set in Set ##
function set_space(space::Space, ::InclusionConstraint{<:LinearImage,
                                                       <:LinearImage})
    return set_space(space, DualSpace)
end
# We can always transform  an ellipsoid to primal or dual space so we can handle
# any space
function set_space(space::Space,
                   ::InclusionConstraint{<:VariableRef,
                                         <:Sets.AbstractEllipsoid{T}}) where T<:Number
    return space
end

# S-procedure: Q ⊆ P <=> xQx ≤ 1 => xPx ≤ 1 <=> xPx ≤ xQx <=> Q - P is PSD
function JuMP.add_constraint(model::JuMP.Model,
                             constraint::InclusionConstraint{<:Sets.EllipsoidAtOrigin,
                                                             <:Sets.EllipsoidAtOrigin},
                             name::String = "")
    Q = constraint.subset.Q
    P = constraint.supset.Q
    @constraint(model, Symmetric(Q - P) in PSDCone())
end

# S-procedure: Q ⊆ P <=> q(x) ≤ 1 => p(x) ≤ 1 <=> p(x) ≤ q(x) <= q - p is SOS
function JuMP.add_constraint(model::JuMP.Model,
                             constraint::InclusionConstraint{<:Sets.ConvexPolynomialSublevelSetAtOrigin,
                                                             <:Sets.ConvexPolynomialSublevelSetAtOrigin},
                             name::String = "")
    q = constraint.subset.p
    p = constraint.supset.p
    @constraint(model, q - p in SOSCone())
end

# S-procedure: Q ⊆ P <=> Q* ⊇ P* <= p - q is SOS
function JuMP.add_constraint(model::JuMP.Model,
                             constraint::InclusionConstraint{<:Union{Sets.DualQuadCone,
                                                                     Sets.DualPolynomialSet},
                                                             <:Union{Sets.DualQuadCone,
                                                                     Sets.DualPolynomialSet}},
                             name::String = "")
    q = constraint.subset.p
    p = constraint.supset.p
    @constraint(model, p - q in SOSCone()) # TODO λ
end


# S ⊆ T <=> polar(T) ⊆ polar(S)
function JuMP.add_constraint(model::JuMP.Model,
                             constraint::InclusionConstraint{<:Union{Sets.PolarEllipsoidAtOrigin,
                                                                     Sets.PolarConvexPolynomialSublevelSetAtOrigin},
                                                             <:Union{Sets.PolarEllipsoidAtOrigin,
                                                                     Sets.PolarConvexPolynomialSublevelSetAtOrigin}},
                             name::String = "")
    S = constraint.subset
    T = constraint.supset
    @constraint(model, Sets.polar(T) ⊆ Sets.polar(S))
end

# See [LTJ18]
function JuMP.add_constraint(model::JuMP.Model,
                             constraint::InclusionConstraint{<:LinearImage{S},
                                                             <:LinearImage{T}},
                             name::String = "") where {S <: Sets.AbstractSet,
                                                       T <: Sets.AbstractSet}
    @constraint(model, apply_map(model, constraint.subset) ⊆ apply_map(model, constraint.supset))
end

## Set in Polyhedron ##
function set_space(space::Space, ::InclusionConstraint{<:VariableRef,
                                                       <:Polyhedra.Rep})
    return set_space(space, DualSpace)
end

# Ellipsoid #
function JuMP.add_constraint(model::JuMP.Model,
                             constraint::InclusionConstraint{<:Sets.AbstractSet{JuMP.VariableRef},
                                                             <:Polyhedra.Rep},
                             name::String = "")
    ◯ = constraint.subset
    □ = constraint.supset
    for hp in hyperplanes(□)
        @constraint(model, ◯ ⊆ hp)
    end
    for hs in halfspaces(□)
        @constraint(model, ◯ ⊆ hs)
    end
end
function JuMP.add_constraint(model::JuMP.Model,
                             constraint::InclusionConstraint{<:Union{Sets.PolarEllipsoidAtOrigin{JuMP.VariableRef},
                                                                     Sets.PolarConvexPolynomialSublevelSetAtOrigin{JuMP.VariableRef}},
                                                             <:Polyhedra.HyperPlane},
                             name::String = "")
    @assert iszero(constraint.supset.β) # Otherwise it is not symmetric around the origin
    @constraint(model, Line(constraint.supset.a) in Sets.polar(constraint.subset))
end
function JuMP.add_constraint(model::JuMP.Model,
                             constraint::InclusionConstraint{<:Union{Sets.DualQuadCone,
                                                                     Sets.DualConvexPolynomialCone},
                                                             <:Polyhedra.HyperPlane},
                             name::String = "")
    val = sublevel_eval(model, constraint.subset, constraint.supset.a,
                        constraint.supset.β)
    @constraint(model, val in MOI.EqualTo(0.0))
end
function JuMP.add_constraint(model::JuMP.Model,
                             constraint::InclusionConstraint{<:Union{Sets.PolarEllipsoidAtOrigin{JuMP.VariableRef},
                                                                     Sets.PolarConvexPolynomialSublevelSetAtOrigin{JuMP.VariableRef}},
                                                             <:Polyhedra.HalfSpace},
                             name::String = "")
    @constraint(model, ScaledPoint(constraint.supset.a, constraint.supset.β) in Sets.polar(constraint.subset))
end
function JuMP.add_constraint(model::JuMP.Model,
                             constraint::InclusionConstraint{<:Union{Sets.DualQuadCone,
                                                                     Sets.DualConvexPolynomialCone},
                                                             <:Polyhedra.HalfSpace},
                             name::String = "")
    val = sublevel_eval(model, constraint.subset, constraint.supset.a,
                        constraint.supset.β)
    @constraint(model, val in MOI.LessThan(0.0))
end

## Polyhedron in Set ##
function set_space(space::Space, ::InclusionConstraint{<:Polyhedra.Rep,
                                                       <:VariableRef})
    return set_space(space, PrimalSpace)
end

function JuMP.add_constraint(model::JuMP.Model,
                             constraint::InclusionConstraint{<:Polyhedra.Rep,
                                                             <:Sets.AbstractSet{JuMP.VariableRef}},
                             name::String = "")
    □ = constraint.subset
    ◯ = constraint.supset
    for line in lines(□)
        @constraint(model, line in ◯)
    end
    for ray in rays(□)
        @constraint(model, ray in ◯)
    end
    for point in points(□)
        @constraint(model, point in ◯)
    end
end
