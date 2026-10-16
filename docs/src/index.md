# SetProg

[SetProg](https://github.com/blegat/SetProg.jl) is a JuMP extension for Set Programming :
optimization with set variables and inclusion/containment constraints.
This package allows the formulation of a mathematical programming involving both classical variables and constraints supported by JuMP and set variables and constraints.

## Variables on fixed domains

The `Domain` wrapper restricts a set variable to a fixed polyhedral region.
For an orthant whose dimension is inferred from the model, use
`MOI.Nonnegatives`:

```julia
using SetProg
model = Model()
@variable(model, S, Domain(Polytope(symmetric = true), MOI.Nonnegatives))
orthant = polyhedron(hrep([-1.0 0.0; 0.0 -1.0], zeros(2)))
unit = SetProg.Sets.Intersection(SetProg.Sets.PolarPoint(ones(2)), orthant)
@constraint(model, S ⊆ unit)
```

Here the unit simplex supplies dimension two and bounds the linear gauge's
coefficients below by one. `MOI.Nonnegatives(n)` specifies the dimension
explicitly. A fixed Polyhedra representation or polyhedron is also accepted as
the second argument, with its dimension checked against the wrapped family.
The wrapper also accepts other families, such as
`Domain(Ellipsoid(symmetric = true), MOI.Nonnegatives)`, and constructs a
`Sets.Intersection`. Available constraints depend on the underlying set type.
Domain-restricted variables currently support primal modeling only.

## Polar variable families

Wrap a family in `polar` to construct the polar of its instantiated set:

```julia
@variable(model, B, polar(Domain(Polytope(symmetric = true), MOI.Nonnegatives)))
```

For positive coefficients `c`, the inner set is the positive simplex
`{x ≥ 0 : c'x ≤ 1}`; its polar has coordinate bounds `xᵢ ≤ cᵢ`.
Intersecting that polar with the nonnegative orthant gives a positive box.
The wrapper forwards dimension and polynomial coordinates, and
`polar(polar(family))` returns `family`.
