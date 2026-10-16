# SetProg

[SetProg](https://github.com/blegat/SetProg.jl) is a JuMP extension for Set Programming :
optimization with set variables and inclusion/containment constraints.
This package allows the formulation of a mathematical programming involving both classical variables and constraints supported by JuMP and set variables and constraints.

## Variables on fixed domains

The `domain` keyword restricts a `Polytope` variable to a fixed polyhedral
region. For an orthant whose dimension is inferred from the model, use
`MOI.Nonnegatives`:

```julia
using SetProg
model = Model()
@variable(model, S, Polytope(domain = MOI.Nonnegatives))
orthant = polyhedron(hrep([-1.0 0.0; 0.0 -1.0], zeros(2)))
unit = SetProg.Sets.Intersection(SetProg.Sets.PolarPoint(ones(2)), orthant)
@constraint(model, S ⊆ unit)
```

Here the unit simplex supplies dimension two and bounds the linear gauge's
coefficients below by one. `MOI.Nonnegatives(n)` specifies the dimension
explicitly. A fixed Polyhedra representation or polyhedron is also accepted as
`domain`, with its dimension checked against the other variable data.
Domain-restricted variables currently support primal modeling only.
