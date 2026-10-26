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
The constraints select the modeling space. `S ⊆ unit` uses the restricted
primal gauge. A normalization written as `polar(unit) ⊆ S` selects dual
modeling, returning the polar of the restricted gauge. For a positive simplex
`unit`, this is the set of points with coordinate bounds `xᵢ ≤ cᵢ`, whose
intersection with the nonnegative orthant is a positive box.

## Sublevel sets without a positivity constraint

Use `psd = false` to omit the PSD constraint on an ellipsoid's coefficient
matrix, or the SOS constraint on a polynomial's Gram matrix. For example,
`Ellipsoid(symmetric = false, psd = false, dimension = 2)` represents
`x'Qx + 2b'x + β ≤ 0`, with free symmetric `Q`, vector `b`, and scalar `β`.
`PolySet(degree = 4, psd = false, variables = x)` represents `p(x) ≤ 1`
with all monomials of degree at most four and a free Gram matrix. These
affine families use primal modeling; homogeneous families keep the existing
space selection. Defaults retain their positivity constraints.

Quadratic and polynomial sublevel sets support S-procedure certificates:

```julia
using SetProg, LinearAlgebra
model = Model()
@variable(model, S, Ellipsoid(psd = false, dimension = 2))
initial = SetProg.Sets.ShiftedEllipsoid(
    LinearAlgebra.Symmetric([1.0 0.0; 0.0 1.0]), zeros(2), -16.0,
)
unsafe = SetProg.Sets.ShiftedEllipsoid(
    LinearAlgebra.Symmetric([-1.0 0.0; 0.0 -1.0]), zeros(2), 36.0,
)
@variable(model, ε >= 0)
@constraint(model, initial ⊆ S, margin = ε)
@constraint(model, [S, unsafe] in Disjoint(), margin = ε)
```

With zero-level polynomials `a` and `b`, inclusion certifies `λa - b ≥ 0`
and disjointness certifies `λa + b ≥ 0` when the first set is fixed. When
the second set is fixed, the multiplier scales `b` instead. A positive
`margin` subtracts `margin * (1 + ‖x‖²)^d` from this certificate, where `2d`
is its degree. Quadratic certificates are PSD matrix constraints; higher
degrees use SOS constraints and require `using SumOfSquares`.

By default a fixed set is scaled by a nonnegative scalar. Pass
`S_procedure_scaling = SProcedureMultiplier(degree = 2)` for an SOS multiplier,
or a nonnegative number for a fixed multiplier. Two variable sets require a
fixed multiplier, defaulting to one. The `Disjoint` constraint defaults to a margin of one, certifying strict
disjointness. A zero-margin certificate
alone can allow the sets to share boundary points.
