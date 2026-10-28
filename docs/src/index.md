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
`PolySet(degree = 0:4, psd = false, variables = x)` represents `p(x) ≤ 1`
with monomials of degrees zero through four and a free Gram matrix.
An integer, such as `degree = 4`, selects homogeneous polynomials.
Degree ranges and collections work with both values of `psd`; their maximum
must be even. `symmetric = true` excludes odd-degree terms. These
affine families use primal modeling; homogeneous families keep the existing
space selection. Defaults retain their positivity constraints.

## Separation certificates

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
fixed multiplier, defaulting to one. The `Disjoint` constraint defaults to a
margin of one, certifying strict disjointness. A zero-margin certificate alone
can allow the sets to share boundary points.

### Polyhedral separation

A polytope family with a fixed partition also supports ellipsoid containment
and separation from an ellipsoid's exterior:

```julia
@variable(model, S, Polytope(symmetric = true, piecewise = partition))
@constraint(model, initial ⊆ S, margin = ε)
@constraint(model, [S, unsafe] in Disjoint(), margin = ε)
```

Here `initial` is a bounded ellipsoid and `unsafe` is a `ShiftedEllipsoid`
whose quadratic part is negative definite: it describes the exterior of a
bounded ellipsoid. The partition must be bounded and contain the origin in its interior. Both
modeling spaces are supported, and other constraints can select between them:
in primal gauge coordinates, initial containment uses second-order cones and
unsafe separation is linear on the fixed cone rays; in dual vertex
coordinates, initial containment is linear and each vertex is constrained by
an SOC. No SOS dependency is needed.

Margins are representation-dependent. In primal coordinates, the gauge is
at most `1 - ε` on the initial set and exceeds the safe ellipsoid's gauge by
`ε * norm(direction)` on each cone ray. In dual coordinates, the margin pads
the initial support and reduces the safe ellipsoid's Cholesky radius. A
positive margin certifies strict separation in either representation.
`certificate_matrix` assembles all conic parts for numerical checks; `dual`
returns one dual vector per part when a constraint generates several cones.

## Lₚ sums and polar Lₚ sums

`LpSum(C, S, p)` is the Firey sum, defined by its support function:
`h(LpSum(C, S, p), x)^p = h(C, x)^p + h(S, x)^p`.
`PolarLpSum(C, S, p)` is its polar counterpart:
`polar(LpSum(polar(C), polar(S), p))`. Its gauge satisfies
`g(PolarLpSum(C, S, p), x)^p = g(C, x)^p + g(S, x)^p`.
The exponent must be finite and at least one. The operands must represent
convex sets containing the origin with nonnegative gauges.

The exponent belongs to the operation, independently of the set family.
For example, exponent two adds squared polyhedral gauges or the quadratic
gauges of ellipsoids; exponent four adds the defining polynomials of
homogeneous quartic PolySets. Exponent one adds ordinary gauges, while
`LpSum(C, S, 1)` is the Minkowski sum.

```julia
using SetProg, LinearAlgebra
C = SetProg.Sets.Ellipsoid(Symmetric([1.0 0.0; 0.0 1.0]))
S = SetProg.Sets.Ellipsoid(Symmetric([2.0 0.0; 0.0 2.0]))
combined = PolarLpSum(C, S, 2)
SetProg.sublevel_eval(combined, [1.0, 2.0]) # 15.0
support_function([1.0, 2.0], LpSum(C, S, 2)) # √7.5
```

`sublevel_eval` for `PolarLpSum` evaluates the sum of the powered gauges;
the represented set is where this sum is at most one. `support_function`
for `LpSum` evaluates its ordinary, degree-one support function. Evaluating
polynomial gauges requires loading SumOfSquares and homogeneous operands.

These expressions support dimension inference, polynomial coordinates,
set-variable loading, and `JuMP.value`. Constraint certificates and controller
synthesis for these operations are not implemented yet.

See Section 5.5 of Gardner, Hug and Weil,
[*Operations between sets in geometry*](https://doi.org/10.4171/JEMS/422).
Their polar-addition convention uses parameter `-p`; these constructors take
positive `p` for both operations.
