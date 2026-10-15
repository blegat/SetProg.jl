# SetProg

[SetProg](https://github.com/blegat/SetProg.jl) is a JuMP extension for Set Programming :
optimization with set variables and inclusion/containment constraints.
This package allows the formulation of a mathematical programming involving both classical variables and constraints supported by JuMP and set variables and constraints.

## Set families

SetProg supports the set families
`Ellipsoid(symmetric = true)`, `PolySet(symmetric = true, degree = 4)`,
and `Polytope`. Polynomial sets require `using SumOfSquares` to load the
optional extension.

Polyhedral configurations share linear pieces and polar operations:

```julia
simplex = Polytope(piecewise = :orthant)
box = Polytope(piecewise = :orthant, polar = true)
weighted_l1 = Polytope(symmetric = true, piecewise = :orthants)
parallelotope = Polytope(symmetric = true, piecewise = :orthants, polar = true)
conic = Polytope(symmetric = true, piecewise = cones)
```

A matrix in `cones` has generating rays as columns. Symbolic partitions infer
the dimension from the model. `polar` describes the requested set family;
it does not impose a modeling space. Constraints select that space as usual.
On `:orthant`, the polar is restricted to the nonnegative orthant.

Set evaluation is available on fitted sets through `V(x)`.
Certificate construction and normalization are provided by PathCompleteCertificates.

```@docs
Polytope
SetProg.gauge_value
```
