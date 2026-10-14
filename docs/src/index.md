# SetProg

[SetProg](https://github.com/blegat/SetProg.jl) is a JuMP extension for Set Programming :
optimization with set variables and inclusion/containment constraints.
This package allows the formulation of a mathematical programming involving both classical variables and constraints supported by JuMP and set variables and constraints.

## Homogeneous gauge families

SetProg supports homogeneous set families for stability certificates:
`Ellipsoid(symmetric = true)`, `PolySet(symmetric = true, degree = 4)`,
`PositiveSimplex`, `PositiveBox`, `Parallelotope`, and `ConicPolytope`.
Ellipsoids and polytopes also support a prescribed piecewise partition.
Polynomial sets require `using SumOfSquares` to load the optional extension.

A fitted set is callable as its homogeneous sublevel function, `V(x)`.
Ellipsoids use `x'Qx`; polynomial sets use their defining polynomial;
the other families use degree-one gauges. The homogeneous degree is returned
by `SetProg.homogeneity_degree`.

```@docs
PositiveSimplex
PositiveBox
Parallelotope
ConicPolytope
SetProg.planar_conic_partition
SetProg.with_dimension
SetProg.add_function_variables!
SetProg.add_normalization!
SetProg.add_decrease!
SetProg.solution_value
SetProg.gauge_value
SetProg.homogeneity_degree
SetProg.decrease_slack
SetProg.supports_dual
SetProg.dual_family
SetProg.dual_value
SetProg.unit_set
```
