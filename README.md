# SetProg

| **Build Status** | **Social** |
|:----------------:|:----------:|
| [![Build Status][build-img]][build-url] | [![Gitter][gitter-img]][gitter-url] |
| [![Codecov branch][codecov-img]][codecov-url] | [<img src="https://upload.wikimedia.org/wikipedia/commons/thumb/a/af/Discourse_logo.png/799px-Discourse_logo.png" width="64">][discourse-url] |

JuMP extension for Set Programming : optimization with set variables and inclusion/containment constraints. This package allows the formulation of a mathematical programming involving both classical variables and constraints supported by JuMP and set variables and constraints.

Three options exist to solve Set Programming:
* [Polyhedral Computation](https://github.com/JuliaPolyhedra/Polyhedra.jl).
* Automatically reformulation into a semidefinite program using [Sum-Of-Squares Programming](https://github.com/JuliaOpt/SumOfSquares.jl) and the S-procedure.
* [Dual Dynamic Programming](https://github.com/JuliaStochOpt/StructDualDynProg.jl).

In fact, the option to use can be automatically chosen depending on the variables created and the objective function set:

| Variable/Objective | Volume of set  | Affine of point |
|--------------------|----------------|-----------------|
| Polyhedron         | Polyhedral     | Dual Dynamic    |
| Ellipsoid/PolySet  | Sum-of-Squares | Sum-of-Squares  |

## Variables

The variables can either be
* an Ellipsoid or more generally the 1-sublevel set of a polynomial of degree `2d`;
* a polyhedron (*not yet implemented*);
* a quadratic cone or more generally the 0-sublevel set of a polynomial of degree `2d` (*not yet implemented*).

```julia
@variable m S Ellipsoid()
@variable m S PolySet(d) # 1-sublevel set of a polynomial of degree 2d
@variable m S PolySet(d, convex=true) # Convex 1-sublevel set of a polynomial of degree 2d
@variable m S PolySet(d, symmetric=true) # 1-sublevel set of a polynomial of degree 2d symmetric around the origin
@variable m S PolySet(d, symmetric=true, point=SetProg.CenterPoint([1, 0])) # 1-sublevel set of a polynomial of degree 2d symmetric around the [1, 0]
```

*not yet implemented*:
```julia
@variable m S Polyhedron()
@variable m S QuadCone()  # Quadratic cone
@variable m S PolyCone(d) # 0-sublevel set of a polynomial of degree 2d
```

## Expressions

The following operations are allowed:

| Operation | Description                   |
|-----------|-------------------------------|
| A\*S      | Linear mapping                |

But more operations are planned to be added:

| Operation | Description                   |
|-----------|-------------------------------|
| S + x     | Translation of `S` by `x`     |
| S1 + S2   | Minkowski sum                 |
| S1 ∩ S2   | Intersection of `S1` and `S2` |
| S1 ∪ S2   | Union of `S1` and `S2`        |
| polar(S)  | Polar of S                    |

## Constraints

The following constraints are implemented

| Operation | Description              |
|-----------|--------------------------|
| x ∈ S     | `x` is contained in `S`  |
| S1 ⊆ S2   | `S1` is included in `S2` |
| S1 ⊇ S2   | `S1` is included in `S2` |

But more are planned to be added:

| Operation | Description              |
|-----------|--------------------------|
| S1 == S2  | `S1` is equal to `S2`    |

## Examples

Consider a polytope
```julia
using Polyhedra
P = HalfSpace([1, 1], 1) ∩ HalfSpace([-1, 0], 0) ∩ HalfSpace([0, -1], 0)
```
Pick an SDP solver (see [here](juliaopt.org) for a list)
```julia
using CSDP # Optimizer
factory = with_optimizer(CSDP.Optimizer)
```

To compute the maximal ellipsoid contained in the polytope `P` defined above (i.e. [Löwner-John ellipsoid](https://github.com/rdeits/LoewnerJohnEllipsoids.jl)):
```julia
using JuMP
model = Model(factory)
@variable(model, S, Ellipsoid())
@constraint(model, S ⊆ P)
@objective(model, Max, nth_root(volume(S)))
optimize!(model)
```

To compute the maximal invariant set contained in a polytope (*not yet implemented*):
```julia
using JuMP
model = Model(factory)
@variable(model, S, Polyhedron())
@constraint(model, S ⊆ P)
@constraint(model, A*S ⊆ S) # Invariance constraint
@objective(model, Max, volume(S))
optimize!(model)
```

To compute the maximal invariant ellipsoid contained in the polytope `P` defined above:
```julia
using JuMP
model = Model(factory)
@variable(model, S, Ellipsoid())
@constraint(model, S ⊆ P)
@constraint(model, A*S ⊆ S) # Invariance constraint
@objective(model, Max, nth_root(volume(S)))
optimize!(model)
```

To compute the maximal algebraic-invariant ellipsoid (i.e. `AS ⊆ ES`) contained in the polytope `P` defined above:
```julia
using JuMP
model = Model(factory)
@variable(model, S, Ellipsoid())
@constraint(model, S ⊆ P)
@constraint(model, A*S ⊆ E*S) # Invariance constraint
@objective(model, Max, L1_heuristic(volume(S), zeros(Polyhedra.fulldim(P))))
optimize!(model)
```

[build-img]: https://travis-ci.org/blegat/SetProg.jl.svg?branch=master
[build-url]: https://travis-ci.org/blegat/SetProg.jl
[codecov-img]: http://codecov.io/github/blegat/SetProg.jl/coverage.svg?branch=master
[codecov-url]: http://codecov.io/github/blegat/SetProg.jl?branch=master

[gitter-url]: https://gitter.im/JuliaPolyhedra/Lobby?utm_source=share-link&utm_medium=link&utm_campaign=share-link
[gitter-img]: https://badges.gitter.im/JuliaPolyhedra/Lobby.svg
[discourse-url]: https://discourse.julialang.org/c/domain/opt
