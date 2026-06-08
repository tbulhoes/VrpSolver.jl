# Custom Resources

Before implementing a custom resource, read the following paper — it describes
the Meta-Solver interface, the required C++ functions, and the theoretical
background behind the design choices:

> [**Bucket Graph Meta-Solver for the Resource Constrained Shortest Path Problem**](https://inria.hal.science/hal-05486295)  
> Ruslan Sadykov, Aurélien Froger, Eduardo Uchoa, Artur Pessoa, Teobaldo Bulhões, Daniel de Araujo  
> *Preprint, 2026. HAL: hal-05486295*

A *custom resource* extends the pricing solver with new feasibility or cost logic
implemented directly in C++, complementing the standard resource consumption mechanism. Instead of the standard
lower-bound / upper-bound / consumption triple, a custom resource carries
arbitrary user-defined parameters for each arc, each vertex, and globally.

Custom resources are an **advanced feature**. Only one custom resource can
be defined per [`VrpGraph`](@ref).

The [CCVRP demo](https://github.com/artalvpes/VRPSolverDemos/tree/main/CCVRP)
is a complete working example — it minimises the weighted sum of customer
arrival times using a custom resource.

## Workflow

### 1. Implement the internal C++ code

The RCSP repository already contains the files with placeholder implementations:

- `Tools/rcsp/include_dev/rcsp_custom_res_impl.hpp`
- `Tools/rcsp/src/rcsp_custom_res_impl.cpp`

Fill them in following the instructions in the header and the companion paper.
The [CCVRP demo](https://github.com/artalvpes/VRPSolverDemos/tree/main/CCVRP/src/meta_solver)
provides a complete reference implementation.

Once the files are filled in, **recompile the BaPCod library**.

### 2. Define matching Julia structs

The Julia side is an interface to the internal C++ code. You must define three
structs — for arc, vertex, and global parameters — whose fields exactly mirror
the three `CustomRes*Parameters` structs you declared in C++, in the same order.
Every field must be a C primitive type (`Cdouble`, `Cint`, `Cfloat`, …). An
empty struct is fine when a parameter level is unused.

In the CCVRP demo, the C++ structs are:

```cpp
struct CustomResArcParameters {
    double t;   // travel time of the arc
    double w;   // demand weight of the arc
};

struct CustomResVertexParameters {};   // no vertex-level parameters

struct CustomResConstParameters {
    double Wmax;   // vehicle capacity
    double Tmax;   // upper bound on total time
};
```

And the matching Julia structs are:

```julia
struct CumulativeResArcParams
    t::Cdouble   # travel time
    w::Cdouble   # demand weight
end

struct CumulativeResVertexParams end   # no vertex-level parameters

struct CumulativeResConstParams
    Wmax::Cdouble   # vehicle capacity
    Tmax::Cdouble   # upper bound on total time
end
```

### 3. Register the types

Call [`@register_custom_res_param_types`](@ref) at the **top level of your
module** (outside any function), before `optimize!` runs:

```julia
@register_custom_res_param_types(CumulativeResArcParams, CumulativeResVertexParams, CumulativeResConstParams)
```

### 4. Add the custom resource to the graph

Use [`add_resource!`](@ref) with `custom=true`. Pass `cost_var` when
`isCostResource()` returns `true` in your C++ implementation:

```julia
cum_res_id = add_resource!(G, custom=true, cost_var=z)
```

### 5. Set parameter values

```julia
# Global parameters (once per graph)
set_const_custom_res_params!(G, cum_res_id, CumulativeResConstParams(Wmax, Tmax))

# Per-arc parameters
for arc_id in arcs
    set_arc_custom_res_params!(G, arc_id, cum_res_id, CumulativeResArcParams(t, w))
end

# Per-vertex parameters (omit if VertexParams is empty)
# set_vertex_custom_res_params!(G, v, cum_res_id, ...)
```

## Constraints and limitations

- **At most one custom resource per graph.** However, the C++ implementation
  is free to embed several logical resources into a single custom resource —
  for example, by storing multiple state fields and combining their extension,
  domination, and concatenation logic inside the required functions.
- **Cannot be `main` or `binary`.**
- **BaPCod must be recompiled** every time the C++ files change.

## API reference

```@docs
@register_custom_res_param_types
```

```@docs
set_arc_custom_res_params!
```

```@docs
set_vertex_custom_res_params!
```

```@docs
set_const_custom_res_params!
```
