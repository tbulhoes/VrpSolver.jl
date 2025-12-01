# Resources

This page discusses rules and limitations of the current version of the VRPSolver w.r.t. the different types of resources. A resource can be created in VRPSolver through the function [`add_resource!`](@ref).

## Monotone x Non-monotone

Resources without any negative consumption are called *monotone*, otherwise they are *non-monotone*. This is automatically detected by the VRPSolver.

## Main x Secondary

A resource can be *main* or *secondary*, where main resources should be monotone. Moreover, there should not exist a cycle in the [`VrpGraph`](@ref) with zero consumption of all main resources. Therefore, unless the `VrpGraph` is acyclic, it is mandatory the existence of at least one main resource. By default, resources created with [`add_resource!`](@ref) are secondary. Main resources are directly linked to the definition of the bucket graphs in the labeling algorithm, which improves significantly the pricing algorithm performance when compared with secondary resources. The first main resource is a special one, since the consumption of the first main resource determines the threshold for the bi-directional labelling. 

## Disposable x Non-disposable

A resource constrained path $p=(v_{\mathrm{source}}^k=v_0,a_1,v_1,\ldots,a_{n-1},v_{n-1},a_n,v_{n}=v_{\mathrm{sink}}^k)$ over a graph $G^k$ should have $n \geq 1$ arcs, $v_j \neq v_{\mathrm{source}}^k$ and $v_j \neq v_{\mathrm{sink}}^k$, $1\leq j\leq n-1$, and is feasible if: 
* for every $r \in R^k$ that is *disposable*, the accumulated resource consumption $S_{j,r}$ at visit $j$, $0\leq j\leq n$, where $S_{0,r}=0$ and $S_{j,r}=\max\{l_{a_j,r},S_{j-1,r}+q_{a_j,r}\}$, does not exceed $u_{a_j,r}$;
* for every $r \in R^k$ that is *non-disposable*, the accumulated resource consumption $S_{j,r}$ at visit $j$, $0\leq j\leq n$, where $S_{0,r}=0$ and $S_{j,r}= S_{j-1,r}+q_{a_j,r}$, lies in the interval $[l_{a_j,r},u_{a_j,r}]$.

Main resources should be disposable.

## Binary x Non-binary

A resource is binary if its accumulated consumption can only be `0` or `1`. Sets of secondary binary resources of the same type (disposable or not disposable) are implemented in a special way, being represented as bitsets in the labels. This greatly decreases the time spent for dominance checks between labels.

Currently, all binary resources must be secondary in VRPSolver and all consumptions must be `-1`, `0` or `1`, and all intervals `[0,0]`, `[0,1]` or `[1,1]`. Moreover, non-disposable binary resources cannot have the interval `[0,1]` set for the node $v_{\mathrm{sink}}^k$ of a graph $G^k$.

## Other remarks

If a `VrpGraph` is acyclic, in principle it is possible to not have any resource. In this case, an artificial main resource should be defined with zero consumption for all arcs. Also in this case bidirectional labelling should not be used (see [Parameters](@ref)).

Currently, VRPSolver has the following limits for resources in a VrpGraph: 
* at most 2 main resources
* at most 20 non-binary resources (including main ones)
* at most 512 binary resources.

