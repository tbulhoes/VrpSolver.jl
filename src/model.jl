"""
    VrpModel()

Create an empty VRPSolver model.

It is the main object of the VRPSolver. It is responsible to keep the RCSP graphs (of type VrpGraph), formulation, packing sets, rounded capacity cuts and other definitions (like branching).

"""
function VrpModel()
    VrpModel(
        Model(),
        VrpGraph[],
        Array{Tuple{VrpGraph,Int},1}[],
        VrpSolver.NoSet,
        VrpSolver.NoSet,
        false,
        Dict{String,Int}(),
        Any[],
        Any[],
        true,
        nothing,
        Dict{String,Any}(),
        CapacityCutInfo[],
        CapacityCutInfo[],
        Array{Array{Tuple{VrpGraph,VrpArc},1},2}(undef, 0, 0),
        Vector{Tuple{Integer,Integer,Bool}}(),
    )
end

function _using_resource_consumption_branching(model::VrpModel)
    return haskey(model.branching_priorities, "_res_cons_branching_")
end

function _using_ryan_and_foster_branching(model::VrpModel)
    return haskey(model.branching_priorities, "_ryanfoster_branching_")
end

function _using_strong_kpath_cuts(model::VrpModel)
    return !isempty(model.strongkpath_cuts_info)
end

function _has_overlapping_sets(model::VrpModel)
    for graph in model.graphs
        for vertex in graph.vertices
            if length(vertex.elem_sets) > 1 || length(vertex.packing_sets) > 1
                return true
            end
        end
        for arc in graph.arcs
            if length(arc.elem_sets) > 1 || length(arc.packing_sets) > 1
                return true
            end
        end
    end
    return false
end

function _has_custom_resource(model::VrpModel)
    for graph in model.graphs
        count(r -> r.is_custom, graph.resources) > 0 && return true
    end
    return false
end

"""
    enable_rank1_cuts!(model::VrpModel)

Enable use of *limited memory rank-1* cuts during the execution.

These cuts depends of the collection of packing sets (whether in arcs or vertices).
Therefore, it will not be applied if you do not define packing sets.
These cuts are potentially very strong, but each added cut can make the pricing subproblems harder.

"""
function enable_rank1_cuts!(model::VrpModel)
    model.use_rank1_cuts = true
end

"""
    disable_rank1_cuts!(model::VrpModel)

Disable use of *limited memory rank-1* cuts during the execution.
"""
function disable_rank1_cuts!(model::VrpModel)
    model.use_rank1_cuts = false
end

"""
    VrpGraph(model::VrpModel, nodes::Array{Int,1}, source::Int, sink::Int, multiplicity::Tuple{Int,Int})

Create a graph for pricing.

In the VRPSolver, the pricing subproblem is modeled as a *Resource Constrained Shortest Problem* (RCSP). 
VRPSolver will produce paths for the RCSP over the `VrpGraph` to generate the columns.  
Each `VrpGraph` has special `source` and `sink` vertices. The source and sink may be distinct vertices, but may also be the same vertex.

This function does not create a complete `VrpGraph`, it is necessary to create the arcs, a set of resources and to define the resource consumption intervals. 
In addition, it is necessary to map some formulation variables to the arcs of this graph.

See also: [`add_arc!`](@ref), [`add_resource!`](@ref), [`set_resource_bounds!`](@ref), [`set_arc_resource_bounds!`](@ref), [`set_arc_consumption!`](@ref), [`add_arc_var_mapping!`](@ref)

# Arguments
- `nodes::Array{Int,1}`: list of nodes (id's) of the `VrpGraph` 
- `source::Int`: id of the source node
- `sink::Int`: id of the sink node
- `multiplicity::Tuple{Int,Int}`: multiplicity of the pricing subproblem, i.e., is given lower and upper bounds, respectively, on the number of paths from this graph in a solution.  
- `standalone_filename::String`: path to a file where the graph will be exported to. 

# Examples
```julia
# let `model` be a VrpModel
# Creates a VrpGraph with 5 nodes, where at most 2 paths in this graph may appear in the solution
graph1 = VrpGraph(model, [1,2,3,4,5], 1, 5, (0,2)) # paths must start at 1 and end at 5
graph2 = VrpGraph(model, [1,2,3,4,5], 1, 1, (0,2)) # another graph whose paths must start and end at 1 (cycles)
```
"""
function VrpGraph(
    model::VrpModel,
    nodes::Array{Int,1},
    source::Int,
    sink::Int,
    multiplicity::Tuple{Int,Int},
    standalone_filename::String = "",
)
    vertices = [
        VrpVertex(
            nodes[i], i, Int[], Int[], Dict{Float64,Float64}(), Int[], Dict{Int,Any}()
        ) for i in 1:length(nodes)
    ]
    user_vertex_id_map = Dict{Int,Int}()
    for i in 1:length(nodes)
        user_vertex_id_map[nodes[i]] = i
    end

    network_source, network_sink = user_vertex_id_map[source], user_vertex_id_map[sink]

    cycle_problem = false
    if source == sink # create a sink node internally?
        network_sink = length(nodes) + 1
        push!(
            vertices,
            VrpVertex(
                maximum(nodes) + 1,
                network_sink,
                Int[],
                Int[],
                Dict{Int,Tuple{Float64,Float64}}(),
                Int[],
                Dict{Int,Any}(),
            ),
        )
        cycle_problem = true
    end

    incoming_arcs = [VrpArc[] for i in 1:length(vertices)]
    VrpGraph(
        -1,
        network_source,
        network_sink,
        vertices,
        VrpArc[],
        Dict{Int,VrpArc}(),
        incoming_arcs,
        VrpResource[],
        multiplicity,
        user_vertex_id_map,
        cycle_problem,
        Dict{Tuple{Int,Int},Tuple{Float64,Float64}}(),
        Dict{Tuple{Int,Int},Tuple{Float64,Float64}}(),
        Dict{Int,Int}(),
        nothing,
        nothing,
        Array{Array{Int,1},1}[],
        standalone_filename,
    )
end

"""
    add_graph!(model::VrpModel, graph::VrpGraph)

Add `VrpGraph` to a `VrpModel`.

This function should be called once for each graph in the model.
"""
function add_graph!(model::VrpModel, graph::VrpGraph)
    graph.id = length(model.graphs) + 1
    push!(model.graphs, graph)

    #first, create automatic resource if needed
    if isempty(filter(r -> r.is_main, graph.resources))
        res_id = add_resource!(graph; main = true, step_size = 1.0)
        graph.resources[res_id].is_automatic = true
        @warn "VRPSolver warning: An automatic main resource for Graph $(graph.id) has been created"
    end

    # second, define intervals on vertices
    for key in keys(graph.res_bounds_vertex)
        vertex, res_id = key
        res = graph.resources[res_id]
        if res.is_binary
            if vertex != graph.source_id
                graph.vertices[vertex].res_bounds[res_id] = graph.res_bounds_vertex[key]
            end
        else
            graph.vertices[vertex].res_bounds[res_id] = graph.res_bounds_vertex[key]
            for arc in graph.incoming_arcs[vertex]
                arc.res_bounds[res_id] = graph.res_bounds_vertex[key]
            end
        end
    end

    # third, define intervals on arcs (overriding, if necessary, the previous definitions)
    for key in keys(graph.res_bounds_arc)
        arc_id, res_id = key
        graph.arcs[arc_id].res_bounds[res_id] = graph.res_bounds_arc[key]
    end

    return graph
end

"""
    show(io::IO, graph::VrpGraph)

Show a VrpGraph.

"""
function show(io::IO, graph::VrpGraph)
    println(io, "===== Graph $(graph.id) =====")
    println(io, "L=$(graph.multiplicity[1]) U=$(graph.multiplicity[2])")
    println(io, "Resources")
    for r in graph.resources
        println(
            io,
            "$(r.id) $(r.is_main ? "main" : "secondary") $(r.is_disposable ? "disposable" : "non-disposable") $(r.is_binary ? "binary" : "")  $(r.is_automatic ? "automatic" : "")",
        )
    end
    bin_res_legend = if length(graph.vertices[1].res_bounds) > 0
        "\nid (bin1, [lb_bin1,ub_bin1]) (bin2, [lb_bin2,ub_bin2]) ... (consumption interval for each binary resource)"
    else
        ""
    end
    println(io, "Nodes$(bin_res_legend)")
    for v in graph.vertices
        if graph.cycle_problem && v.id == graph.sink_id
            continue
        end
        print(io, "$(v.user_id) ")
        for (i, rb) in enumerate(sort(collect(keys(v.res_bounds))))
            print(io, "($rb, [$(v.res_bounds[rb][1]),$(v.res_bounds[rb][2])]) ")
        end
        println(io)
    end
    println(
        io,
        "Arcs\n(i,j) arc_id consumption_res1 consumption_res2 ... (consumption interval for each non-binary resource) [list of mapped variables]",
    )
    net_vertex_id_map = Dict(value => key for (key, value) in graph.user_vertex_id_map) # from user_id to net_id
    for a in graph.arcs
        i = net_vertex_id_map[a.tail]
        j = if (graph.cycle_problem && a.head == graph.sink_id)
            net_vertex_id_map[graph.source_id]
        else
            net_vertex_id_map[a.head]
        end
        print(io, "($i,$j) $(a.id) ")
        for rc in a.res_consumption
            print(io, "$rc ")
        end
        print(io, "(")
        for r in 1:length(a.res_consumption)
            if !haskey(a.res_bounds, r)
                print(io, "[,] ")
            else
                (lb, ub) = a.res_bounds[r]
                print(io, "[$lb, $ub] ")
            end
        end
        print(io, ") [")
        for (id, mv) in enumerate(a.vars)
            coeff = mv[2] != 1 ? "$(mv[2])*" : ""
            if id < length(a.vars)
                print(io, string(coeff * "$(mv[1]),"))
            else
                print(io, string(coeff * "$(mv[1])"))
            end
        end
        print(io, "]\n")
    end
    println(io, "====================")
    flush(stdout)
end

"""
    add_resource!(graph::VrpGraph; main = false, binary = false, disposable = true, step_size = 0.0)

Add a resource to the VrpGraph `graph`. 

# Optional arguments 
- `main::Bool`: indicates that the resource is main or secondary. 
- `binary::Bool`: indicates that the resource is binary, i.e., its accumulated consumption can only be `0` or `1`. 
- `disposable::Bool`: indicates that the resource is disposable or non-disposable.
- `step_size::Float64`: only for main resources, this advanced parameter is used for determining the resource consumption intervals that define each bucket on the labeling algorithm during the pricing. 
if step_size is not given, it is determined automatically for main resources, based on the parameter `RCSPnumberOfBucketsPerVertex`.

# Examples
```julia
# let `graph` be a VrpGraph
r1 = add_resource!(graph) # create a secondary, disposable, non-binary resource
r2 = add_resource!(graph, main=true, disposable=false) # create a main, non-disposable, non-binary resource
```

"""
function add_resource!(
    graph::VrpGraph;
    main = false,
    binary = false,
    disposable = true,
    custom = false,
    cost_var::Union{JuMP.VariableRef,Nothing} = nothing,
    step_size = 0.0,
)
    main && binary && error("VRPSolver error: binary resource cannot be main resource")
    binary && disposable && error("VRPSolver error: binary resource cannot be disposable")
    custom &&
        (main || binary) &&
        error("VRPSolver error: custom resource cannot be main or binary")
    custom &&
        (count(r -> r.is_custom, graph.resources) > 0) &&
        error("VRPSolver error: only one custom resource can be defined")

    res_id = length(graph.resources) + 1
    push!(
        graph.resources,
        VrpResource(
            res_id, main, binary, disposable, custom, false, step_size, nothing, cost_var
        ),
    )
    if binary
        for vertex in graph.vertices
            vertex.res_bounds[res_id] = (0.0, 1.0)
        end
    end
    for arc in graph.arcs
        push!(arc.res_consumption, 0.0)
        if !binary
            arc.res_bounds[res_id] = (0.0, 0.0)
        end
    end
    return res_id
end

"""
    set_resource_bounds!(graph::VrpGraph, vertex::Int, res_id::Int, lb::Float64, ub::Float64)

Set the resource bounds to a vertex of the VrpGraph `graph`.

Defining the interval ``[lb,ub]`` for res_id at vertex is equivalent to defining the same interval for every incoming arc of vertex.

# Arguments
- `vertex::Int`: vertex id in the VrpGraph `graph`.
- `res_id::Int`: resource id in the VrpGraph `graph`.
- `lb::Float64`: lower bound on the resource consumption at the vertex.
- `ub::Float64`: upper bound on the resource consumption at the vertex.

"""
function set_resource_bounds!(
    graph::VrpGraph, vertex::Int, res_id::Int, lb::Float64, ub::Float64
)
    _check_id(res_id, 1, length(graph.resources))
    _set_resource_bounds!(graph, graph.user_vertex_id_map[vertex], res_id, lb, ub)
end

function set_resource_bounds!(graph::VrpGraph, vertex::Int, res_id::Int, lb, ub)
    set_resource_bounds!(graph, vertex, res_id, Float64(lb), Float64(ub))
end

function _set_resource_bounds!(
    graph::VrpGraph, vertex::Int, res_id::Int, lb::Float64, ub::Float64
)
    if !graph.cycle_problem && vertex == graph.source_id && (lb != 0.0 || ub != 0.0)
        @warn(
            "VRPSolver warning: Interval set for resource $(res_id) on source node (when source != sink) is ignored (by definition is [0.0,0.0])",
        )
    end
    if graph.resources[res_id].is_binary
        if (lb, ub) != (0, 0) && (lb, ub) != (0, 1) && (lb, ub) != (1, 1)
            error(
                "VRPSolver error: binary resources only supports the intervals [0,0], [0,1] and [1,1].",
            )
        end
        if !graph.resources[res_id].is_disposable &&
            (
                (vertex == graph.sink_id) ||
                (graph.cycle_problem && (vertex == graph.source_id))
            ) &&
            (lb != ub)
            error(
                "VRPSolver error: non-disposable binary resources cannot have the interval [0,1] set for the sink node",
            )
        end
    end
    graph.res_bounds_vertex[vertex, res_id] = (lb, ub) # store interval for res_id on vertex (to be defined in add_graph)
    if graph.cycle_problem && (vertex == graph.source_id) # add bounds for graph.sink_id (because the user will not set)
        _set_resource_bounds!(graph, graph.sink_id, res_id, lb, ub)
    end
end

"""
    add_elem_set_to_vertex_init_ng_neighbourhood!(model::VrpModel, graph::VrpGraph, vertex_id::Int, es_id::Int)

Add elementarity set to the vertex ng-set. 

This is an explicit way to set initial ng-neighbourhoods, which has priority over the definition with [`define_elementarity_sets_distance_matrix!`](@ref).
In fact, distance matrix is taken into account only if user-defined ng-sets are empty.

# Arguments
- `model::VrpModel`: model
- `graph::VrpGraph`: graph
- `vertex_id::Int`: vertex id
- `es_id::Int`: elementarity set id. The valid ids are ``\\{1,2,\\dots,|\\mathcal{P}^V|,|\\mathcal{P}^V|+1,|\\mathcal{P}^V|+2,\\dots,|\\mathcal{P}^V|+|\\mathcal{E}^V|\\}`` (considering automatic and additional elem. sets)
"""
function add_elem_set_to_vertex_init_ng_neighbourhood!(
    model::VrpModel, graph::VrpGraph, vertex_id::Int, es_id::Int
)
    _check_id(es_id, 1, length(model.packing_sets) + length(graph.elem_sets))
    _check_vertex_id(graph, vertex_id)
    vertex = graph.vertices[graph.user_vertex_id_map[vertex_id]]
    if !(es_id in vertex.ng_set)
        push!(vertex.ng_set, es_id)
    end
end

"""
    add_elem_set_to_arc_init_ng_neighbourhood!(model::VrpModel, graph::VrpGraph, arc_id::Int, es_id::Int)

Add elementarity set to the arc ng-set.

This is an explicit way to set initial ng-neighbourhoods, which has priority over the definition with [`define_elementarity_sets_distance_matrix!`](@ref).
In fact, distance matrix is taken into account only if user-defined ng-sets are empty.


# Arguments
- `model::VrpModel`: model
- `graph::VrpGraph`: graph
- `arc_id::Int`: arc id
- `es_id::Int`: elementarity set id. The valid ids are ``\\{1,2,\\dots,|\\mathcal{P}|,|\\mathcal{P}|+1,|\\mathcal{P}|+2,\\dots,|\\mathcal{P}|+|\\mathcal{E}|\\}`` (considering automatic and additional elem. sets)
"""
function add_elem_set_to_arc_init_ng_neighbourhood!(
    model::VrpModel, graph::VrpGraph, arc_id::Int, es_id::Int
)
    _check_id(arc_id, 1, length(graph.arcs))
    _check_id(es_id, 1, length(model.packing_sets) + length(graph.elem_sets))
    arc = graph.arcs[arc_id]
    if !(es_id in arc.ng_set)
        push!(arc.ng_set, es_id)
    end
end

"""
    add_arc_var_mapping!(graph::VrpGraph, arc_id::Int, vars::Array{Tuple{JuMP.VariableRef, Float64},1})

Define variable mapping for an existing arc.

# Arguments

- `arc_id::Int`: id of the arc
- `vars::Array{Tuple{JuMP.VariableRef, Float64},1}`: variables to be mapped to the arc. It is a set of pairs of variable and coefficient. There are additional methods where `vars` is `Array{JuMP.VariableRef,1}` and `JuMP.VariableRef` which consider all coefficients as `1.0`.

# Examples
```julia
add_arc_var_mapping!(graph, arc_id, [(x1,2.0), (x2,2.0)]) # map 2x1 and 2x2 to the arc in `graph` with id 1
add_arc_var_mapping!(graph, arc_id, [x1, x2]) # map x1 and x2
add_arc_var_mapping!(graph, arc_id, x) # map x
```
"""
function add_arc_var_mapping!(
    graph::VrpGraph, arc_id::Int, vars::Array{Tuple{JuMP.VariableRef,Float64},1}
)
    _check_id(arc_id, 1, length(graph.arcs))
    for (user_var, coeff) in vars
        if is_binary(user_var) == :Bin
            error(
                "VRPSolver error: mapping a binary variable is not allowed. Please redefine it as integer or continuous.",
            )
        end
        for (var, _) in graph.arcs[arc_id].vars
            if var == user_var
                error(
                    "VRPSolver error: variable $(user_var) is mapped more than once to arc $(arc_id) of graph $(graph.id)",
                )
            end
        end
        if coeff <= 0.0
            error(
                "VRPSolver error: mapping a variable to an arc with a non-positive coefficient is not allowed",
            )
        end
        push!(graph.arcs[arc_id].vars, (user_var, coeff))
    end
end
function add_arc_var_mapping!(graph::VrpGraph, arc_id::Int, var::JuMP.VariableRef)
    add_arc_var_mapping!(graph, arc_id, [(var, 1.0)])
end
function add_arc_var_mapping!(graph::VrpGraph, arc_id::Int, vars::Array{JuMP.VariableRef,1})
    add_arc_var_mapping!(graph, arc_id, [(var, 1.0) for var in vars])
end

"""
    add_arc!(graph::VrpGraph, tail::Int, head::Int, vars::Array{Tuple{JuMP.VariableRef, Float64},1} = Tuple{JuMP.VariableRef, Float64}[])

Add arc `(tail,head)` to `graph` and return the arc id. 

Adding parallel arcs is allowed, since they will have different identifiers in `graph`.

# Optional argument
- `vars::Array{Tuple{JuMP.VariableRef, Float64},1}`: variables to be mapped to the arc `(tail,head)`. It is a set of pairs of variable and coefficient. There are additional methods where `vars` is `Array{JuMP.VariableRef,1}` and `JuMP.VariableRef` which consider all coefficients as `1.0`. 

# Examples
```julia
# let `x1` and `x2` two decision variables
add_arc!(graph, 1, 2, [(x1,1.0), (x2,2.0)]) # add (1,2) mapped to x1 with coefficient 1 and to x2 with coefficient 2 
add_arc!(graph, 2, 1, [x1, x2]) # add (2,1) mapped to x1 and x2
arc_id = add_arc!(graph, 1, 2, x1) # add (1,2) mapped to x and get the arc id
arc_id = add_arc!(graph, 1, 2) # add (1,2) without mapped variable
```
"""
function add_arc!(
    graph::VrpGraph,
    tail::Int,
    head::Int,
    vars::Array{Tuple{JuMP.VariableRef,Float64},1} = Tuple{JuMP.VariableRef,Float64}[],
)
    arc = []
    if graph.cycle_problem && (graph.user_vertex_id_map[head] == graph.source_id)
        arc = VrpArc(
            length(graph.arcs) + 1,
            graph.user_vertex_id_map[tail],
            graph.sink_id,
            Int[],
            Int[],
            _default_consumptions(graph),
            _default_resbounds(graph),
            vars,
            Int[],
            Dict{Int,Any}(),
        )
    else
        arc = VrpArc(
            length(graph.arcs) + 1,
            graph.user_vertex_id_map[tail],
            graph.user_vertex_id_map[head],
            Int[],
            Int[],
            _default_consumptions(graph),
            _default_resbounds(graph),
            vars,
            Int[],
            Dict{Int,Any}(),
        )
    end
    push!(graph.incoming_arcs[arc.head], arc)
    push!(graph.arcs, arc)
    graph.arc_id_to_arc[arc.id] = arc
    return arc.id
end
function add_arc!(graph::VrpGraph, tail::Int, head::Int, var::JuMP.VariableRef)
    add_arc!(graph, tail, head, [(var, 1.0)])
end
function add_arc!(graph::VrpGraph, tail::Int, head::Int, vars::Array{JuMP.VariableRef,1})
    add_arc!(graph, tail, head, [(var, 1.0) for var in vars])
end

_default_consumptions(graph::VrpGraph) = [0.0 for i in 1:length(graph.resources)]

function _default_resbounds(graph::VrpGraph)
    bounds = Dict{Int,Tuple{Float64,Float64}}()
    for res in graph.resources
        if !res.is_binary
            bounds[res.id] = (0.0, 0.0)
        end
    end
    return bounds
end

"""
    set_arc_consumption!(graph::VrpGraph, arc_id::Int, res_id::Int, value::Union{Int,Float64})

Set the arc consumption for a specific resource.

# Arguments
- `graph::VrpGraph`: graph to be considered
- `arc_id::Int`: arc to be considered
- `res_id::Int`: resource id to define consumption 
- `value::Union{Int,Float64}`: consumption value which can be integer or real

# Example
```julia
set_arc_consumption!(graph, 3, 1, 2.5) # define a consumption of 2.5 for the resource 1 when passing by the arc 3 
```
"""
function set_arc_consumption!(
    graph::VrpGraph, arc_id::Int, res_id::Int, value::Union{Int,Float64}
)
    _check_id(arc_id, 1, length(graph.arcs))
    graph.resources[res_id].is_binary &&
        (value != -1 && value != 0 && value != 1) &&
        error("VRPSolver error: arc consumption for binary resources must be -1, 0, or 1")
    graph.arcs[arc_id].res_consumption[res_id] = Float64(value)
end

"""
    get_arc_consumption(graph::VrpGraph, arc_id::Int, res_id::Int)

Get the arc consumption value for a specific resource.
"""
function get_arc_consumption(graph::VrpGraph, arc_id::Int, res_id::Int)
    graph.arcs[arc_id].res_consumption[res_id]
end

"""
    set_arc_resource_bounds!(graph::VrpGraph, arc_id::Int, res_id::Int, lb::Float64, ub::Float64)

Set the resource bounds for an arc of the VrpGraph `graph`.

# Arguments
- `arc::Int`: arc id in the VrpGraph `graph`.
- `res_id::Int`: resource id in the VrpGraph `graph`.
- `lb::Float64`: lower bound on the resource consumption at the vertex.
- `ub::Float64`: upper bound on the resource consumption at the vertex.

"""
function set_arc_resource_bounds!(
    graph::VrpGraph, arc_id::Int, res_id::Int, lb::Float64, ub::Float64
)
    _check_id(res_id, 1, length(graph.resources))
    res = graph.resources[res_id]
    if res.is_binary
        error(
            "VRPSolver error: Resource bounds on arcs for binary resources are not yet implemented. Please use the vertex-based function set_resource_bounds!.",
        )
    else
        graph.res_bounds_arc[arc_id, res_id] = (lb, ub)
    end
end

function set_arc_resource_bounds!(graph::VrpGraph, arc_id::Int, res_id::Int, lb, ub)
    set_arc_resource_bounds!(graph, arc_id, res_id, Float64(lb), Float64(ub))
end

# """
#     set_arc_custom_res_params!(graph::VrpGraph, arc_id::Int, res_id::Int, value::Union{Int,Float64})

# Set the parameters attached to a given arc for a specific customized resource.

# # Arguments
# - `graph::VrpGraph`: graph to be considered
# - `arc_id::Int`: arc to be considered
# - `res_id::Int`: resource id to define parameters 
# - `values::T`: struct containing all parameter values

# # Example
# ```julia
# struct MyParams
#     first::Cint
#     secoond::Cdouble
# end

# set_arc_custom_res_params!(graph, 3, 1, MyParams(4, 2.5)) # set the parameters values 4 and 2.5 for the resource 1 when passing by the arc 3 
# ```
# """
function set_arc_custom_res_params!(
    graph::VrpGraph, arc_id::Int, res_id::Int, value::T
) where {T}
    _check_id(res_id, 1, length(graph.resources))
    !graph.resources[res_id].is_custom &&
        error("VRPSolver error: resource $(res_id) is not a custom resource")
    _check_id(arc_id, 1, length(graph.arcs))
    graph.arcs[arc_id].custom_data[res_id] = value
end

# """
#     set_vertex_custom_res_params!(graph::VrpGraph, vertex::Int, res_id::Int, value::T)

# Set the parameters attached to a given vertex for a specific customized resource.

# Defining the interval ``[lb,ub]`` for res_id at vertex is equivalent to defining the same interval for every incoming arc of vertex.

# # Arguments
# - `graph::VrpGraph`: graph to be considered
# - `vertex::Int`: vertex id in the VrpGraph `graph`.
# - `res_id::Int`: resource id in the VrpGraph `graph`.
# - `values::T`: struct containing all parameter values
# """
function set_vertex_custom_res_params!(
    graph::VrpGraph, vertex::Int, res_id::Int, value::T
) where {T}
    _check_id(res_id, 1, length(graph.resources))
    !graph.resources[res_id].is_custom &&
        error("VRPSolver error: resource $(res_id) is not a custom resource")
    v_id = graph.user_vertex_id_map[vertex]
    graph.vertices[v_id].custom_data[res_id] = value # store the value for res_id on vertex
    if graph.cycle_problem && (v_id == graph.source_id) # also store for graph.sink_id (because the user will not set)
        graph.vertices[graph.sink_id].custom_data[res_id] = value
    end
end

# """
#     set_const_custom_res_params!(graph::VrpGraph, res_id::Int, value::Union{Int,Float64})

# Set the constant parameters for a specific customized resource.

# # Arguments
# - `graph::VrpGraph`: graph to be considered
# - `res_id::Int`: resource id to define parameters 
# - `values::T`: struct containing all parameter values
# """
function set_const_custom_res_params!(graph::VrpGraph, res_id::Int, value::T) where {T}
    _check_id(res_id, 1, length(graph.resources))
    !graph.resources[res_id].is_custom &&
        error("VRPSolver error: resource $(res_id) is not a custom resource")
    graph.resources[res_id].custom_data = value
end

"""
    get_arc_set(graph::VrpGraph)

Return the set of arcs of a VrpGraph as an array of triples, where each triple contains the arc (first two elements) and its id (third element).

For example, if the function returns `[(1,2,1), (2,1,2)]`, it must be interpreted as `graph` having two arcs: `(1,2)` with id `1` and `(2,1)` with id `2`.
"""
function get_arc_set(graph::VrpGraph)
    arcs = []
    if graph.cycle_problem
        for arc in graph.arcs
            i = (arc.tail == graph.sink_id) ? graph.source_id : arc.tail
            j = (arc.head == graph.sink_id) ? graph.source_id : arc.head
            push!(arcs, (graph.vertices[i].user_id, graph.vertices[j].user_id, arc.id))
        end
    else
        arcs = [
            (graph.vertices[arc.tail].user_id, graph.vertices[arc.head].user_id, arc.id) for
            arc in graph.arcs
        ]
    end
    return arcs
end

"""
    set_arc_packing_sets!(user_model::VrpModel, collection::Array{Array{Tuple{VrpGraph,Int}, 1}, 1})

Define a collection of packing sets on arcs. For each defined packing set, VRPSolver automatically will create an equivalent elementarity set on arcs.

The collection must be a set of mutually disjoint subsets of arcs. Not all arcs need to belong to some packing set.
The index of a packing set in the array defines its packing set id.

# Examples
```julia
# let `G` be an array of VrpGraph and `model` the VrpModel
# A packing set is a list of pairs of VrpGraph and arc id
ps_1 = [(G[1],1),(G[2],1),(G[3],1)] # packing set composed by arcs of different graphs
ps_2 = [(G[2],2),(G[2],4)] # packing set composed by arcs of the same graph
ps_3 = [(G[2],3),(G[3],2)] # another packing set
set_arc_packing_sets!(model, [ps_1, ps_2, ps_3]) # passing the collection of packing sets to the model
```
"""
function set_arc_packing_sets!(
    user_model::VrpModel, collection::Array{Array{Tuple{VrpGraph,Int},1},1}
)
    if user_model.elem_sets_type != VrpSolver.NoSet
        error("VRPSolver error: Packing sets must be defined before elementarity sets")
    end
    if !isempty(user_model.cap_cuts_info)
        error(
            "VRPSolver error: Packing sets must be defined before capacity cuts separators"
        )
    end
    if !isempty(user_model.strongkpath_cuts_info)
        error(
            "VRPSolver error: Packing sets must be defined before capacity cuts separators"
        )
    end

    _reset_packing_sets(user_model)
    user_model.packing_sets = collection
    user_model.packing_sets_type = VrpSolver.ArcSet
    for ps_id in 1:length(collection)
        for (graph, arc_id) in collection[ps_id]
            _add_arc_to_packing_set(user_model, graph, arc_id, ps_id)
        end
    end
end

function _reset_packing_sets(user_model::VrpModel)
    empty!(user_model.packing_sets)
    user_model.packing_sets_type = VrpSolver.NoSet
    for graph in user_model.graphs
        for vertex in graph.vertices
            empty!(vertex.packing_sets)
            empty!(vertex.ng_set)
        end
        for arc in graph.arcs
            empty!(arc.packing_sets)
            empty!(arc.ng_set)
        end
    end
end

function _add_arc_to_packing_set(
    model::VrpModel, graph::VrpGraph, arc_id::Int, packing_set_id::Int
)
    _check_id(arc_id, 1, length(graph.arcs))
    _check_id(packing_set_id, 1, length(model.packing_sets))
    (packing_set_id in graph.arcs[arc_id].packing_sets) &&
        error("VRPSolver error: an arc cannot be added more than once to a packing set")
    push!(graph.arcs[arc_id].packing_sets, packing_set_id)
end

"""
    set_additional_arc_elementarity_sets!(user_model::VrpModel, collection::Array{Tuple{VrpGraph,Array{Int,1}}, 1})

Define an additional collection of elementarity sets on arcs.

The collection must be a set of mutually disjoint subsets of arcs. Not all arcs need to belong to some elementarity set.
The index ``i`` of a elementarity set in the array defines its elementarity set id as ``|\\mathcal{P}|+i`` because there are ``|\\mathcal{P}|`` 
automatic elementarity sets (one for each packing set) created with [`set_arc_packing_sets!`](@ref).

# Examples
```julia
# let `G` be an array of VrpGraph and `model` the VrpModel
# An elementarity set is a pair of VrpGraph and a list of arc ids
es_1 = (G[1],[1,2,4]) # elem. set on graph G[1] composed by the arcs 1, 2, and 4  
es_2 = (G[2],[2,3,5,9]) # elem. set on G[2] with 4 arcs
es_3 = (G[2],[1,4,6,7]) # another elem. set on G[2]
set_additional_arc_elementarity_sets!(model, [es_1, es_2, es_3]) # passing the collection of elem. sets to the model
```
"""
function set_additional_arc_elementarity_sets!(
    user_model::VrpModel, collection::Array{Tuple{VrpGraph,Array{Int,1}},1}
)
    if user_model.packing_sets_type == VrpSolver.VertexSet
        error("Vertex packing sets and arc elementarity sets are not compatible")
    end
    _reset_elem_sets(user_model)
    user_model.elem_sets_type = VrpSolver.ArcSet
    for es_id in 1:length(collection)
        graph, arc_ids = collection[es_id]
        push!(graph.elem_sets, arc_ids)
        for arc_id in arc_ids
            _add_arc_to_elementarity_set(user_model, graph, arc_id, es_id)
        end
    end
    return length(user_model.packing_sets)
end

function _reset_elem_sets(user_model::VrpModel)
    user_model.elem_sets_type = VrpSolver.NoSet
    for graph in user_model.graphs
        empty!(graph.elem_sets)
        graph.es_dist_matrix = nothing
        for vertex in graph.vertices
            empty!(vertex.elem_sets)
            empty!(vertex.ng_set)
        end
        for arc in graph.arcs
            empty!(arc.elem_sets)
            empty!(arc.ng_set)
        end
    end
end

function _add_arc_to_elementarity_set(
    model::VrpModel, graph::VrpGraph, arc_id::Int, es_id::Int
)
    _check_id(arc_id, 1, length(graph.arcs))
    _check_id(es_id, 1, length(graph.elem_sets))
    (es_id in graph.arcs[arc_id].elems_sets) && error(
        "VRPSolver error: an arc cannot be added more than once to an elementarity set"
    )
    push!(graph.arcs[arc_id].elem_sets, es_id)
end

"""
    set_vertex_packing_sets!(user_model::VrpModel, collection::Array{Array{Tuple{VrpGraph,Int}, 1}, 1})

Define a collection of packing sets on vertices. For each defined packing set, VRPSolver automatically will create an equivalent elementarity set on vertices.

The collection must be a set of subsets of vertices. The subsets are not necessarily disjoint. 
Not all vertices need to belong to some packing set. The index of a packing set in the array defines its packing set id.

# Examples
```julia
# let `G` be an array of VrpGraph and `model` the VrpModel
# A packing set is a list of pairs of VrpGraph and vertex id
ps_1 = [(G[1],1),(G[2],1),(G[3],1)] # packing set composed by vertices of different graphs
ps_2 = [(G[2],2),(G[2],4)] # packing set composed by vertices of the same graph
ps_3 = [(G[2],3),(G[3],2)] # another packing set
set_vertex_packing_sets!(model, [ps_1, ps_2, ps_3]) # passing the collection of packing sets to the model
```
"""
function set_vertex_packing_sets!(
    user_model::VrpModel,
    collection::Array{Array{Tuple{VrpGraph,Int},1},1},
    define_covering_sets::Bool = false,
)
    if user_model.elem_sets_type != NoSet
        error("VRPSolver error:  Packing sets must be defined before elementarity sets")
    end
    if !isempty(user_model.cap_cuts_info)
        error(
            "VRPSolver error: Packing sets must be defined before capacity cuts separators"
        )
    end
    if !isempty(user_model.strongkpath_cuts_info)
        error(
            "VRPSolver error: Packing sets must be defined before capacity cuts separators"
        )
    end

    _reset_packing_sets(user_model)
    user_model.packing_sets = collection
    user_model.packing_sets_type = VertexSet
    user_model.define_covering_sets = define_covering_sets
    n = length(collection)
    for ps_id in 1:n
        for (graph, vertex_id) in collection[ps_id]
            _add_vertex_to_packing_set(user_model, graph, vertex_id, ps_id)
        end
    end

    _has_overlapping_sets(user_model) && return nothing

    # function to compute the packing set pair connected by the arc of a pair
    # (graph, arc) where vectices not associated to packing sets are assigned to
    # the dummy packing set id n+1
    function get_packing_set_pair(gr_arc::Tuple{VrpGraph,VrpArc})::Tuple{Int,Int}
        graph = gr_arc[1]
        arc = gr_arc[2]
        head = if !isempty(graph.vertices[arc.head].packing_sets)
            first(graph.vertices[arc.head].packing_sets)
        else
            n + 1
        end
        tail = if !isempty(graph.vertices[arc.tail].packing_sets)
            first(graph.vertices[arc.tail].packing_sets)
        else
            n + 1
        end
        if head < tail
            return head, tail
        else
            return tail, head
        end
    end

    # build data structures to access lists of arcs by packing set pairs
    # and by mapped variables (to be used next)
    user_model.arcs_by_packing_set_pairs = [Tuple{VrpGraph,VrpArc}[] for _ in 1:n, _ in 1:n]
    mapped_arcs_by_vars = Dict{JuMP.VariableRef,Array{Tuple{VrpGraph,VrpArc},1}}()
    for graph in user_model.graphs
        for arc in graph.arcs
            head, tail = get_packing_set_pair((graph, arc))
            if tail <= n
                push!(user_model.arcs_by_packing_set_pairs[head, tail], (graph, arc))
            end
            for (var, _) in arc.vars
                if haskey(mapped_arcs_by_vars, var)
                    push!(mapped_arcs_by_vars[var], (graph, arc))
                else
                    mapped_arcs_by_vars[var] = [(graph, arc)]
                end
            end
        end
    end

    # determine, for each packing set pairs, which arcs are not covered by
    # appropriate variables (variables mapped to arcs connecting the same packing
    # set pair, in any direction)
    for (_, gr_arcs) in mapped_arcs_by_vars
        head, tail = get_packing_set_pair(gr_arcs[1])
        appropriate = (head != tail) && (tail <= n)
        for gr_arc in gr_arcs[2:end]
            if (head, tail) != get_packing_set_pair(gr_arc)
                appropriate = false
            end
        end
        if appropriate
            ps_gr_arcs = user_model.arcs_by_packing_set_pairs[head, tail]
            filter!(x -> !(x in gr_arcs), ps_gr_arcs)
        end
    end
end

function _add_vertex_to_packing_set(
    model::VrpModel, graph::VrpGraph, vertex_id::Int, packing_set_id::Int
)
    _check_vertex_id(graph, vertex_id)
    _check_id(packing_set_id, 1, length(model.packing_sets))
    vertexAlgId = graph.user_vertex_id_map[vertex_id]
    (packing_set_id in graph.vertices[vertexAlgId].packing_sets) &&
        error("VRPSolver error: a vertex cannot be added more than once to a packing set")
    push!(graph.vertices[vertexAlgId].packing_sets, packing_set_id)
end

"""
    set_additional_vertex_elementarity_sets!(user_model::VrpModel, collection::Array{Tuple{VrpGraph,Array{Int,1}}, 1})

Define an additional collection of elementarity sets on vertices.

The collection must be a set of mutually disjoint subsets of vertices. Not all vertices need to belong to some elementarity set.
The index ``i`` of a elementarity set in the array defines its elementarity set id as ``|\\mathcal{P}^V|+i`` because there are ``|\\mathcal{P}^V|`` 
automatic elementarity sets (one for each packing set) created with [`set_vertex_packing_sets!`](@ref).

# Examples
```julia
# let `G` be an array of VrpGraph and `model` the VrpModel
# An elementarity set is a pair of VrpGraph and a list of vertex ids
es_1 = (G[1],[1,2,4]) # elem. set on graph G[1] composed by the vertices 1, 2, and 4  
es_2 = (G[2],[2,3,5,9]) # elem. set on G[2] with 4 vertices
es_3 = (G[2],[1,4,6,7]) # another elem. set on G[2]
set_additional_vertex_elementarity_sets!(model, [es_1, es_2, es_3]) # passing the collection of elem. sets to the model
```
"""
function set_additional_vertex_elementarity_sets!(
    user_model::VrpModel, collection::Array{Tuple{VrpGraph,Array{Int,1}},1}
)
    if user_model.packing_sets_type == VrpSolver.ArcSet
        error(
            "VRPSolver error: Arc packing sets and vertex elementarity sets are not compatible",
        )
    end
    _reset_elem_sets(user_model)
    user_model.elem_sets_type = VrpSolver.VertexSet
    n = length(collection)
    for es_id in 1:n
        (graph, vertex_ids) = collection[es_id]
        push!(graph.elem_sets, vertex_ids)
        for vertex_id in vertex_ids
            _add_vertex_to_elem_set(user_model, graph, vertex_id, es_id)
        end
    end
end

function _add_vertex_to_elem_set(
    model::VrpModel, graph::VrpGraph, vertex_id::Int, es_id::Int
)
    _check_vertex_id(graph, vertex_id)
    _check_id(es_id, 1, length(graph.elem_sets))
    vertexAlgId = graph.user_vertex_id_map[vertex_id]
    (es_id in graph.vertices[vertexAlgId].elem_sets) && error(
        "VRPSolver error: a vertex cannot be added more than once to an elementarity set",
    )
    push!(graph.vertices[vertexAlgId].elem_sets, es_id)
end

"""
    set_branching_priority!(model::VrpModel, var_container_name::String, priority::Int)

Set the branching priority for a decision variable. 

Branching is performed on a variable with priority ``k`` only when there is no possible branching for variables with priority higher than ``k``.

# Arguments
- `var_container_name::String`: JuMP decision variable name.
- `priority::Int`: the priority of the variable, which is higher when this value is increased.

# Example
```julia
model = VrpModel()
@variable(model.formulation, x[i in 1:10] >= 0, Int)
@variable(model.formulation, y[i in 1:20] >= 0, Int)
... 
set_branching_priority!(model, "x", 2) # x has higher priority
set_branching_priority!(model, "y", 1) # than y
```
"""
function set_branching_priority!(model::VrpModel, var_container_name::String, priority::Int)
    model.branching_priorities[var_container_name] = priority
end

"""
    function set_branching_priority!(model::VrpModel, exps_family, name::String, priority::Int)

Set the branching priority for a set of JuMP expressions.

# Arguments
- `exps_family`: set of expressions created with JuMP macro @expression 
- `name::String`: JuMP expression name.
- `priority::Int`: the priority of the expressions, which is higher when this value is increased.

"""
function set_branching_priority!(model::VrpModel, exps_family, name::String, priority::Int)
    push!(model.branching_exp_families, (exps_family, name, priority))
end

"""
    function set_branching_priority!(model::VrpModel, exp::JuMP.GenericAffExpr{Float64,JuMP.VariableRef}, name::String, priority::Int)

Set the branching priority for a JuMP expression.

# Arguments
- `exp`: expression created with JuMP macro @expression 
- `name::String`: JuMP expression name.
- `priority::Int`: the priority of the expression, which is higher when this value is increased.

# Example
```julia
model = VrpModel()
@variable(model.formulation, x[i in 1:10] >= 0, Int)
@expression(model.formulation, exp, sum(x[i] for i in 1:5))
... 
set_branching_priority!(model, "x", 2) # x has higher priority
set_branching_priority!(model, exp, "exp", 1) # than exp
```
"""
function set_branching_priority!(
    model::VrpModel,
    exp::JuMP.GenericAffExpr{Float64,JuMP.VariableRef},
    name::String,
    priority::Int,
)
    push!(model.branching_exps, (exp, name, priority))
end

"""
    function enable_resource_consumption_branching!(model::VrpModel, priority::Int)

Enable the accumulated resource consumption branching. 
Given a packing set ``S \\in \\mathcal{P}``, a main resource ``r \\in R^k_M``, and a certain threshold value ``t^*``: in the left child make ``u_{a,r}=t^*``,
for all ``a \\in S``; in the right child make ``l_{a,r}=t^*``, for all ``a \\in S``.
The packing set, the main resource, and the threshold are automatically chosen by VRPSolver during the execution.
The branching is not likely to be complete, in the sense that some fractional ``\\lambda`` solutions can not be eliminated by it. 
However, it does not increase the pricing difficulty and it may work well in practice, postponing (and even avoiding) the use of a branching that makes pricing harder.

# Arguments
- `priority::Int`: the branching priority, which is higher when this value is increased.

# Example
```julia
model = VrpModel()
... 
enable_resource_consumption_branching!(model, 1) 
```
"""
function enable_resource_consumption_branching!(model::VrpModel, priority::Int)
    model.branching_priorities["_res_cons_branching_"] = priority
end

"""
    function enable_packset_ryanfoster_branching!(model::VrpModel, priority::Int)

Enable the Ryan and Foster branching.
Given two distinct packing sets ``S`` and ``S'`` in ``\\mathcal{P}``, let ``P(S,S') \\subseteq P`` be the subset of the paths that contain arcs in both ``S`` and ``S'``.
The branching is over the value of ``\\sum_{p\\in P(S,S')} \\lambda_p``, either 0 or 1.
The packing sets are automatically chosen by VRPSolver during the execution.
This branching is still to be avoided if possible, because it makes the pricing harder. However, using that scheme leads to more balanced search trees.

# Arguments
- `priority::Int`: the branching priority, which is higher when this value is increased.

# Example
```julia
model = VrpModel()
... 
enable_packset_ryanfoster_branching!(model, 2) 
```
"""
function enable_packset_ryanfoster_branching!(model::VrpModel, priority::Int)
    model.branching_priorities["_ryanfoster_branching_"] = priority
end

"""
    function add_permanent_ryanfoster_constraint!(model::VrpModel, firstPackSetId::Integer, secondPackSetId::Integer, together::Bool)

Set a permanent Ryan and Foster constraint for a pair of packing sets.

# Arguments
- `firstPackSetId::Int`: first packing set id.
- `secondPackSetId::Int`: second packing set id.
- `together::Bool`: say whether they should be together in the same path or should not be together.

# Example
```julia
model = VrpModel()
... 
add_permanent_ryanfoster_constraint!(model, 1, 2, true) 
add_permanent_ryanfoster_constraint!(model, 3, 4, false) 
```
"""
function add_permanent_ryanfoster_constraint!(
    model::VrpModel, firstPackSetId::Integer, secondPackSetId::Integer, together::Bool
)
    _check_id(firstPackSetId, 1, length(model.packing_sets))
    _check_id(secondPackSetId, 1, length(model.packing_sets))
    (firstPackSetId == secondPackSetId) && error(
        "VRPSolver error: packing sets must be different in a permanent Ryan and Foster constraint",
    )
    push!(model.ryanfoster_constraints, (firstPackSetId, secondPackSetId, together))
end

"""
    define_elementarity_sets_distance_matrix!(model::VrpModel, graph::VrpGraph, matrix::Array{Array{Float64,1},1})

Define distance matrix between the elementarity sets for a specific graph.

This is a convenient way of defining the initial ng-sets for pricing over ``\\mathcal{E}``-ng-paths: the ng-set of each vertex (arc) is defined as being the ng-size elementarity sets 
that are closer to its own elementarity set, according to the given distance matrix. The element `matrix[i][j]` represents the distance between the elementarity set ``i`` to the 
elementarity set ``j``.
An index of an elementarity set is defined during its creation: automatic elementarity sets (one for each packing set) has the indexes
 ``\\{1,2,\\dots,|\\mathcal{P}|\\}``, whereas the additional elementarity sets has the indexes ``\\{|\\mathcal{P}|+1,|\\mathcal{P}|+2,\\dots\\}``
  (these indexes consider ``|\\mathcal{P}^V|`` instead of ``|\\mathcal{P}|`` for packing sets on vertices). 
To define the distance between two elementarity sets, 
it is recommended to consider some metric which involves all elements (vertices or arcs) of both elementarity sets.

"""
function define_elementarity_sets_distance_matrix!(
    model::VrpModel, graph::VrpGraph, matrix::Array{Array{Float64,1},1}
)
    n = length(model.packing_sets) + length(graph.elem_sets)
    size(matrix, 1) != n && error("VRPSolver error: invalid matrix dimension")
    for i in 1:n
        size(matrix[i], 1) != n && error("VRPSolver error: invalid matrix dimension")
    end
    graph.es_dist_matrix = matrix
end

"""
    add_cut_callback!(user_model::VrpModel, callback::Any, constr_name::String)

Add user cut callback for separation during the execution.

Cuts must be added through the function [`add_dynamic_constr!`](@ref).

# Arguments
- `model::VrpModel`: model to be added the cut callback.
- `callback::Any`: function with the separation algorithm 
- `constr_name::String`: a nome for the cut callback 

# Example
```julia
# let model be a VrpModel and x a set of variables for edges
function edge_ub_callback()
   for (i,j) in E
     e = (i,j)
      if i != 0 && get_value(model.optimizer, x[e]) > 1.001
         println("Adding edge ub cut for e = ", e)
         add_dynamic_constr!(model.optimizer, [x[e]], [1.0], <=, 1.0, "edge_ub")
      end
   end 
end 
add_cut_callback!(model, edge_ub_callback, "edge_ub")
```
"""
function add_cut_callback!(user_model::VrpModel, callback::Any, constr_name::String)
    user_model.callbacks[constr_name] = callback
end

"""
    add_capacity_cut_separator!(model::VrpModel, demands::Array{Tuple{Array{Tuple{VrpGraph,Int}, 1},Float64},1}, capacity::Float64)

Define a *rounded capacity cut* (RCC) separator over a collection of packing sets defined on vertices.
RCC separators cannot be used if the packings sets are defined on arcs.

# Arguments
- `model::VrpModel`: model to be added the RCC separator.
- `demands::Array{Tuple{Array{Tuple{VrpGraph,Int}, 1},Float64},1}`: array of pairs of packing set and demand.
- `capacity::Float64`: capacity considered in the RCC separator.

# Examples
```julia
# Let PS be an array of `n` packing sets in vertices, where PS[i] is the i-th packing set
# Let d[i] be the demand associated with the i-th packing set
# Let `model` be a VrpModel and Q the capacity to be used in the RCC
add_capacity_cut_separator!(model, [(PS[i], d[i]) for i in 1:n], Q) # add a RCC separator
```
"""
function add_capacity_cut_separator!(
    model::VrpModel,
    demands::Array{Tuple{Array{Tuple{VrpGraph,Int},1},Float64},1},
    capacity::Float64,
    two_path_cuts_res_id::Int = -1,
)
    _has_overlapping_sets(model) && error(
        "VRPSolver error: capacity cut separators are incompatible with overlapping (packing/elementarity) sets",
    )

    for (ps_set, _) in demands
        !(ps_set in model.packing_sets) && error(
            "VRPSolver error: collection that is not a packing set was used in a capacity cut separator." *
            " Only the packing set collections can be used for add_capacity_cut_separator!",
        )
    end

    # create and map variables to all uncovered arcs connecting packing set pairs
    id_demands = [0 for _ in 1:length(model.packing_sets)]
    if !isempty(model.arcs_by_packing_set_pairs)
        for (ps_set, d) in demands
            ps_id = findall(x -> x == ps_set, model.packing_sets)
            id_demands[ps_id[1]] = Int(d)
        end

        num_missing_arcs = 0
        uncovered = Tuple{Int,Int}[]
        dims_psp = size(model.arcs_by_packing_set_pairs)
        arcs_by_psp = model.arcs_by_packing_set_pairs
        for head in 1:dims_psp[1], tail in (head + 1):dims_psp[2]
            if (id_demands[head] > 0) &&
                (id_demands[tail] > 0) &&
                !isempty(arcs_by_psp[head, tail])
                push!(uncovered, (head, tail))
                num_missing_arcs += length(arcs_by_psp[head, tail])
            end
        end
        if length(uncovered) > 0
            @warn "VrpSolver: adding $(length(uncovered)) internal variables mapping to $num_missing_arcs arcs for use by capacity cuts"
        end
        if num_missing_arcs > 0
            @variable(model.formulation, RCCsepX[ps_pair in uncovered], Int)
            for (head, tail) in uncovered
                for (graph, arc) in arcs_by_psp[head, tail]
                    add_arc_var_mapping!(graph, arc.id, RCCsepX[(head, tail)])
                end
                arcs_by_psp[head, tail] = Array{Tuple{VrpGraph,VrpArc},1}(undef, 0)
            end
        end
    end

    push!(
        model.cap_cuts_info,
        CapacityCutInfo(id_demands, Int(capacity), two_path_cuts_res_id),
    )
end

"""
    add_strongkpath_cut_separator!(model::VrpModel, demands::Array{Tuple{Array{Tuple{VrpGraph,Int}, 1},Float64},1}, capacity::Float64)

Define a *strong k-path* (SKP) separator over a collection of packing sets defined on vertices.
SKP separators cannot be used if the packings sets are defined on arcs.

# Arguments
- `model::VrpModel`: model to be added the RCC separator.
- `demands::Array{Tuple{Array{Tuple{VrpGraph,Int}, 1},Float64},1}`: array of pairs of packing set and demand.
- `capacity::Float64`: capacity considered in the RCC separator.

# Examples
```julia
# Let PS be an array of `n` packing sets in vertices, where PS[i] is the i-th packing set
# Let d[i] be the demand associated with the i-th packing set
# Let `model` be a VrpModel and Q the capacity to be used in the RCC
add_capacity_cut_separator!(model, [(PS[i], d[i]) for i in 1:n], Q) # add a RCC separator
```
"""
function add_strongkpath_cut_separator!(
    model::VrpModel,
    demands::Array{Tuple{Array{Tuple{VrpGraph,Int},1},Float64},1},
    capacity::Float64,
)
    for (ps_set, _) in demands
        !(ps_set in model.packing_sets) && error(
            "VRPSolver error: collection that is not a packing set was used in a strong k-path cut separator." *
            " Only the packing set collections can be used for add_strongkpath_cut_separator!",
        )
    end

    _has_overlapping_sets(model) && error(
        "VRPSolver error: strong k-path cut separators are incompatible with overlapping (packing/elementarity) sets",
    )

    # need to define covering sets for strong k-path cuts
    model.define_covering_sets = true

    # create and map variables to all uncovered arcs connecting packing set pairs
    id_demands = [0 for i in 1:length(model.packing_sets)]
    if !isempty(model.arcs_by_packing_set_pairs)
        for (ps_set, d) in demands
            ps_id = findall(x -> x == ps_set, model.packing_sets)
            id_demands[ps_id[1]] = Int(d)
        end

        num_missing_arcs = 0
        uncovered = Tuple{Int,Int}[]
        dims_psp = size(model.arcs_by_packing_set_pairs)
        arcs_by_psp = model.arcs_by_packing_set_pairs
        for head in 1:dims_psp[1], tail in (head + 1):dims_psp[2]
            if (id_demands[head] > 0) &&
                (id_demands[tail] > 0) &&
                !isempty(arcs_by_psp[head, tail])
                push!(uncovered, (head, tail))
                num_missing_arcs += length(arcs_by_psp[head, tail])
            end
        end
        if length(uncovered) > 0
            @warn "VrpSolver: adding $(length(uncovered)) internal variables mapping to $num_missing_arcs arcs for use by strong k-path cuts"
        end
        if num_missing_arcs > 0
            @variable(model.formulation, RCCsepX[ps_pair in uncovered], Int)
            for (head, tail) in uncovered
                for (graph, arc) in arcs_by_psp[head, tail]
                    add_arc_var_mapping!(graph, arc.id, RCCsepX[(head, tail)])
                end
                arcs_by_psp[head, tail] = Array{Tuple{VrpGraph,VrpArc},1}(undef, 0)
            end
        end
    end

    push!(model.strongkpath_cuts_info, CapacityCutInfo(id_demands, Int(capacity), -1))
end

function _should_use_meta_solver(user_model::VrpModel)
    features_A = [
        "Strong k-path cuts", "Ryan and Foster branching", "Resource consumption branching"
    ]
    features_B = ["Custom resource", "Overlapping (packing/elementarity) sets"]

    found_A = []
    if _using_resource_consumption_branching(user_model)
        push!(found_A, "Resource consumption branching")
    end
    if _using_ryan_and_foster_branching(user_model)
        push!(found_A, "Ryan and Foster branching")
    end
    if _using_strong_kpath_cuts(user_model)
        push!(found_A, "Strong k-path cuts")
    end

    found_B = []
    if _has_custom_resource(user_model)
        push!(found_B, "Custom resource")
    end
    if _has_overlapping_sets(user_model)
        push!(found_B, "Overlapping (packing/elementarity) sets")
    end

    if !isempty(found_A) && !isempty(found_B)
        error("""
VRPSolver error: Model cannot include features from both sets simultaneously.
This is a current limitation of VRPSolver and will be addressed in future versions.

Set A: $(join(features_A, ", "))
Set B: $(join(features_B, ", "))

Features from Set A found in your model: $(join(found_A, ", "))
Features from Set B found in your model: $(join(found_B, ", "))
""")
    end

    return isempty(found_A)
end
