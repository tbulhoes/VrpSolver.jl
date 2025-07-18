module VrpSolver

using JuMP
using MathOptInterface: MathOptInterface
const MOI = MathOptInterface

import Base.show
using Printf

include("bapcod_interface.jl")

export VrpModel, VrpGraph, VrpOptimizer
export add_resource!,
    set_resource_bounds!,
    add_elem_set_to_vertex_init_ng_neighbourhood!,
    add_arc!,
    add_arc_var_mapping!,
    set_arc_consumption!,
    set_arc_resource_bounds!,
    get_arc_consumption,
    add_elem_set_to_arc_init_ng_neighbourhood!,
    get_arc_set,
    set_arc_packing_sets!,
    set_vertex_packing_sets!,
    set_additional_arc_elementarity_sets!,
    set_additional_vertex_elementarity_sets!,
    add_graph!,
    define_elementarity_sets_distance_matrix!,
    add_capacity_cut_separator!,
    add_strongkpath_cut_separator!,
    set_branching_priority!,
    enable_rank1_cuts!,
    disable_rank1_cuts!,
    enable_resource_consumption_branching!,
    enable_packset_ryanfoster_branching!,
    set_cutoff!,
    get_objective_value,
    get_value,
    get_values,
    get_number_of_positive_paths,
    get_path_arcs,
    add_cut_callback!,
    add_dynamic_constr!,
    show,
    get_complete_formulation,
    print_enum_paths,
    add_permanent_ryanfoster_constraint!

@enum SetType NoSet = 0 ArcSet = 1 VertexSet = 2

mutable struct VrpVertex
    user_id::Int
    id::Int
    packing_set::Int
    elem_set::Int
    res_bounds::Dict{Int,Tuple{Float64,Float64}} #only for binary resources
    ng_set::Array{Int,1}
end

mutable struct VrpArc
    id::Int
    tail::Int
    head::Int
    packing_set::Int
    elem_set::Int
    res_consumption::Array{Float64,1}
    res_bounds::Dict{Int,Tuple{Float64,Float64}} #only for non-binary resources
    vars::Array{Tuple{JuMP.VariableRef,Float64},1}
    ng_set::Array{Int,1}
end

mutable struct VrpResource
    id::Int
    is_main::Bool
    is_binary::Bool
    is_disposable::Bool
    is_automatic::Bool
    step_size::Float64
end

mutable struct VrpGraph
    id::Int
    source_id::Int
    sink_id::Int
    vertices::Vector{VrpVertex}
    arcs::Vector{VrpArc}
    arc_id_to_arc::Dict{Int,VrpArc}
    incoming_arcs::Vector{Vector{VrpArc}}
    resources::Array{VrpResource,1}
    multiplicity::Tuple{Int,Int}
    user_vertex_id_map::Dict{Int,Int}
    cycle_problem::Bool # true when source=sink for the user (vertices[sink_id] is only an internal vertex)
    res_bounds_vertex::Dict{Tuple{Int,Int},Tuple{Float64,Float64}} # stores all intervals defined by the user on vertices
    res_bounds_arc::Dict{Tuple{Int,Int},Tuple{Float64,Float64}} # stores all intervals defined by the user on arcs
    arc_bapcod_id_to_id::Dict{Int,Int}
    net
    es_dist_matrix
    elem_sets::Array{Array{Int,1},1}
    standalone_filename::String
end

mutable struct CapacityCutInfo
    demands::Array{Int,1}
    capacity::Int
    two_path_cuts_res_id::Int
end

mutable struct VrpModel
    formulation::JuMP.Model
    graphs::Array{VrpGraph,1}
    packing_sets::Array{Array{Tuple{VrpGraph,Int},1},1}
    packing_sets_type::SetType
    elem_sets_type::SetType
    define_covering_sets::Bool
    branching_priorities::Dict{String,Int}
    branching_exp_families::Array{Any,1}
    branching_exps::Array{Any,1}
    use_rank1_cuts::Bool
    optimizer::Any
    callbacks::Dict{String,Any}
    cap_cuts_info::Array{CapacityCutInfo,1}
    strongkpath_cuts_info::Array{CapacityCutInfo,1}
    arcs_by_packing_set_pairs::Array{Array{Tuple{VrpGraph,VrpArc}},2}
    ryanfoster_constraints::Vector{Tuple{Integer,Integer,Bool}} # (firstPackSetId,secondPackSetId,together)
end

mutable struct DynamicConstrInfo
    vars::Array{JuMP.VariableRef,1}
    coeffs::Array{Float64,1}
    sense::Function
    rhs::Float64
end

mutable struct CallbackInfo
    constr_name::String
    to_add_constrs::Array{DynamicConstrInfo,1}
end

struct OptimizerColsInfo
    uservar_to_colids::Dict{JuMP.VariableRef,Array{Int,1}}
    uservar_to_problem_type::Dict{JuMP.VariableRef,Symbol}
    cols_uservar::Vector{JuMP.VariableRef}
    cols_names::Vector{String}
    cols_problems::Vector{Tuple{Symbol,Int}}
    cols_lbs::Vector{Float64}
    cols_ubs::Vector{Float64}
    cols_costs::Vector{Float64}
    cols_types::Vector{Char}
end

struct SpSol
    graph_id::Int
    multiplicity::Int
    user_vars_in_sol::Dict{JuMP.VariableRef,Float64}
    arc_seq::Vector{VrpArc}
end

mutable struct VrpOptimizer
    user_model::VrpModel
    bapcod_model
    param_file::String
    instance_name::String
    callbacks::Dict{String,CallbackInfo}
    spsols_in_sol::Vector{SpSol}
    unmapped_vars_in_sol::Dict{JuMP.VariableRef,Float64}
    sol_defined::Bool
    integer_objective::Bool
    initUB::Float64
    stats::Dict{} # execution statistics
    baptreedot_file::String
    optimizer_cols_info::OptimizerColsInfo
end

contains(p, s) = !isnothing(findnext(s, p, 1))

function check_id(id::Int, min_id::Int, max_id::Int)
    (id < min_id || id > max_id) && error("Unknown id $id.")
end

function check_vertex_id(graph::VrpGraph, id::Int)
    !(id in keys(graph.user_vertex_id_map)) && error("Unknown vertex id $id.")
end

function resource_id_in_bapcod(res::VrpResource, graph::VrpGraph)
    if !res.is_binary
        return res.id
    end
    b_id = 1
    for res2 in graph.resources
        if res2.id == res.id
            return b_id
        elseif res2.is_binary
            b_id += 1
        end
    end
end

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
        NoSet,
        NoSet,
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
        VrpVertex(nodes[i], i, -1, -1, Dict{Float64,Float64}(), Int[]) for
        i in 1:length(nodes)
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
                -1,
                -1,
                Dict{Int,Tuple{Float64,Float64}}(),
                Int[],
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
    graph::VrpGraph; main = false, binary = false, disposable = true, step_size = 0.0
)
    main && binary && error("VRPSolver error: binary resource cannot be main resource")
    binary && disposable && error("VRPSolver error: binary resource cannot be disposable")

    res_id = length(graph.resources) + 1
    push!(graph.resources, VrpResource(res_id, main, binary, disposable, false, step_size))
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

function set_resource_bounds_aux!(
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
        set_resource_bounds_aux!(graph, graph.sink_id, res_id, lb, ub)
    end
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
    check_id(res_id, 1, length(graph.resources))
    set_resource_bounds_aux!(graph, graph.user_vertex_id_map[vertex], res_id, lb, ub)
end

function set_resource_bounds!(graph::VrpGraph, vertex::Int, res_id::Int, lb, ub)
    set_resource_bounds!(graph, vertex, res_id, Float64(lb), Float64(ub))
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
    check_id(es_id, 1, length(model.packing_sets) + length(graph.elem_sets))
    check_vertex_id(graph, vertex_id)
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
    check_id(arc_id, 1, length(graph.arcs))
    check_id(es_id, 1, length(model.packing_sets) + length(graph.elem_sets))
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
- `vars::Array{Tuple{JuMP.VariableRef, Float64},1}`: variables to be mapped to the arc. It is a set of pairs of variable and coefficient. There are extensions where `vars` is `Array{JuMP.VariableRef,1}` and `JuMP.VariableRef` which consider all coefficients as `1.0`.

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
    check_id(arc_id, 1, length(graph.arcs))
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

default_consumptions(graph::VrpGraph) = [0.0 for i in 1:length(graph.resources)]

function default_resbounds(graph::VrpGraph)
    bounds = Dict{Int,Tuple{Float64,Float64}}()
    for res in graph.resources
        if !res.is_binary
            bounds[res.id] = (0.0, 0.0)
        end
    end
    return bounds
end

"""
    add_arc!(graph::VrpGraph, tail::Int, head::Int, vars::Array{Tuple{JuMP.VariableRef, Float64},1} = Tuple{JuMP.VariableRef, Float64}[])

Add arc `(tail,head)` to `graph` and return the arc id. 

Adding parallel arcs is allowed, since they will have different identifiers in `graph`.

# Optional argument
- `vars::Array{Tuple{JuMP.VariableRef, Float64},1}`: variables to be mapped to the arc `(tail,head)`. It is a set of pairs of variable and coefficient. There are extensions where `vars` is `Array{JuMP.VariableRef,1}` and `JuMP.VariableRef` which consider all coefficients as `1.0`. 

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
            -1,
            -1,
            default_consumptions(graph),
            default_resbounds(graph),
            vars,
            Int[],
        )
    else
        arc = VrpArc(
            length(graph.arcs) + 1,
            graph.user_vertex_id_map[tail],
            graph.user_vertex_id_map[head],
            -1,
            -1,
            default_consumptions(graph),
            default_resbounds(graph),
            vars,
            Int[],
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
    check_id(arc_id, 1, length(graph.arcs))
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
    check_id(res_id, 1, length(graph.resources))
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

function is_preprocessed_arc(graph::VrpGraph, arc::VrpArc)
    if !graph.cycle_problem && (arc.head == graph.source_id || arc.tail == graph.sink_id)
        return true
    end
    return false
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

function reset_packing_sets(user_model::VrpModel)
    empty!(user_model.packing_sets)
    user_model.packing_sets_type = NoSet
    for graph in user_model.graphs
        for vertex in graph.vertices
            vertex.packing_set = -1
            empty!(vertex.ng_set)
        end
        for arc in graph.arcs
            arc.packing_set = -1
            empty!(arc.ng_set)
        end
    end
end

function reset_elem_sets(user_model::VrpModel)
    user_model.elem_sets_type = NoSet
    for graph in user_model.graphs
        empty!(graph.elem_sets)
        graph.es_dist_matrix = nothing
        for vertex in graph.vertices
            vertex.elem_set = -1
            empty!(vertex.ng_set)
        end
        for arc in graph.arcs
            arc.elem_set = -1
            empty!(arc.ng_set)
        end
    end
end

function add_arc_to_packing_set(
    model::VrpModel, graph::VrpGraph, arc_id::Int, packing_set_id::Int
)
    check_id(arc_id, 1, length(graph.arcs))
    check_id(packing_set_id, 1, length(model.packing_sets))
    (graph.arcs[arc_id].packing_set != -1) &&
        error("VRPSolver error: an arc cannot belong to more than 1 packing set")
    graph.arcs[arc_id].packing_set = packing_set_id
    (graph.arcs[arc_id].elem_set != -1) &&
        error("VRPSolver error: an arc cannot belong to more than 1 elementarity set")
end

function add_arc_to_elementarity_set(
    model::VrpModel, graph::VrpGraph, arc_id::Int, es_id::Int
)
    check_id(arc_id, 1, length(graph.arcs))
    check_id(es_id, 1, length(graph.elem_sets) + length(model.packing_sets))
    (graph.arcs[arc_id].elem_set != -1) &&
        error("VRPSolver error: an arc cannot belong to more than 1 elementarity set")
    (graph.arcs[arc_id].packing_set != -1) &&
        error("VRPSolver error: an arc cannot belong to more than 1 elementarity set")
    graph.arcs[arc_id].elem_set = es_id
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
    if user_model.elem_sets_type != NoSet
        error("VRPSolver error: Packing sets cannot be defined after elementarity sets")
    end
    reset_packing_sets(user_model)
    user_model.packing_sets = collection
    user_model.packing_sets_type = ArcSet
    for ps_id in 1:length(collection)
        for (graph, arc_id) in collection[ps_id]
            add_arc_to_packing_set(user_model, graph, arc_id, ps_id)
        end
    end
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
    if user_model.packing_sets_type == VertexSet
        error("Vertex packing sets and arc elementarity sets are not compatible")
    end
    reset_elem_sets(user_model)
    user_model.elem_sets_type = ArcSet
    for es_id in 1:length(collection)
        graph, arc_ids = collection[es_id]
        push!(graph.elem_sets, arc_ids)
        for arc_id in arc_ids
            add_arc_to_elementarity_set(user_model, graph, arc_id, es_id)
        end
    end
    return length(user_model.packing_sets)
end

function add_vertex_to_packing_set(
    model::VrpModel, graph::VrpGraph, vertex_id::Int, packing_set_id::Int
)
    check_vertex_id(graph, vertex_id)
    check_id(packing_set_id, 1, length(model.packing_sets))
    vertexAlgId = graph.user_vertex_id_map[vertex_id]
    (graph.vertices[vertexAlgId].packing_set != -1) &&
        error("VRPSolver error: a vertex cannot belong to more than 1 packing set")
    graph.vertices[vertexAlgId].packing_set = packing_set_id
    (graph.vertices[vertexAlgId].elem_set != -1) &&
        error("VRPSolver error: a vertex cannot belong to more than 1 elementarity set")
end

function add_vertex_to_elem_set(
    model::VrpModel, graph::VrpGraph, vertex_id::Int, es_id::Int
)
    check_vertex_id(graph, vertex_id)
    check_id(es_id, 1, length(graph.elem_sets) + length(model.packing_sets))
    vertexAlgId = graph.user_vertex_id_map[vertex_id]
    (graph.vertices[vertexAlgId].elem_set != -1) &&
        error("VRPSolver error: a vertex cannot belong to more than 1 elementarity set")
    graph.vertices[vertexAlgId].elem_set = es_id
    (graph.vertices[vertexAlgId].packing_set != -1) &&
        error("VRPSolver error: a vertex cannot belong to more than 1 elementarity set")
end

"""
    set_vertex_packing_sets!(user_model::VrpModel, collection::Array{Array{Tuple{VrpGraph,Int}, 1}, 1})

Define a collection of packing sets on vertices. For each defined packing set, VRPSolver automatically will create an equivalent elementarity set on vertices.

The collection must be a set of mutually disjoint subsets of vertices. Not all vertices need to belong to some packing set.
The index of a packing set in the array defines its packing set id.

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
        error("VRPSolver error: Packing sets cannot be defined after elementarity sets")
    end
    reset_packing_sets(user_model)
    user_model.packing_sets = collection
    user_model.packing_sets_type = VertexSet
    user_model.define_covering_sets = define_covering_sets
    n = length(collection)
    for ps_id in 1:n
        for (graph, vertex_id) in collection[ps_id]
            add_vertex_to_packing_set(user_model, graph, vertex_id, ps_id)
        end
    end

    # function to compute the packing set pair connected by the arc of a pair
    # (graph, arc) where vectices not associated to packing sets are assigned to
    # the dummy packing set id n+1
    function get_packing_set_pair(gr_arc::Tuple{VrpGraph,VrpArc})::Tuple{Int,Int}
        graph = gr_arc[1]
        arc = gr_arc[2]
        head = graph.vertices[arc.head].packing_set
        head = (head < 1) ? n + 1 : head
        tail = graph.vertices[arc.tail].packing_set
        tail = (tail < 1) ? n + 1 : tail
        if head < tail
            return head, tail
        else
            return tail, head
        end
    end

    # build data structures to access lists of arcs by packing set pairs
    # and by mapped variables (to be used next)
    user_model.arcs_by_packing_set_pairs = [Tuple{VrpGraph,VrpArc}[] for i in 1:n, j in 1:n]
    mapped_arcs_by_vars = Dict{JuMP.VariableRef,Array{Tuple{VrpGraph,VrpArc},1}}()
    for graph in user_model.graphs
        for arc in graph.arcs
            head, tail = get_packing_set_pair((graph, arc))
            if tail <= n
                push!(user_model.arcs_by_packing_set_pairs[head, tail], (graph, arc))
            end
            for (var, val) in arc.vars
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
    if user_model.packing_sets_type == ArcSet
        error(
            "VRPSolver error: Arc packing sets and vertex elementarity sets are not compatible",
        )
    end
    reset_elem_sets(user_model)
    user_model.elem_sets_type = VertexSet
    n = length(collection)
    for es_id in 1:n
        (graph, vertex_ids) = collection[es_id]
        push!(graph.elem_sets, vertex_ids)
        for vertex_id in vertex_ids
            add_vertex_to_elem_set(user_model, graph, vertex_id, es_id)
        end
    end
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
    check_id(firstPackSetId, 1, length(model.packing_sets))
    check_id(secondPackSetId, 1, length(model.packing_sets))
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
    size(matrix, 1) != n && error("VRPSolver error: wrong matrix dimension")
    for i in 1:n
        size(matrix[i], 1) != n && error("VRPSolver error: wrong matrix dimension")
    end
    graph.es_dist_matrix = matrix
end

function extract_user_var_to_graphs(user_model::VrpModel)
    var_to_graphs = Dict{JuMP.VariableRef,Array{Int,1}}()
    for graph_id in 1:length(user_model.graphs)
        for arc_id in 1:length(user_model.graphs[graph_id].arcs)
            for (var, _) in user_model.graphs[graph_id].arcs[arc_id].vars
                if !haskey(var_to_graphs, var)
                    var_to_graphs[var] = [graph_id]
                elseif !(graph_id in var_to_graphs[var])
                    push!(var_to_graphs[var], graph_id)
                end
            end
        end
    end
    return var_to_graphs
end

function split_var_name(var)
    var_name = name(var)
    var_container_name = contains(var_name, "[") ? split(var_name, "[")[1] : var_name
    var_id = contains(var_name, "[") ? split(split(var_name, "[")[2], "]")[1] : "1"
    var_id = replace(var_id, r"\(" => s"")
    var_id = replace(var_id, r"\)" => s"")
    var_id = tuple([parse(Int64, x) for x in split(var_id, ",") if x != ""]...)
    return (var_container_name, var_id)
end

function generate_pricing_networks(
    user_model::VrpModel, bapcod_model, optimizer_cols_info::OptimizerColsInfo
)
    isempty(user_model.graphs) && error("VRPSolver error: no graph defined")

    c_register_subproblems(
        bapcod_model, [(spid, :DW_SP) for spid in 0:(length(user_model.graphs) - 1)]
    )

    for graph in user_model.graphs
        nbPackSets = length(user_model.packing_sets)
        nbElemSets = nbPackSets + length(graph.elem_sets)
        nbCovSets = 0
        if user_model.define_covering_sets
            nbCovSets = nbPackSets
        end
        c_net_ptr = new_network!(
            bapcod_model,
            graph.id - 1,
            :DW_SP,
            length(graph.vertices),
            nbPackSets,
            nbElemSets,
            nbCovSets,
        )
        wbcr_set_source(c_net_ptr, graph.source_id - 1)
        wbcr_set_sink(c_net_ptr, graph.sink_id - 1)

        graph.net = c_net_ptr
        #resources and vertices
        for resource in graph.resources
            if resource.is_binary
                if !resource.is_disposable
                    wbcr_set_special_as_nondisposable_resource(
                        c_net_ptr, resource_id_in_bapcod(resource, graph) - 1
                    )
                end
            else
                wbcr_new_resource(c_net_ptr, resource.id - 1)
                if resource.is_main
                    wbcr_set_as_main_resource(c_net_ptr, resource.id - 1, 0.0)
                end
                if !resource.is_disposable
                    wbcr_set_as_nondisposable_resource(c_net_ptr, resource.id - 1)
                end
            end
            for vertex in graph.vertices
                if resource.is_binary
                    res_seq_id = resource_id_in_bapcod(resource, graph)
                    wbcr_set_vertex_special_consumption_lb(
                        c_net_ptr,
                        vertex.id - 1,
                        res_seq_id - 1,
                        vertex.res_bounds[resource.id][1],
                    )
                    wbcr_set_vertex_special_consumption_ub(
                        c_net_ptr,
                        vertex.id - 1,
                        res_seq_id - 1,
                        vertex.res_bounds[resource.id][2],
                    )
                else
                    if !haskey(vertex.res_bounds, resource.id)
                        wbcr_set_vertex_consumption_lb(
                            c_net_ptr, vertex.id - 1, resource.id - 1, -1e12
                        )
                        wbcr_set_vertex_consumption_ub(
                            c_net_ptr, vertex.id - 1, resource.id - 1, 1e12
                        )
                    else
                        wbcr_set_vertex_consumption_lb(
                            c_net_ptr,
                            vertex.id - 1,
                            resource.id - 1,
                            vertex.res_bounds[resource.id][1],
                        )
                        wbcr_set_vertex_consumption_ub(
                            c_net_ptr,
                            vertex.id - 1,
                            resource.id - 1,
                            vertex.res_bounds[resource.id][2],
                        )
                    end
                end
            end
        end

        #adding vertices to packing and elem sets
        for vertex in graph.vertices
            if vertex.packing_set != -1
                wbcr_add_vertex_to_packing_set(
                    c_net_ptr, vertex.id - 1, vertex.packing_set - 1
                )
                wbcr_attach_elementarity_set_to_node(
                    c_net_ptr, vertex.id - 1, vertex.packing_set - 1
                )
                if user_model.define_covering_sets
                    wbcr_add_vertex_to_covering_set(
                        c_net_ptr, vertex.id - 1, vertex.packing_set - 1
                    )
                end
            elseif vertex.elem_set != -1
                wbcr_attach_elementarity_set_to_node(
                    c_net_ptr, vertex.id - 1, vertex.elem_set - 1
                )
            end
            #defining the neighbourhood of the vertex
            for es_id in vertex.ng_set
                wbcr_add_vertex_to_mem_of_elementarity_set(
                    c_net_ptr, vertex.id - 1, es_id - 1
                )
            end
        end
        #arcs
        for arc in graph.arcs
            if is_preprocessed_arc(graph, arc)
                continue
            end
            arc_bapcod_id = wbcr_new_arc(c_net_ptr, arc.tail - 1, arc.head - 1, 0.0)
            graph.arc_bapcod_id_to_id[arc_bapcod_id] = arc.id
            for resource in graph.resources
                if resource.is_binary
                    res_seq_id = resource_id_in_bapcod(resource, graph)
                    wbcr_set_edge_special_consumption_value(
                        c_net_ptr,
                        arc_bapcod_id,
                        res_seq_id - 1,
                        arc.res_consumption[resource.id],
                    )
                else
                    wbcr_set_edge_consumption_value(
                        c_net_ptr,
                        arc_bapcod_id,
                        resource.id - 1,
                        arc.res_consumption[resource.id],
                    )
                    wbcr_set_arc_consumption_lb(
                        c_net_ptr,
                        arc_bapcod_id,
                        resource.id - 1,
                        arc.res_bounds[resource.id][1],
                    )
                    wbcr_set_arc_consumption_ub(
                        c_net_ptr,
                        arc_bapcod_id,
                        resource.id - 1,
                        arc.res_bounds[resource.id][2],
                    )
                end
            end
            # adding arc to packing_set
            if arc.packing_set != -1
                wbcr_add_edge_to_packing_set(c_net_ptr, arc_bapcod_id, arc.packing_set - 1)
                wbcr_attach_elementarity_set_to_edge(
                    c_net_ptr, arc_bapcod_id, arc.packing_set - 1
                )
            elseif arc.elem_set != -1
                wbcr_attach_elementarity_set_to_edge(
                    c_net_ptr, arc_bapcod_id, nbPackSets + arc.elem_set - 1
                )
            end
            # defining the neighbourhood of the arc
            for es_id in arc.ng_set
                wbcr_add_arc_to_mem_of_elementarity_set(c_net_ptr, arc_bapcod_id, es_id - 1)
            end
            #arc variables
            for (user_var, coeff) in arc.vars
                for colid in optimizer_cols_info.uservar_to_colids[user_var]
                    if optimizer_cols_info.cols_problems[colid + 1][2] == graph.id - 1
                        wbcr_attach_bcvar_to_arc(
                            c_net_ptr, arc_bapcod_id, bapcod_model, colid, coeff
                        )
                    end
                end
            end
        end

        # elem sets distance matrix
        if !isnothing(graph.es_dist_matrix)
            wbcr_set_elementarity_sets_distance_matrix(
                c_net_ptr, graph.es_dist_matrix, nbElemSets
            )
        end

        #ryan and foster constraints
        for (ps1, ps2, tog) in user_model.ryanfoster_constraints
            wbcr_add_permanent_ryanfoster_constraint(c_net_ptr, ps1 - 1, ps2 - 1, tog)
        end

        new_oracle!(
            c_net_ptr, bapcod_model, :DW_SP, graph.id - 1, graph.standalone_filename
        )
    end

    c_set_sp_multiplicities(
        bapcod_model,
        [
            (spid, :DW_SP, user_model.graphs[spid + 1].multiplicity...) for
            spid in 0:(length(user_model.graphs) - 1)
        ],
    )
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
    for (ps_set, d) in demands
        !(ps_set in model.packing_sets) && error(
            "Collection that is not a packing set was used in a capacity cut separator." *
            " Only the packing set collections can be used for add_capacity_cut_separator",
        )
    end

    # create and map variables to all uncovered arcs connecting packing set pairs
    if !isempty(model.arcs_by_packing_set_pairs)
        id_demands = [0 for i in 1:length(model.packing_sets)]
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

function add_capacity_cut_separators_to_optimizer(optimizer::VrpOptimizer)
    user_model = optimizer.user_model
    for cap_cut_info in user_model.cap_cuts_info
        demands = [0 for i in 1:length(user_model.packing_sets)]
        for (ps_set, d) in cap_cut_info.demands
            ps_id = findall(x -> x == ps_set, user_model.packing_sets)
            demands[ps_id[1]] = Int(d)
        end
        add_rcsp_capacity_cuts!(
            optimizer.formulation,
            Int(cap_cut_info.capacity),
            demands;
            is_facultative = false,
            root_priority_level = 3.0,
            two_path_cuts_res_id = cap_cut_info.two_path_cuts_res_id,
        )
    end
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
            "VRPSolver error: Collection that is not a packing set was used in a strong k-path cut separator." *
            " Only the packing set collections can be used for add_strongkpath_cut_separator",
        )
    end

    # create and map variables to all uncovered arcs connecting packing set pairs
    if !isempty(model.arcs_by_packing_set_pairs)
        id_demands = [0 for i in 1:length(model.packing_sets)]
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

    push!(model.strongkpath_cuts_info, CapacityCutInfo(demands, capacity, -1))
end

function add_strongkpath_cut_separators_to_optimizer(optimizer::VrpOptimizer)
    user_model = optimizer.user_model
    for cap_cut_info in user_model.strongkpath_cuts_info
        demands = [0 for i in 1:length(user_model.packing_sets)]
        for (ps_set, d) in cap_cut_info.demands
            ps_id = findall(x -> x == ps_set, user_model.packing_sets)
            demands[ps_id[1]] = Int(d)
        end
        add_rcsp_strongkpath_cuts!(
            optimizer.formulation,
            Int(cap_cut_info.capacity),
            demands;
            is_facultative = false,
            root_priority_level = 1.0,
        )
    end
end

function extract_optimizer_cols_info(user_model::VrpModel)
    user_var_to_graphs = extract_user_var_to_graphs(user_model)
    user_form = user_model.formulation
    uservar_to_colids = Dict{JuMP.VariableRef,Vector{Int}}()
    uservar_to_problem_type = Dict{JuMP.VariableRef,Symbol}()
    cols_uservar = JuMP.VariableRef[]
    cols_names = String[]
    cols_problems = Tuple{Symbol,Int}[]
    cols_lbs = Float64[]
    cols_ubs = Float64[]
    cols_costs = Float64[]
    cols_types = Char[]
    obj = objective_function(user_form)
    nextcolid = 0
    user_vars = all_variables(user_form)
    for user_var in user_vars
        colsids = Int[]
        if !haskey(user_var_to_graphs, user_var) # unmapped variable
            push!(colsids, nextcolid)
            push!(cols_problems, (:DW_MASTER, 0))
            push!(cols_lbs, has_lower_bound(user_var) ? lower_bound(user_var) : -Inf)
            push!(cols_ubs, has_upper_bound(user_var) ? upper_bound(user_var) : Inf)
            push!(cols_names, name(user_var))
            push!(cols_uservar, user_var)
            push!(cols_costs, get(obj.terms, user_var, 0.0))
            push!(
                cols_types,
                if is_integer(user_var)
                    'I'
                else
                    'C'
                end,
            )
            uservar_to_problem_type[user_var] = :DW_MASTER
            nextcolid += 1
        else
            for graph_id in user_var_to_graphs[user_var] # mapped variable
                push!(colsids, nextcolid)
                push!(cols_problems, (:DW_SP, graph_id - 1))
                push!(cols_names, name(user_var) * "_$(graph_id-1)")
                push!(cols_lbs, 0.0)
                push!(cols_ubs, Inf)
                if has_lower_bound(user_var) && lower_bound(user_var) != 0.0
                    error(
                        "VRPSolver error: lower bound of mapped variables must be zero. Add explicit master constraints to impose a different lower bound",
                    )
                end
                if has_upper_bound(user_var)
                    error(
                        "VRPSolver error: upper bound of mapped variables must be infinity. Add explicit master constraints to impose a different upper bound",
                    )
                end
                push!(cols_uservar, user_var)
                push!(cols_costs, get(obj.terms, user_var, 0.0))
                push!(
                    cols_types,
                    if is_integer(user_var)
                        'I'
                    else
                        'C'
                    end,
                )
                nextcolid += 1
            end
            uservar_to_problem_type[user_var] = :DW_SP
        end
        if !isempty(colsids)
            uservar_to_colids[user_var] = colsids
        end
    end
    return OptimizerColsInfo(
        uservar_to_colids,
        uservar_to_problem_type,
        cols_uservar,
        cols_names,
        cols_problems,
        cols_lbs,
        cols_ubs,
        cols_costs,
        cols_types,
    )
end

function build_optimizer_vars_and_constrs(
    user_model::VrpModel, bapcod_model_ptr, optimizer_cols_info::OptimizerColsInfo
)
    user_form = user_model.formulation
    nconstrs = sum(
        values(num_constraints(user_form; count_variable_in_set_constraints = false))
    )
    ncols = length(optimizer_cols_info.cols_problems)
    init_model!(bapcod_model_ptr, nconstrs, ncols)

    # creating variables
    user_vars = all_variables(user_form)
    optimizer_vars = Tuple{Symbol,Int,Symbol,Int}[]
    user_vars = all_variables(user_form)
    for user_var in user_vars
        for colid in optimizer_cols_info.uservar_to_colids[user_var]
            push!(
                optimizer_vars,
                (
                    Symbol(optimizer_cols_info.cols_names[colid + 1]),
                    colid,
                    optimizer_cols_info.cols_problems[colid + 1]...,
                ),
            )
        end
    end
    c_register_vars(
        bapcod_model_ptr,
        optimizer_cols_info.cols_lbs,
        optimizer_cols_info.cols_ubs,
        optimizer_cols_info.cols_costs,
        optimizer_cols_info.cols_types,
        optimizer_vars,
    )

    constr_idx = 0
    clbs = Float64[]
    cubs = Float64[]
    constrs = Tuple{Symbol,Int,Symbol,Int}[]
    constrs_refs = ConstraintRef[]
    for set_type in (MOI.LessThan{Float64}, MOI.GreaterThan{Float64}, MOI.EqualTo{Float64})
        append!(constrs_refs, all_constraints(user_form, AffExpr, set_type))
    end

    for constr_ref in constrs_refs
        set = JuMP.constraint_object(constr_ref).set
        if set isa MathOptInterface.LessThan
            push!(clbs, -Inf)
            push!(cubs, set.upper)
        elseif set isa MathOptInterface.GreaterThan
            push!(clbs, set.lower)
            push!(cubs, Inf)
        elseif set isa MathOptInterface.EqualTo
            push!(clbs, set.value)
            push!(cubs, set.value)
        else
            error("VRPSolver error: cannot recognize constraint sense.")
        end
        push!(constrs, (Symbol(name(constr_ref)), constr_idx, :DW_MASTER, -1))
        constr_idx += 1
    end
    nb_cols = sum(length(v) for (_, v) in optimizer_cols_info.uservar_to_colids)
    rows_id_vec = [Int[] for _ in 1:nb_cols]
    nonzeros_vec = [Float64[] for _ in 1:nb_cols]
    for (constr_idx, constr_ref) in enumerate(constrs_refs)
        terms = JuMP.constraint_object(constr_ref).func.terms
        for (user_var, coeff) in terms
            for id in optimizer_cols_info.uservar_to_colids[user_var]
                push!(rows_id_vec[id + 1], constr_idx - 1)
                push!(nonzeros_vec[id + 1], coeff)
            end
        end
    end
    nonzeros = vcat(nonzeros_vec...)
    rows_id = vcat(rows_id_vec...)
    starts = Int[]
    pos = 0
    for v in nonzeros_vec
        push!(starts, pos)
        pos += length(v)
    end
    push!(starts, pos)
    c_register_cstrs(bapcod_model_ptr, starts, rows_id, nonzeros, clbs, cubs, constrs)
end

function has_integer_objective(user_model::VrpModel, optimizer_cols_info::OptimizerColsInfo)
    user_form = user_model.formulation
    user_vars = all_variables(user_form)
    obj = objective_function(user_form)
    for user_var in user_vars
        # checking if the coefficient in the objective function is integral
        var_cost = get(obj.terms, user_var, 0.0)
        if modf(var_cost)[1] != 0.0
            return false
        end
        # checking if unmapped variables are integer
        if optimizer_cols_info.uservar_to_problem_type[user_var] == :DW_MASTER &&
            !is_integer(user_var)
            return false
        end
    end

    # if some variable is mapped with a non-integer coeff, return false
    for graph in user_model.graphs
        for arc in graph.arcs
            for (_, coeff) in arc.vars
                if modf(coeff)[1] != 0.0
                    return false
                end
            end
        end
    end

    return true
end

function set_branching_priorities_in_optimizer(
    user_model::VrpModel, bapcod_model_ptr, optimizer_cols_info::OptimizerColsInfo
)
    user_vars = all_variables(user_model.formulation)

    # defining the priority of single variables
    cols_priority_infos = Tuple{Symbol,Symbol,Int,Float64}[]
    for (var_container_name, priority) in user_model.branching_priorities
        for user_var in user_vars
            length(optimizer_cols_info.uservar_to_colids[user_var]) > 1 && continue # we will branch on an expression aggregating the mapped variables
            (container_name, _) = split_var_name(user_var)
            if var_container_name == container_name
                for colid in optimizer_cols_info.uservar_to_colids[user_var]
                    push!(
                        cols_priority_infos,
                        (
                            Symbol(optimizer_cols_info.cols_names[colid + 1]),
                            optimizer_cols_info.cols_problems[colid + 1]...,
                            priority,
                        ),
                    )
                end
            end
        end
    end
    if !isempty(cols_priority_infos)
        c_vars_branching_priorities(bapcod_model_ptr, cols_priority_infos)
    end

    ##########
    # now we add branching on expressions
    ##########

    # first, we define the expressions corresponding to the aggregation of mapped variables.
    # This is an automatic step, meaning that the user does not need to mind about this aggregation.
    for (var_container_name, priority) in user_model.branching_priorities
        expr_array_id = nothing
        nb_exprs_in_array = 0
        for user_var in user_vars
            (container_name, _) = split_var_name(user_var)
            if var_container_name == container_name
                colsids = get(optimizer_cols_info.uservar_to_colids, user_var, [])
                length(colsids) <= 1 && continue
                if isnothing(expr_array_id)
                    expr_array_id = register_branching_expression(
                        bapcod_model_ptr, "aggr_" * container_name, Float64(priority)
                    )
                end
                coeffs = [1.0 for _ in 1:length(colsids)]
                add_branching_expression(
                    bapcod_model_ptr, expr_array_id, nb_exprs_in_array + 1, colsids, coeffs
                )
                nb_exprs_in_array += 1
            end
        end
    end

    # second, we add the expression families defined by the user
    for (exp_family, name, priority) in user_model.branching_exp_families
        exp_array_id = register_branching_expression(
            bapcod_model_ptr, name, Float64(priority)
        )
        for index in keys(exp_family)
            expr = exp_family[index]
            if expr.constant != 0.0
                error(
                    "VRPSolver error: constant part of branching expression must be equal to zero",
                )
            end
            colsids, coeffs = Int[], Float64[]
            for user_var_idx in 1:length(expr.terms.keys)
                user_var = expr.terms.keys[user_var_idx]
                coeff = expr.terms.vals[user_var_idx]
                for colid in optimizer_cols_info.uservar_to_colids[user_var]
                    push!(colsids, colid)
                    push!(coeffs, coeff)
                end
            end
            if isempty(colsids)
                @warn "VRPSolver warning: branching expression with name $(name) and index $(index.I[1]) is empty"
            else
                add_branching_expression(
                    bapcod_model_ptr, exp_array_id, index.I[1], colsids, coeffs
                )
            end
        end
    end

    # finally, we add the single expressions defined the user
    for (expr, name, priority) in user_model.branching_exps
        if expr.constant != 0.0
            error(
                "VRPSolver error: constant part of branching expression must be equal to zero",
            )
        end
        colsids, coeffs = Int[], Float64[]
        for user_var_idx in 1:length(expr.terms.keys)
            user_var = expr.terms.keys[user_var_idx]
            coeff = expr.terms.vals[user_var_idx]
            for colid in optimizer_cols_info.uservar_to_colids[user_var]
                push!(colsids, colid)
                push!(coeffs, coeff)
            end
        end
        if isempty(colsids)
            @warn "VRPSolver warning: branching expression with name $(name) is empty"
        else
            exp_array_id = register_branching_expression(
                bapcod_model_ptr, name, Float64(priority)
            )
            add_branching_expression(bapcod_model_ptr, exp_array_id, (1,), colsids, coeffs)
        end
    end
end

"""
    VrpOptimizer(user_model::VrpModel, param_file::String, instance_name = ""; baptreedot="BaPTree.dot")

Build an optimizer for a VrpModel.

# Arguments
- `model::VrpModel`: model to be considered
- `param_file::String`: path for the VRPSolver parameters file

# Optional arguments
- `instance_name::String`: the instance name to be shown in the results line (line with execution statistics).
- `baptreedot::String`: path to the file to output the BaP Tree in dot format.
"""
function VrpOptimizer(
    user_model::VrpModel,
    param_file::String,
    instance_name = "";
    baptreedot = "BaPTree.dot",
    fixed_params = [
        "",
        "--MaxNbOfStagesInColGenProcedure",
        "3",
        "--colGenSubProbSolMode",
        "3",
        "--MipSolverMultiThread",
        "1",
        "--ApplyStrongBranchingEvaluation",
        "true",
        "-t",
        baptreedot,
    ],
)
    optimizer_cols_info = extract_optimizer_cols_info(user_model)
    integer_objective = has_integer_objective(user_model, optimizer_cols_info)

    bapcod_model_ptr = new!(
        param_file,
        true,
        integer_objective,
        false, # CHECK
        11,
        fixed_params,
    )

    #creating optimizer formulation
    build_optimizer_vars_and_constrs(user_model, bapcod_model_ptr, optimizer_cols_info)

    #creating bapcod networks
    generate_pricing_networks(user_model, bapcod_model_ptr, optimizer_cols_info)

    #branching priorities
    set_branching_priorities_in_optimizer(user_model, bapcod_model_ptr, optimizer_cols_info)

    for rcc in user_model.cap_cuts_info
        wbcr_add_generic_capacity_cut(bapcod_model_ptr, rcc.capacity, rcc.demands)
    end

    if user_model.use_rank1_cuts
        wbc_add_generic_lim_mem_one_cut(bapcod_model_ptr)
    end

    if haskey(user_model.branching_priorities, "_res_cons_branching_")
        priority = user_model.branching_priorities["_res_cons_branching_"]
        c_add_elemset_resource_cons_branching(bapcod_model_ptr, priority)
    end
    if haskey(user_model.branching_priorities, "_ryanfoster_branching_")
        priority = user_model.branching_priorities["_ryanfoster_branching_"]
        c_add_packset_ryan_and_foster_branching(bapcod_model_ptr, priority)
    end

    optimizer = VrpOptimizer(
        user_model,
        bapcod_model_ptr,
        param_file,
        instance_name,
        Dict{String,CallbackInfo}(),
        Tuple{Int,Dict{JuMP.VariableRef,Float64}}[],
        Dict{JuMP.VariableRef,Float64}(),
        false,
        integer_objective,
        -1,
        Dict(),
        baptreedot,
        optimizer_cols_info,
    )
    user_model.optimizer = optimizer

    # CHECK
    # add_strongkpath_cut_separators_to_optimizer(optimizer)

    return optimizer
end

"""
    optimize!(optimizer::VrpOptimizer)

Solve a VRPSolver problem.

It returns a pair with the status of the execution and a flag indicating whether a solution was found, respectively.

# Example
```julia
# let model be a VrpModel
optimizer = VrpOptimizer(model, "path_to_config/config.cfg")
(status, solution_found) = optimize!(optimizer)
```
"""
function JuMP.optimize!(optimizer::VrpOptimizer)
    # Save access to the optimizer in the JuMP model to allow for getting variable values without
    # having to pass it as an argument.
    optimizer.user_model.formulation.ext[:optimizer] = optimizer

    if !isempty(optimizer.user_model.callbacks)
        register_lazycb!(optimizer.bapcod_model, optimizer)
        for constr_name in keys(optimizer.user_model.callbacks)
            optimizer.callbacks[constr_name] = CallbackInfo(
                constr_name, DynamicConstrInfo[]
            )
        end
    end

    sol_ptr = new_sol!()
    status = c_optimize(optimizer.bapcod_model, sol_ptr)
    has_solution = register_solutions(optimizer, sol_ptr)

    println(
        "statistics_cols: instance & :Optimal & cutoff & :bcRecRootDb & :bcTimeRootEval & :bcCountNodeProc & :bcRecBestDb & :bcRecBestInc & :bcTimeMain \\\\",
    )
    print("statistics: $(optimizer.instance_name) & ")
    print("$(status == :Optimal ? 1 : 0) & ")
    print("$(optimizer.initUB) & ")
    @printf("%.2f & ", getstatistic(optimizer.bapcod_model, :bcRecRootDb))
    @printf("%.2f & ", getstatistic(optimizer.bapcod_model, :bcTimeRootEval) / 100)
    print("$(getstatistic(optimizer.bapcod_model, :bcCountNodeProc)) & ")
    @printf("%.2f & ", getstatistic(optimizer.bapcod_model, :bcRecBestDb))
    if has_solution
        if optimizer.integer_objective
            print(
                "$(Int(floor(getstatistic(optimizer.bapcod_model, :bcRecBestInc) + 0.5))) & ",
            )
        else
            @printf("%.2f & ", getstatistic(optimizer.bapcod_model, :bcRecBestInc))
        end
        optimizer.stats[:bcRecBestInc] = getstatistic(optimizer.bapcod_model, :bcRecBestInc)
    else
        print("-- & ")
    end
    @printf("%.2f \\\\\n", getstatistic(optimizer.bapcod_model, :bcTimeMain) / 100)
    flush(stdout)

    for statkey in bcvalues
        optimizer.stats[statkey] = getstatistic(optimizer.bapcod_model, statkey)
    end
    for statkey in bctimers
        optimizer.stats[statkey] = getstatistic(optimizer.bapcod_model, statkey)
    end
    for statkey in bccounters
        optimizer.stats[statkey] = getstatistic(optimizer.bapcod_model, statkey)
    end

    return status, has_solution
end

function register_solutions(optimizer::VrpOptimizer, bapcodsol)
    bapcod_model = optimizer.bapcod_model
    user_vars = all_variables(optimizer.user_model.formulation)
    spsols_in_sol = SpSol[]
    optimizer_cols_info = optimizer.optimizer_cols_info
    unmapped_vars_in_sol = Dict{JuMP.VariableRef,Float64}()

    # emptying the solution stored in the optimizer
    empty!(optimizer.unmapped_vars_in_sol)
    empty!(optimizer.spsols_in_sol)

    # getting master solution
    if c_start(bapcodsol, bapcod_model) != 1
        optimizer.sol_defined = false
        return false
    end

    # getting unmapped vars values
    for user_var in user_vars
        if optimizer_cols_info.uservar_to_problem_type[user_var] == :DW_MASTER
            colid = optimizer_cols_info.uservar_to_colids[user_var][1]
            unmapped_vars_in_sol[user_var] = c_getValueOfVar(bapcod_model, bapcodsol, colid)
        end
    end
    optimizer.unmapped_vars_in_sol = unmapped_vars_in_sol

    # retrieving subproblem sols
    status = c_next(bapcodsol) # skip master solution
    while status == 1
        subproblem_id = Int(c_getProblemFirstId(bapcodsol))
        mult = c_getMultiplicity(bapcodsol)
        graph = optimizer.user_model.graphs[subproblem_id + 1]
        arcs_ids = [
            graph.arc_bapcod_id_to_id[arc_bap_id] for arc_bap_id in c_getArcs(bapcodsol)
        ]
        spsol = SpSol(
            graph.id,
            mult,
            get_path_uservar_map(arcs_ids, graph),
            [graph.arcs[arc_id] for arc_id in arcs_ids],
        )
        push!(spsols_in_sol, spsol)
        status = c_next(bapcodsol)
    end
    optimizer.spsols_in_sol = spsols_in_sol

    optimizer.sol_defined = true

    return true
end

"""
    set_cutoff!(optimizer::VrpOptimizer, ub::Float64)

Set an upper bound (primal bound) for the execution.

"""
function set_cutoff!(optimizer::VrpOptimizer, ub::Float64)
    set_obj_ub!(optimizer.bapcod_model, ub)
    optimizer.initUB = ub
    absolute_ub = abs(ub)
    if (absolute_ub < 10000.0)
        set_art_cost_value!(optimizer.bapcod_model, 10000.0)
    else
        set_art_cost_value!(optimizer.bapcod_model, absolute_ub)
    end
end
set_cutoff!(optimizer::VrpOptimizer, ub::Int) = set_cutoff!(optimizer, Float64(ub))

"""
    get_objective_value(optimizer::VrpOptimizer)

Get the objective function value after optimization.

"""
function get_objective_value(optimizer::VrpOptimizer)
    !optimizer.sol_defined && error("VRPSolver error: solution is undefined")
    user_form = optimizer.user_model.formulation
    user_vars = all_variables(user_form)
    obj_value = 0.0
    obj = objective_function(user_form)
    for user_var in user_vars
        obj_coeff = get(obj.terms, user_var, 0.0)
        var_val = get_value(optimizer, user_var)
        obj_value += obj_coeff * var_val
    end
    # CHECK: is this really needed?
    if optimizer.integer_objective
        return round(obj_value)
    else
        return obj_value
    end
end

"""
    get_value(optimizer::VrpOptimizer, user_var::JuMP.VariableRef)

Get the value for a decision variable after optimization.

"""
function get_value(optimizer::VrpOptimizer, user_var::JuMP.VariableRef)
    !optimizer.sol_defined && error("VRPSolver error: solution is undefined")
    optimizer_cols_info = optimizer.optimizer_cols_info
    if optimizer_cols_info.uservar_to_problem_type[user_var] == :DW_MASTER
        return get(optimizer.unmapped_vars_in_sol, user_var, 0.0)
    else
        val = 0.0
        for path_id in 1:length(optimizer.spsols_in_sol)
            val += get_value(optimizer, user_var, path_id)
        end
        return val
    end
end
function JuMP.value(user_var::JuMP.VariableRef)
    return get_value(user_var.model.ext[:optimizer], user_var)
end

"""
    get_values(optimizer::VrpOptimizer, user_vars::Array{JuMP.VariableRef,1})

Get the values for an array of decision variables after optimization.
   

"""
function get_values(optimizer::VrpOptimizer, user_vars::Array{JuMP.VariableRef,1})
    return [get_value(optimizer, user_var) for user_var in user_vars]
end

"""
    get_value(optimizer::VrpOptimizer, user_var::JuMP.VariableRef, path_id::Int)

Get the value for a decision variable due to a specific path.

`path_id` shoul be a value between 1 and the number of positive paths.
"""
function get_value(optimizer::VrpOptimizer, user_var::JuMP.VariableRef, path_id::Int)
    !optimizer.sol_defined && error("VRPSolver error: solution is undefined")
    !(1 <= path_id <= length(optimizer.spsols_in_sol)) &&
        error("VrpSolver error: invalid path id")
    return get(optimizer.spsols_in_sol[path_id].user_vars_in_sol, user_var, 0.0)
end

"""
    get_value(optimizer::VrpOptimizer, path_id::Int)

Get the value for a path variable (λ variable created internally due to the mapping).

"""
function get_value(optimizer::VrpOptimizer, path_id::Int)
    !optimizer.sol_defined && error("VRPSolver error: solution is undefined")
    !(1 <= path_id <= length(optimizer.spsols_in_sol)) &&
        error("VrpSolver error: invalid path id")
    return optimizer.spsols_in_sol[path_id][1]
end

"""
    get_values(optimizer::VrpOptimizer, user_vars::Array{JuMP.VariableRef,1}, path_id::Int)

Get the values for an array of decision variables that would be obtained from mapping only a single specific path variable (with value 1) to them. 
Necessary for identifying which paths are part of the solution in some models.

"""
function get_values(
    optimizer::VrpOptimizer, user_vars::Array{JuMP.VariableRef,1}, path_id::Int
)
    return [get_value(optimizer, user_var, path_id) for user_var in user_vars]
end

"""
    get_number_of_positive_paths(optimizer::VrpOptimizer)

Get the number of paths (lambda variables) with positive value in the solution. Those paths will be numbered from 1 to get_number_of_positive_paths for the purpose 
of retrieving the value of the lambda variables and for identifying the path.

"""
function get_number_of_positive_paths(optimizer::VrpOptimizer)
    !optimizer.sol_defined && error("VRPSolver error: solution is undefined")
    return length(optimizer.spsols_in_sol)
end

"""
    get_path_arcs(optimizer::VrpOptimizer, path_id::Int)

Returns a tuple where the first element is the `VrpGraph` and the second element is the sequence of arcs associated with the path.

"""
function get_path_arcs(optimizer::VrpOptimizer, path_id::Int)
    !optimizer.sol_defined && error("VRPSolver error: solution is undefined")
    !(1 <= path_id <= length(optimizer.spsols_in_sol)) &&
        error("VrpSolver error: invalid path id")

    spsol = optimizer.spsols_in_sol[path_id]
    graph = optimizer.user_model.graphs[spsol.graph_id]
    return (graph, spsol.arc_seq)
end

"""
    add_dynamic_constr!(optimizer::VrpOptimizer, vars, coeffs, sense, rhs, constr_name::String)

Add user cut (dynamic constraint) to the formulation.

It should be called inside a cut callback function registered with [`add_cut_callback!`](@ref).

# Arguments
- `vars`: array of variables of the cut
- `coeffs`: array of coefficients, where `coeffs[i]` is the coefficient of `vars[i]`  
- `sense`: sense of the cut: <=, >= or ==. 
- `rhs`: right-hand side value
- `constr_name`: a name for the cut

# Example

```julia
# let model be a VrpModel and x a set of variables
function my_callback()
   # add cut 1.0*x[1] + 2.0*x[2] <= 1.0
   add_dynamic_constr!(model.optimizer, [x[1],x[2]], [1.0,2.0], <=, 1.0, "my_cut")
end 
add_cut_callback!(model, my_callback, "my_callback")
```
"""
function add_dynamic_constr!(
    optimizer::VrpOptimizer, vars, coeffs, sense, rhs, constr_name::String
)
    constr_info = DynamicConstrInfo(vars, coeffs, sense, rhs)
    push!(optimizer.callbacks[constr_name].to_add_constrs, constr_info)
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

function get_enum_paths(user_model::VrpModel, paramfile::String)
    optimizer = VrpOptimizer(user_model, paramfile)
    paths = []
    bcsol = get_enumerated_sp_sols(optimizer.bapcod_model)
    status = c_next(bcsol)
    while status == 1
        graph_id = c_getProblemFirstId(bcsol) + 1
        bapcodarcids = c_getArcs(bcsol)
        graph = user_model.graphs[graph_id]
        path = Int[]
        for bid in bapcodarcids
            arc_id = graph.arc_bapcod_id_to_id[bid]
            push!(path, arc_id)
        end
        push!(paths, (graph, path))
        status = c_next(bcsol)
    end
    return paths
end

function print_enum_paths(model::VrpModel, paramfile::String)
    paths = get_enum_paths(model, paramfile)
    for (id, (graph, arcs)) in enumerate(paths)
        print("path $(id) graph $(graph.id): ")
        for arcid in arcs
            print(arcid, " ")
        end
        println()
    end
end

"""
    print_enum_paths(paths)
   
Print the enumerated paths for all graphs.

Warning 1: the enumeration procedure only produces ``\\mathcal{E}``-elementarity-paths. Moreover, for all paths that visit exactly the same elementarity sets, only a single least cost path is produced.

# Example
```julia
# Let `model` be a VrpModel and `path_to_params_file` the path for the parameters file
enum_paths, complete_form = get_complete_formulation(model, path_to_params_file)
complete_form.solver = CplexSolver() # set MIP solver (it can be another one than CPLEX)
print_enum_paths(enum_paths)
```
"""
function print_enum_paths(paths)
    println("\n\nEnumerated paths (v1->(arc_id)->v2->...):")
    net_vertex_id_map = Dict{Int,Dict{Int,Int}}()
    for (id, (graph, arcs)) in enumerate(paths)
        if !haskey(net_vertex_id_map, graph.id)
            net_vertex_id_map[graph.id] = Dict(
                value => key for (key, value) in graph.user_vertex_id_map
            ) # from user_id to net_id
        end
        print("path $(id) graph $(graph.id): ")
        arc = graph.arc_id_to_arc[arcs[1]]
        i = net_vertex_id_map[graph.id][arc.tail]
        j = if (graph.cycle_problem && arc.head == graph.sink_id)
            net_vertex_id_map[graph.id][graph.source_id]
        else
            net_vertex_id_map[graph.id][arc.head]
        end
        print("$i-($(arcs[1]))->$j")
        prev = j
        for arcid in arcs[2:end]
            arc = graph.arc_id_to_arc[arcid]
            j = if (graph.cycle_problem && arc.head == graph.sink_id)
                net_vertex_id_map[graph.id][graph.source_id]
            else
                net_vertex_id_map[graph.id][arc.head]
            end
            print("-($arcid)->$j")
            prev = j
        end
        println()
    end
    println()
end

"""
    get_complete_formulation(model::VrpModel, paramfile::String)

Get the complete formulation, which includes mapping constraints with λ variables (for paths).

The enumerated paths and the complete formulation are returned. It has all enumerated paths for all graphs and a JuMP.Model. 
This function can be seen as a tool for debugging the VRPSolver model, for example, when applied to small instances to check the correctness
of the final model produced. The user may solve it with another solver option in JuMP 0.18 (e.g. for Gurobi, Coin CBC or CPLEX).

Warning 1: the enumeration procedure only produces ``\\mathcal{E}``-elementarity-paths. Moreover, for all paths that visit exactly the same elementarity sets, only a single least cost path is produced.

Warning 2: if some user cuts are essential for the correctness of the model (i.e., they are not only used to improve the linear relaxation), the corresponding cut callbacks should be called for solving the complete formulation.
    
# Example
```julia
# Let `model` be a VrpModel and `path_to_params_file` the path for the parameters file
enum_paths, complete_form = get_complete_formulation(model, path_to_params_file)
complete_form.solver = CplexSolver() # set MIP solver (it can be another one than CPLEX)
print_enum_paths(enum_paths)
println(complete_form)
solve(complete_form)
println("Objective value: ", getobjectivevalue(complete_form))
```
"""
function get_complete_formulation(model::VrpModel, paramfile::String)
    paths = get_enum_paths(model, paramfile)

    optimizer_cols_info = extract_optimizer_cols_info(model)

    user_form = model.formulation
    formulation, orig_to_copied_uservar = copy_jump_model(user_form, optimizer_cols_info)

    # creating lamba variables (path variables)
    lambda_vars = []
    for i in 1:length(paths)
        v = @variable(formulation, base_name = "λ[$i]", lower_bound = 0, integer = true)
        push!(lambda_vars, v)
    end

    # mapping constraints
    paths_uservar_maps = []
    for (graph, path) in paths
        push!(paths_uservar_maps, get_path_uservar_map(path, graph))
    end
    for (orig_user_var, copied_user_var) in orig_to_copied_uservar
        if optimizer_cols_info.uservar_to_problem_type[orig_user_var] == :DW_SP
            coeffs = [1]
            vars = [copied_user_var]
            #lambda vars
            for (id, _) in enumerate(paths)
                if haskey(paths_uservar_maps[id], orig_user_var)
                    push!(coeffs, -paths_uservar_maps[id][orig_user_var])
                    push!(vars, lambda_vars[id])
                end
            end
            expr = JuMP.AffExpr()
            for (v, c) in zip(vars, coeffs)
                expr += c * v
            end
            @constraint(formulation, expr in MOI.EqualTo(0.0))
        end
    end

    # convexity constraints
    for graph in model.graphs
        expr = JuMP.AffExpr()
        for (id, (g, _)) in enumerate(paths)
            if g == graph
                expr += lambda_vars[id]
            end
        end
        @constraint(formulation, expr in MOI.GreaterThan(graph.multiplicity[1]))
        @constraint(formulation, expr in MOI.LessThan(graph.multiplicity[2]))
    end

    return paths, formulation
end

function copy_jump_model(original::JuMP.Model, optimizer_cols_info::OptimizerColsInfo)
    copied_model = Model()

    function substitute_affexpr(expr::JuMP.AffExpr, varmap)
        new_expr = JuMP.AffExpr(constant(expr))
        for (coef, v) in linear_terms(expr)
            new_expr += coef * varmap[v]
        end
        return new_expr
    end

    orig_to_copied_uservar = Dict{VariableRef,VariableRef}()
    for orig_uservar in all_variables(original)
        orig_to_copied_uservar[orig_uservar] = @variable(
            copied_model,
            base_name = name(orig_uservar),
            lower_bound = has_lower_bound(orig_uservar) ? lower_bound(orig_uservar) : -Inf,
            upper_bound = has_upper_bound(orig_uservar) ? upper_bound(orig_uservar) : Inf,
            integer = is_integer(orig_uservar)
        )
    end

    obj_expr = objective_function(original)
    copied_obj = substitute_affexpr(obj_expr, orig_to_copied_uservar)
    set_objective(copied_model, objective_sense(original), copied_obj)

    for set_type in (MOI.LessThan{Float64}, MOI.GreaterThan{Float64}, MOI.EqualTo{Float64})
        for c in all_constraints(original, AffExpr, set_type)
            cref = constraint_object(c)
            expr = cref.func
            sense = cref.set
            copied_expr = substitute_affexpr(expr, orig_to_copied_uservar)
            @constraint(copied_model, copied_expr in sense)
        end
    end

    return copied_model, orig_to_copied_uservar
end

function get_path_uservar_map(path::Vector{Int}, graph::VrpGraph)
    uservar_map = Dict{JuMP.VariableRef,Float64}()
    for arc_id in path
        arc = graph.arcs[arc_id]
        for (uservar, coef) in arc.vars
            if !haskey(uservar_map, uservar)
                uservar_map[uservar] = coef
            else
                uservar_map[uservar] += coef
            end
        end
    end
    return uservar_map
end

end
