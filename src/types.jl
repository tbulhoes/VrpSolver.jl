@enum SetType NoSet = 0 ArcSet = 1 VertexSet = 2

mutable struct VrpVertex
    user_id::Int
    id::Int
    packing_sets::Vector{Int}
    elem_sets::Vector{Int}
    res_bounds::Dict{Int,Tuple{Float64,Float64}} #only for binary resources
    ng_set::Array{Int,1}
    custom_data::Dict{Int,Any}
end

mutable struct VrpArc
    id::Int
    tail::Int
    head::Int
    packing_sets::Vector{Int}
    elem_sets::Vector{Int}
    res_consumption::Array{Float64,1}
    res_bounds::Dict{Int,Tuple{Float64,Float64}} #only for non-binary resources
    vars::Array{Tuple{JuMP.VariableRef,Float64},1}
    ng_set::Array{Int,1}
    custom_data::Dict{Int,Any}
end

mutable struct VrpResource
    id::Int
    is_main::Bool
    is_binary::Bool
    is_disposable::Bool
    is_custom::Bool
    is_automatic::Bool
    step_size::Float64
    custom_data::Any
    cost_var::Union{JuMP.VariableRef,Nothing}
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
    multiplicity::Float64
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
    sol_from_model::Bool
    integer_objective::Bool
    initUB::Float64
    stats::Dict{} # execution statistics
    baptreedot_file::String
    optimizer_cols_info::OptimizerColsInfo
end
