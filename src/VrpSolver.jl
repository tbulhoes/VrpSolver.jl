module VrpSolver

using JuMP
using MathOptInterface: MathOptInterface
const MOI = MathOptInterface

import Base.show
using Printf

include("bapcod_interface.jl")
include("types.jl")
include("helpers.jl")
include("model.jl")
include("optimizer.jl")

export VrpModel, VrpGraph, VrpOptimizer
export add_resource!,
    set_resource_bounds!,
    add_elem_set_to_vertex_init_ng_neighbourhood!,
    add_arc!,
    add_arc_var_mapping!,
    set_arc_consumption!,
    set_arc_resource_bounds!,
    get_arc_consumption,
    set_arc_custom_res_params!,
    set_vertex_custom_res_params!,
    set_const_custom_res_params!,
    @register_custom_res_param_types,
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

end
