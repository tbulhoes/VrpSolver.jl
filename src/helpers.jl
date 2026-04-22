_contains(p, s) = !isnothing(findnext(s, p, 1))

function _check_id(id::Int, min_id::Int, max_id::Int)
    (id < min_id || id > max_id) && error("Unknown id $id.")
end

function _check_vertex_id(graph::VrpGraph, id::Int)
    !(id in keys(graph.user_vertex_id_map)) && error("Unknown vertex id $id.")
end

function split_var_name(var)
    var_name = name(var)
    var_container_name = _contains(var_name, "[") ? split(var_name, "[")[1] : var_name
    var_id = _contains(var_name, "[") ? split(split(var_name, "[")[2], "]")[1] : "1"
    var_id = replace(var_id, r"\(" => s"")
    var_id = replace(var_id, r"\)" => s"")
    var_id = tuple([parse(Int64, x) for x in split(var_id, ",") if x != ""]...)
    return (var_container_name, var_id)
end

function get_obj_terms(jump_model::Model)
    obj = objective_function(jump_model)
    obj_terms = if isa(obj, JuMP.GenericAffExpr)
        !iszero(obj.constant) && error(
            "VRPSolver error: constant part of the objective function must be zero"
        )
        obj.terms
    elseif isa(obj, JuMP.VariableRef)
        Dict(obj => 1.0)
    else
        error(
            "VRPSolver error: supported types for the objective function are: GenericAffExpr and VariableRef",
        )
    end
    obj_terms
end
