
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
        "--RCSPrankOneCutsTypeToSeparate",
        "1",      # FIXME: Covering cuts makes VrpSolver crash
        "-t",
        baptreedot,
        "--ComputeDissagregateSpSol",
        "false",
    ],
)
    optimizer_cols_info = _extract_optimizer_cols_info(user_model)
    integer_objective = _has_integer_objective(user_model, optimizer_cols_info)
    use_meta_solver = _should_use_meta_solver(user_model)
    push!(fixed_params, "--RCSPuseMetaSolver")
    push!(fixed_params, use_meta_solver ? "1" : "0")

    bapcod_model_ptr = new!(
        param_file,
        true,
        integer_objective,
        false, # CHECK
        length(fixed_params),
        fixed_params,
    )

    #creating optimizer formulation
    _build_optimizer_vars_and_constrs(user_model, bapcod_model_ptr, optimizer_cols_info)

    #creating bapcod networks
    _generate_pricing_networks(user_model, bapcod_model_ptr, optimizer_cols_info)

    #branching priorities
    _set_branching_priorities_in_optimizer(
        user_model, bapcod_model_ptr, optimizer_cols_info
    )

    for rcc in user_model.cap_cuts_info
        wbcr_add_generic_capacity_cut(
            bapcod_model_ptr, rcc.capacity, rcc.demands, rcc.two_path_cuts_res_id
        )
    end

    for rcc in user_model.strongkpath_cuts_info
        wbcr_add_generic_strongkpath_cut(bapcod_model_ptr, rcc.capacity, rcc.demands)
    end

    if user_model.use_rank1_cuts
        wbc_add_generic_lim_mem_one_cut(bapcod_model_ptr)
    end

    if haskey(user_model.branching_priorities, "_res_cons_branching_")
        priority = user_model.branching_priorities["_res_cons_branching_"]
        c_add_elemset_resource_cons_branching(bapcod_model_ptr, Float64(priority))
    end
    if haskey(user_model.branching_priorities, "_ryanfoster_branching_")
        priority = user_model.branching_priorities["_ryanfoster_branching_"]
        c_add_packset_ryan_and_foster_branching(bapcod_model_ptr, Float64(priority))
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
        false,
        integer_objective,
        -1,
        Dict(),
        baptreedot,
        optimizer_cols_info,
    )
    user_model.optimizer = optimizer

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
    has_solution = _register_solutions(optimizer, sol_ptr)

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
            print("$(Int(round(getstatistic(optimizer.bapcod_model, :bcRecBestInc)))) & ")
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

function _extract_user_var_to_graphs(user_model::VrpModel)
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
        for res in user_model.graphs[graph_id].resources
            if !isnothing(res.cost_var)
                if !haskey(var_to_graphs, res.cost_var)
                    var_to_graphs[res.cost_var] = [graph_id]
                elseif !(graph_id in var_to_graphs[res.cost_var])
                    push!(var_to_graphs[res.cost_var], graph_id)
                end
            end
        end
    end
    return var_to_graphs
end

function _is_preprocessed_arc(graph::VrpGraph, arc::VrpArc)
    if !graph.cycle_problem && (arc.head == graph.source_id || arc.tail == graph.sink_id)
        return true
    end
    return false
end

function _generate_pricing_networks(
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
                wbcr_new_special_resource(
                    c_net_ptr, resource.id - 1, resource.is_disposable
                )
            elseif resource.is_custom
                colids = get(
                    optimizer_cols_info.uservar_to_colids, resource.cost_var, Int[]
                )
                colid = isempty(colids) ? -1 : colids[1]
                wbcr_new_custom_resource(c_net_ptr, resource.id - 1, bapcod_model, colid)
                wbcr_set_const_custom_res_params(
                    c_net_ptr, resource.id - 1, resource.custom_data
                )
            else
                wbcr_new_standard_resource(
                    c_net_ptr,
                    resource.id - 1,
                    resource.is_disposable,
                    resource.is_main,
                    resource.step_size,
                )
            end
            for vertex in graph.vertices
                if resource.is_binary
                    wbcr_set_vertex_special_res_params(
                        c_net_ptr,
                        vertex.id - 1,
                        resource.id - 1,
                        Int(vertex.res_bounds[resource.id][1]),
                        Int(vertex.res_bounds[resource.id][2]),
                    )
                elseif resource.is_custom
                    if haskey(vertex.custom_data, resource.id)
                        wbcr_set_vertex_custom_res_params(
                            c_net_ptr,
                            vertex.id - 1,
                            resource.id - 1,
                            vertex.custom_data[resource.id],
                        )
                    end
                else
                    if !haskey(vertex.res_bounds, resource.id)
                        wbcr_set_vertex_standard_res_params(
                            c_net_ptr, vertex.id - 1, resource.id - 1, -1e12, 1e12
                        )
                    else
                        wbcr_set_vertex_standard_res_params(
                            c_net_ptr,
                            vertex.id - 1,
                            resource.id - 1,
                            vertex.res_bounds[resource.id][1],
                            vertex.res_bounds[resource.id][2],
                        )
                    end
                end
            end
        end

        #adding vertices to packing and elem sets
        for vertex in graph.vertices
            for ps_id in vertex.packing_sets
                wbcr_add_vertex_to_packing_set(c_net_ptr, vertex.id - 1, ps_id - 1)
                wbcr_attach_elementarity_set_to_node(c_net_ptr, vertex.id - 1, ps_id - 1)
                if user_model.define_covering_sets
                    wbcr_add_vertex_to_covering_set(c_net_ptr, vertex.id - 1, ps_id - 1)
                end
            end
            for es_id in vertex.elem_sets
                wbcr_attach_elementarity_set_to_node(
                    c_net_ptr, vertex.id - 1, nbPackSets + es_id - 1
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
            if _is_preprocessed_arc(graph, arc)
                continue
            end
            arc_bapcod_id = wbcr_new_arc(c_net_ptr, arc.tail - 1, arc.head - 1, 0.0)
            graph.arc_bapcod_id_to_id[arc_bapcod_id] = arc.id
            for resource in graph.resources
                if resource.is_binary
                    wbcr_set_arc_special_res_params(
                        c_net_ptr,
                        arc_bapcod_id,
                        resource.id - 1,
                        Int(arc.res_consumption[resource.id]),
                    )
                elseif resource.is_custom
                    if haskey(arc.custom_data, resource.id)
                        wbcr_set_arc_custom_res_params(
                            c_net_ptr,
                            arc_bapcod_id,
                            resource.id - 1,
                            arc.custom_data[resource.id],
                        )
                    end
                else
                    wbcr_set_arc_standard_res_params(
                        c_net_ptr,
                        arc_bapcod_id,
                        resource.id - 1,
                        arc.res_bounds[resource.id][1],
                        arc.res_bounds[resource.id][2],
                        arc.res_consumption[resource.id],
                    )
                end
            end
            # adding arc to packing_set
            for ps_id in arc.packing_sets
                wbcr_add_edge_to_packing_set(c_net_ptr, arc_bapcod_id, ps_id - 1)
                wbcr_attach_elementarity_set_to_edge(c_net_ptr, arc_bapcod_id, ps_id - 1)
                if user_model.define_covering_sets
                    wbcr_add_edge_to_covering_set(c_net_ptr, arc_bapcod_id, ps_id - 1)
                end
            end
            for es_id in arc.elem_sets
                wbcr_attach_elementarity_set_to_edge(
                    c_net_ptr, arc_bapcod_id, nbPackSets + es_id - 1
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

function _extract_optimizer_cols_info(user_model::VrpModel)
    _check_resources_vars(user_model)
    user_var_to_graphs = _extract_user_var_to_graphs(user_model)
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
    obj_terms = get_obj_terms(user_form)
    nextcolid = 0
    user_vars = all_variables(user_form)
    for user_var in user_vars
        colsids = Int[]
        if !haskey(user_var_to_graphs, user_var) # unmapped variable
            push!(colsids, nextcolid)
            push!(cols_problems, (:DW_MASTER, 0))
            if !is_binary(user_var) &&
                (!has_lower_bound(user_var) || lower_bound(user_var) < 0.0)
                error(
                    "VRPSolver error: lower bound of unmapped variables must be nonnegative.",
                )
            end
            lb = is_binary(user_var) ? 0.0 : lower_bound(user_var)
            push!(cols_lbs, lb)
            ub = if is_binary(user_var)
                1.0
            elseif has_upper_bound(user_var)
                upper_bound(user_var)
            else
                Inf
            end
            push!(cols_ubs, ub)
            push!(cols_names, name(user_var))
            push!(cols_uservar, user_var)
            push!(cols_costs, get(obj_terms, user_var, 0.0))
            push!(
                cols_types,
                if is_integer(user_var) || is_binary(user_var)
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
                push!(cols_costs, get(obj_terms, user_var, 0.0))
                push!(
                    cols_types,
                    if is_integer(user_var) || is_binary(user_var)
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

function _build_optimizer_vars_and_constrs(
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
    nconstrs != length(constrs_refs) &&
        error("VRPSolver error: nonlinear constraints are not supported.")

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
        !iszero(JuMP.constraint_object(constr_ref).func.constant) && error(
            "VRPSolver error: AffExpr with nonzero constant is not allowed when defining constraints.",
        )
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

function _has_integer_objective(
    user_model::VrpModel, optimizer_cols_info::OptimizerColsInfo
)
    user_form = user_model.formulation
    user_vars = all_variables(user_form)
    obj_terms = get_obj_terms(user_form)
    for user_var in user_vars
        # checking if the coefficient in the objective function is integral
        var_cost = get(obj_terms, user_var, 0.0)
        if modf(var_cost)[1] != 0.0
            return false
        end
        # checking if unmapped variables are integer
        if optimizer_cols_info.uservar_to_problem_type[user_var] == :DW_MASTER &&
            !(is_integer(user_var) || is_binary(user_var))
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

function _set_branching_priorities_in_optimizer(
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
        for index in
            keys((exp_family isa Containers.SparseAxisArray) ? exp_family.data : exp_family)
            expr = exp_family[index]
            if expr isa VariableRef
                colsids, coeffs = Int[], Float64[]
                for colid in optimizer_cols_info.uservar_to_colids[expr]
                    push!(colsids, colid)
                    push!(coeffs, 1.0)
                end
            else
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
            end
            if isempty(colsids)
                @warn "VRPSolver warning: branching expression with name $(name) and index $(index) is empty"
            else
                add_branching_expression(
                    bapcod_model_ptr, exp_array_id, index, colsids, coeffs
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

# user_var must be a mapped or a resource variable
function _get_sp_var_value(
    optimizer::VrpOptimizer, bapcod_model, user_var::JuMP.VariableRef, bapcod_spsol
)
    optimizer_cols_info = optimizer.optimizer_cols_info
    subproblem_id = Int(c_getProblemFirstId(bapcod_spsol))

    var_value = 0.0
    for colid in optimizer_cols_info.uservar_to_colids[user_var]
        if optimizer_cols_info.cols_problems[colid + 1][2] == subproblem_id
            var_value += c_getValueOfVar(bapcod_model, bapcod_spsol, colid)
        end
    end
    return var_value
end

function _get_sp_sols(optimizer::VrpOptimizer, bapcodsol)
    bapcod_model = optimizer.bapcod_model
    user_vars = all_variables(optimizer.user_model.formulation)
    spsols_in_sol = SpSol[]
    optimizer_cols_info = optimizer.optimizer_cols_info

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
            Dict{JuMP.VariableRef,Float64}(),
            [graph.arcs[arc_id] for arc_id in arcs_ids],
        )

        for user_var in user_vars
            if optimizer_cols_info.uservar_to_problem_type[user_var] == :DW_MASTER
                continue
            end
            user_var_val = _get_sp_var_value(optimizer, bapcod_model, user_var, bapcodsol)
            if user_var_val > 0
                spsol.user_vars_in_sol[user_var] = mult * user_var_val
            end
        end

        push!(spsols_in_sol, spsol)
        status = c_next(bapcodsol)
    end

    return spsols_in_sol
end

function _register_solutions(optimizer::VrpOptimizer, bapcodsol; from_model = true)
    bapcod_model = optimizer.bapcod_model
    user_vars = all_variables(optimizer.user_model.formulation)
    optimizer_cols_info = optimizer.optimizer_cols_info
    unmapped_vars_in_sol = Dict{JuMP.VariableRef,Float64}()

    # emptying the solution stored in the optimizer
    empty!(optimizer.unmapped_vars_in_sol)
    empty!(optimizer.spsols_in_sol)

    if from_model
        # getting master solution if not already in `bapcodsol`
        if c_start(bapcodsol, bapcod_model) != 1
            optimizer.sol_defined = false
            return false
        end
    end

    # getting unmapped vars values
    for user_var in user_vars
        if optimizer_cols_info.uservar_to_problem_type[user_var] == :DW_MASTER
            colid = optimizer_cols_info.uservar_to_colids[user_var][1]
            unmapped_vars_in_sol[user_var] = c_getValueOfVar(bapcod_model, bapcodsol, colid)
        end
    end
    optimizer.unmapped_vars_in_sol = unmapped_vars_in_sol

    optimizer.spsols_in_sol = _get_sp_sols(optimizer, bapcodsol)
    optimizer.sol_from_model = from_model
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
    obj_terms = get_obj_terms(user_form)
    for user_var in user_vars
        obj_coeff = get(obj_terms, user_var, 0.0)
        var_val = get_value(optimizer, user_var)
        obj_value += obj_coeff * var_val
    end
    if optimizer.integer_objective && optimizer.sol_from_model
        return round(obj_value)
    else
        return obj_value
    end
end

"""
    get_value(optimizer::VrpOptimizer, user_var::JuMP.VariableRef)

Get the value for a decision variable after optimization or from a cut callback.

"""
function get_value(optimizer::VrpOptimizer, user_var::JuMP.VariableRef)
    !optimizer.sol_defined && error("VRPSolver error: solution is undefined")
    optimizer_cols_info = optimizer.optimizer_cols_info
    if optimizer_cols_info.uservar_to_problem_type[user_var] == :DW_MASTER
        return get(optimizer.unmapped_vars_in_sol, user_var, 0.0)
    else
        val = 0.0
        for path_id in 1:length(optimizer.spsols_in_sol)
            val += _get_value(optimizer, user_var, path_id)
        end
        return val
    end
end
function JuMP.value(user_var::JuMP.VariableRef)
    if (haskey(user_var.model.ext, :optimizer))
        return get_value(user_var.model.ext[:optimizer], user_var)
    else
        return MOI.get(owner_model(user_var), MOI.VariablePrimal(1), user_var)
    end
end

"""
    get_values(optimizer::VrpOptimizer, user_vars::Array{JuMP.VariableRef,1})

Get the values for an array of decision variables after optimization or from a cut callback.
   

"""
function get_values(optimizer::VrpOptimizer, user_vars::Array{JuMP.VariableRef,1})
    return [get_value(optimizer, user_var) for user_var in user_vars]
end

"""
    get_value(optimizer::VrpOptimizer, user_var::JuMP.VariableRef, path_id::Int)

Get the value for a decision variable due to a specific path of the solution. This function cannot be called
from a cut callback as VrpSolver only supports robust cuts.

`path_id` shoul be a value between 1 and the number of positive paths.
"""
function get_value(optimizer::VrpOptimizer, user_var::JuMP.VariableRef, path_id::Int)
    !optimizer.sol_defined && error("VRPSolver error: solution is undefined")
    !optimizer.sol_from_model && error(
        "VRPSolver error: only the values of original variables can be retrived in callbacks",
    )
    return _get_value(optimizer, user_var, path_id)
end

function _get_value(optimizer::VrpOptimizer, user_var::JuMP.VariableRef, path_id::Int)
    !(1 <= path_id <= length(optimizer.spsols_in_sol)) &&
        error("VrpSolver error: invalid path id")
    path = optimizer.spsols_in_sol[path_id]
    return get(path.user_vars_in_sol, user_var, 0.0)
end

"""
    get_value(optimizer::VrpOptimizer, path_id::Int)

Get the value for a path variable (λ variable created internally due to the mapping). This function cannot 
be called from a cut callback as VrpSolver only supports robust cuts.

"""
function get_value(optimizer::VrpOptimizer, path_id::Int)
    !optimizer.sol_defined && error("VRPSolver error: solution is undefined")
    !optimizer.sol_from_model && error(
        "VRPSolver error: only the values of original variables can be retrived in callbacks",
    )
    !(1 <= path_id <= length(optimizer.spsols_in_sol)) &&
        error("VrpSolver error: invalid path id")
    return optimizer.spsols_in_sol[path_id].multiplicity
end

"""
    get_values(optimizer::VrpOptimizer, user_vars::Array{JuMP.VariableRef,1}, path_id::Int)

Get the values for an array of decision variables due to a specific path of the solution. This function cannot be called
from a cut callback as VrpSolver only supports robust cuts.

"""
function get_values(
    optimizer::VrpOptimizer, user_vars::Array{JuMP.VariableRef,1}, path_id::Int
)
    return [get_value(optimizer, user_var, path_id) for user_var in user_vars]
end

"""
    get_number_of_positive_paths(optimizer::VrpOptimizer)

Get the number of paths (lambda variables) with positive value in the solution. Those paths will be numbered from 1 to get_number_of_positive_paths for the purpose 
of retrieving the value of the lambda variables and for identifying the path. This function cannot be called
from a cut callback as VrpSolver only supports robust cuts.

"""
function get_number_of_positive_paths(optimizer::VrpOptimizer)
    !optimizer.sol_defined && error("VRPSolver error: solution is undefined")
    !optimizer.sol_from_model && error(
        "VRPSolver error: only the values of original variables can be retrived in callbacks",
    )
    return length(optimizer.spsols_in_sol)
end

"""
    get_path_arcs(optimizer::VrpOptimizer, path_id::Int)

Returns a tuple where the first element is the `VrpGraph` and the second element is the sequence of arcs associated with the path.
This function cannot be called from a cut callback as VrpSolver only supports robust cuts.
"""
function get_path_arcs(optimizer::VrpOptimizer, path_id::Int)
    !optimizer.sol_defined && error("VRPSolver error: solution is undefined")
    !optimizer.sol_from_model && error(
        "VRPSolver error: only the values of original variables can be retrived in callbacks",
    )
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

function get_enum_paths(user_model::VrpModel, paramfile::String)
    optimizer = VrpOptimizer(
        user_model,
        paramfile;
        fixed_params = ["", "--RCSPmaxNumOfLabelsInHeurEnumeration", "1000"],
    )
    bcsol = get_enumerated_sp_sols(optimizer.bapcod_model)
    return _get_sp_sols(optimizer, bcsol)
end

"""
    print_enum_paths(model, paths)
   
Print the enumerated paths for all graphs.

Warning 1: the enumeration procedure only produces ``\\mathcal{E}``-elementarity-paths. Moreover, for all paths that visit exactly the same elementarity sets, only a single least cost path is produced.

# Example
```julia
# Let `model` be a VrpModel and `path_to_params_file` the path for the parameters file
enum_paths, complete_form = get_complete_formulation(model, path_to_params_file)
complete_form.solver = CplexSolver() # set MIP solver (it can be another one than CPLEX)
print_enum_paths(model, enum_paths)
```
"""
function print_enum_paths(user_model::VrpModel, paths::Vector{SpSol})
    println("\n\nEnumerated paths (v1->(arc_id)->v2->...):")
    net_vertex_id_map = Dict{Int,Dict{Int,Int}}()
    for (id, path) in enumerate(paths)
        graph = user_model.graphs[path.graph_id]
        arcs = path.arc_seq
        if !haskey(net_vertex_id_map, graph.id)
            net_vertex_id_map[graph.id] = Dict(
                value => key for (key, value) in graph.user_vertex_id_map
            ) # from user_id to net_id
        end
        print("path $(id) graph $(graph.id): ")
        arc = arcs[1]
        i = net_vertex_id_map[graph.id][arc.tail]
        j = if (graph.cycle_problem && arc.head == graph.sink_id)
            net_vertex_id_map[graph.id][graph.source_id]
        else
            net_vertex_id_map[graph.id][arc.head]
        end
        print("$i-($(arcs[1].id))->$j")
        prev = j
        for arc in arcs[2:end]
            j = if (graph.cycle_problem && arc.head == graph.sink_id)
                net_vertex_id_map[graph.id][graph.source_id]
            else
                net_vertex_id_map[graph.id][arc.head]
            end
            print("-($(arc.id))->$j")
            prev = j
        end
        println()
    end
    println()
end

"""
    get_complete_formulation(model::VrpModel, paramfile::String)

Get the complete formulation, which includes mapping constraints with λ variables (for paths). The λ variables are, by default,
integer, but integrality can be relaxed through the optional parameter `relax_lambda`.

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
print_enum_paths(model, enum_paths)
println(complete_form)
solve(complete_form)
println("Objective value: ", getobjectivevalue(complete_form))
```
"""
function get_complete_formulation(
    model::VrpModel, paramfile::String, relax_lambda::Bool = false
)
    paths = get_enum_paths(model, paramfile)

    optimizer_cols_info = _extract_optimizer_cols_info(model)

    user_form = model.formulation
    formulation, orig_to_copied_uservar = _copy_jump_model(user_form, optimizer_cols_info)

    # creating lamba variables (path variables)
    lambda_vars = []
    for i in 1:length(paths)
        v = @variable(
            formulation, base_name = "λ[$i]", lower_bound = 0, integer = !relax_lambda
        )
        push!(lambda_vars, v)
    end

    # mapping constraints
    for (orig_user_var, copied_user_var) in orig_to_copied_uservar
        optimizer_cols_info.uservar_to_problem_type[orig_user_var] == :DW_MASTER && continue
        coeffs = [1]
        vars = [copied_user_var]
        #lambda vars
        for (id, path) in enumerate(paths)
            if haskey(path.user_vars_in_sol, orig_user_var)
                push!(coeffs, -path.user_vars_in_sol[orig_user_var])
                push!(vars, lambda_vars[id])
            end
        end
        expr = JuMP.AffExpr()
        for (v, c) in zip(vars, coeffs)
            expr += c * v
        end
        @constraint(formulation, expr in MOI.EqualTo(0.0))
    end

    # convexity constraints
    for graph in model.graphs
        expr = JuMP.AffExpr()
        for (id, path) in enumerate(paths)
            if path.graph_id == graph.id
                expr += lambda_vars[id]
            end
        end
        @constraint(formulation, expr in MOI.GreaterThan(graph.multiplicity[1]))
        @constraint(formulation, expr in MOI.LessThan(graph.multiplicity[2]))
    end

    return paths, formulation
end

function _copy_jump_model(original::JuMP.Model, optimizer_cols_info::OptimizerColsInfo)
    copied_model = Model()

    function substitute_affexpr(expr::JuMP.AffExpr, varmap)
        new_expr = JuMP.AffExpr(constant(expr))
        for (coef, v) in linear_terms(expr)
            new_expr += coef * varmap[v]
        end
        return new_expr
    end

    function substitute_affexpr(var::JuMP.VariableRef, varmap)
        new_expr = JuMP.AffExpr(0.0)
        new_expr += varmap[var]
        return new_expr
    end

    orig_to_copied_uservar = Dict{VariableRef,VariableRef}()
    for orig_uservar in all_variables(original)
        orig_to_copied_uservar[orig_uservar] = @variable(
            copied_model,
            base_name = name(orig_uservar),
            lower_bound = has_lower_bound(orig_uservar) ? lower_bound(orig_uservar) : -Inf,
            upper_bound = has_upper_bound(orig_uservar) ? upper_bound(orig_uservar) : Inf,
            integer = is_integer(orig_uservar),
            binary = is_binary(orig_uservar)
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
