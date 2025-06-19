bapcod_path = get(ENV, "BAPCOD_RCSP_LIB", "")

## ccall redef
macro bcm_ccall(func, args...)
    f = "bcInterfaceModel_$(func)"
    args = map(esc, args)
    return quote
        ccall(($f, $bapcod_path), $(args...))
    end
end

macro bcr_ccall(func, args...)
    f = "bcRCSP_$(func)"
    args = map(esc, args)
    return quote
        ccall(($f, $bapcod_path), $(args...))
    end
end

macro bcs_ccall(func, args...)
    f = "bcInterfaceSolve_$(func)"
    args = map(esc, args)
    return quote
        ccall(($f, $bapcod_path), $(args...))
    end
end

macro bcsol_ccall(func, args...)
    f = "bcSolution_$(func)"
    args = map(esc, args)
    return quote
        ccall(($f, $bapcod_path), $(args...))
    end
end

function new!(
    param_file::String,
    print_param::Bool,
    int_obj::Bool,
    int_valued_bound::Bool,
    argc::Integer,
    argv::Array{String,1},
)
    @bcm_ccall("new", Ptr{Cvoid}, (Ptr{UInt8}, UInt8, UInt8, UInt8, Cint, Ptr{Ptr{UInt8}}),
        param_file, print_param, int_obj, int_valued_bound, Cint(argc), argv)
end

function init_model!(modelptr::Ptr{Cvoid}, nbrows::Integer, nbcols::Integer)
    @bcm_ccall("initModel", Cvoid, (Ptr{Cvoid}, Cint, Cint),
        modelptr, Cint(nbrows), Cint(nbcols))
end

function set_art_cost_value!(mptr::Ptr{Cvoid}, acv::Float64)
    @bcm_ccall("setArtCostValue", Cvoid, (Ptr{Cvoid}, Cdouble), mptr, Cdouble(acv))
end

function set_obj_ub!(mptr::Ptr{Cvoid}, ub::Float64)
    @bcm_ccall("setObjUb", Cvoid, (Ptr{Cvoid}, Cdouble), mptr, Cdouble(ub))
end

function register_sub_problem!(mptr::Ptr{Cvoid}, subproblemtype::Cint, spmid::Vector{Cint})
    @bcm_ccall("registerSubProblem", Cint, (Ptr{Cvoid}, Cint, Ptr{Cint}),
        mptr, subproblemtype, spmid)
end

function sptype_to_int(sp_type::Symbol)
    if sp_type == :MIP
        return 0
    elseif sp_type == :DW_MASTER || sp_type == :B_MASTER
        return 1
    elseif sp_type == :DW_SP
        return 2
    elseif sp_type == :B_SP
        return 3
    elseif sp_type == :ALL
        return 4
    else
        error(
            "Cannot recognize problem type : $sp_type. It must be :DW_MASTER or :DW_SP for Dantzig-Wolfe decomposition and :B_MASTER or :B_SP for Benders decomposition.",
        )
    end
end

function toArray(a)
    if isa(a, Integer)
        return [a]
    end
    arr = Vector{Int}()
    for i in a
        arr = vcat(arr, toArray(i))
    end
    return arr
end

# CHECK: acho que nao precisamos disso
function createMultiIndex(array_mid::Vector{Cint}, array::Vector{Int})
    length_array = length(array)
    length_array_mid = length(array_mid)
    (length(array) > 8) && error("BaPCod does not support multi-index with more than 8 indices.")
    for i in 1:length_array_mid
        if i <= length_array
            array_mid[i] = array[i]
        else
            array_mid[i] = (i == 1) ? 0 : -1
        end
    end
end

from_index_to_BaPCodindex(id, array_mid::Vector{Cint}) = createMultiIndex(array_mid, toArray(id))

function c_register_subproblems(mptr::Ptr{Cvoid}, spids)
    for (spid, sptype) in spids
        spmid = Array{Cint,1}(undef, 8)
        from_index_to_BaPCodindex(spid, spmid)
        subproblemtype = Cint(sptype_to_int(sptype))
        register_sub_problem!(mptr, subproblemtype, spmid)
    end
end

function register_var!(
    mptr::Ptr{Cvoid},
    name::Symbol,
    column_id::Integer,
    sp_type::Integer,
    sp_mid::Vector{Cint},
    var_mid::Vector{Cint},
)
    @bcm_ccall("registerVar", Cint, (Ptr{Cvoid}, Ptr{Cchar}, Cint, Cint, Ptr{Cint}, Ptr{Cint}),
        mptr, name, Cint(column_id), sp_type, sp_mid, var_mid)
end

function init_vars!(mptr::Ptr{Cvoid}, l::Vector{Cdouble}, u::Vector{Cdouble}, c::Vector{Cdouble}, t::Vector{Cchar})
    @bcm_ccall("initVars", Cvoid, (Ptr{Cvoid}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
        mptr, l, u, c)
    @bcm_ccall("setVarType", Cint, (Ptr{Cvoid}, Ptr{Cchar}, Cint),
        mptr, t, Cint(length(t)))
end

function c_register_vars(mptr::Ptr{Cvoid}, lbs::Vector{Float64}, ubs::Vector{Float64}, costs::Vector{Float64}, cols_types::Vector{Char}, vars_decomposition)
    var_bcid = Array{Cint,1}(undef, 8)
    sp_bcid = Array{Cint,1}(undef, 8)

    for (column_id, (name, v_id, sp_type, sp_id)) in enumerate(vars_decomposition)
        # BaPCod needs an index
        from_index_to_BaPCodindex(v_id, var_bcid)
        from_index_to_BaPCodindex(sp_id, sp_bcid)
        # Register the variable
        register_var!(mptr, name, column_id - 1, sptype_to_int(sp_type), sp_bcid, var_bcid)
    end
    init_vars!(mptr, [Cdouble(lb) for lb in lbs], [Cdouble(ub) for ub in ubs], [Cdouble(cost) for cost in costs], [Cchar(t) for t in cols_types])
end

function init_cstrs!(
    mptr::Ptr{Cvoid},
    starts::Vector{Cint},
    rows_id::Vector{Cint},
    nonzeros::Vector{Cdouble},
    lb::Vector{Cdouble},
    ub::Vector{Cdouble},
)
    @bcm_ccall("initCstrs", Cint, (Ptr{Cvoid}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
        mptr, starts, rows_id, nonzeros, lb, ub)
end

function register_cstr!(
    mptr::Ptr{Cvoid},
    name::Symbol,
    row_id::Cint,
    sp_type::Integer,
    sp_mid::Vector{Cint},
    cstr_mid::Vector{Cint},
)
    @bcm_ccall("registerCstr", Cint, (Ptr{Cvoid}, Ptr{Cchar}, Cint, Cint, Ptr{Cint}, Ptr{Cint}),
        mptr, string(name), row_id, sp_type, sp_mid, cstr_mid)
end

struct CMatrix
    starts::Array{Cint,1}
    rows_id::Array{Cint,1}
    nonzeros::Array{Cdouble,1}
end

function c_register_cstrs(mptr::Ptr{Cvoid}, starts::Array{Int,1}, rows_id::Array{Int,1}, nonzeros::Array{Float64,1}, lbs::Array{Float64,1}, ubs::Array{Float64,1}, cstrs_decomposition)
    cstr_bcid = Array{Cint,1}(undef, 8)
    sp_bcid = Array{Cint,1}(undef, 8)

    for (row_id, (name, c_id, sp_type, sp_id)) in enumerate(cstrs_decomposition)
        # BaPCod needs an idnex
        from_index_to_BaPCodindex(c_id, cstr_bcid)
        from_index_to_BaPCodindex(sp_id, sp_bcid)
        # Register the constraint
        register_cstr!(mptr, name, Cint(row_id - 1), sptype_to_int(sp_type), sp_bcid, cstr_bcid)
    end
    init_cstrs!(mptr, [Cint(s) for s in starts], [Cint(r) for r in rows_id], [Cdouble(nz) for nz in nonzeros], [Cdouble(lb) for lb in lbs], [Cdouble(ub) for ub in ubs])
end

function sub_problem_mult!(mptr::Ptr{Cvoid}, mult_lb::Cint, mult_ub::Cint, sp_type::Integer, sp_bcid::Vector{Cint})
    @bcm_ccall("subProblemMult", Cint, (Ptr{Cvoid}, Cint, Cint, Cint, Ptr{Cint}),
        mptr, mult_lb, mult_ub, sp_type, sp_bcid)
end

function c_set_sp_multiplicities(mptr::Ptr{Cvoid}, sp_mult)
    sp_bcid = Array{Cint,1}(undef, 8)
    for (sp_id, sp_type, mult_lb, mult_ub) in sp_mult
        from_index_to_BaPCodindex(sp_id, sp_bcid)
        status = sub_problem_mult!(mptr, Cint(mult_lb), Cint(mult_ub), sptype_to_int(sp_type), sp_bcid)
        (status != 1) &&
            error("Cannot set multiplicity on the subproblem with the index $sp_id. Make sure it exists.")
    end
end

function c_add_packset_ryan_and_foster_branching(
    modelptr::Ptr{Cvoid},
    priority::Float64,
)
    @bcr_ccall("addPackSetRyanAndFosterBranching", Cint, (Ptr{Cvoid}, Cdouble),
        modelptr, Cdouble(priority))
end

function c_add_elemset_resource_cons_branching(
    modelptr::Ptr{Cvoid},
    priority::Float64,
)
    @bcr_ccall("addElemSetResourceConsumptionBranching", Cint, (Ptr{Cvoid}, Cdouble),
        modelptr, Cdouble(priority))
end

function set_var_priority_in_master!(
    modelptr::Ptr{Cvoid},
    varname::Symbol,
    sp_bctype::Integer,
    sp_bcid::Vector{Cint},
    priority::Cdouble,
)
    @bcm_ccall("setVarPriorityInMaster", Cint, (Ptr{Cvoid}, Ptr{Cchar}, Cint, Ptr{Cint}, Cdouble),
        modelptr, varname, sp_bctype, sp_bcid, priority)
end

function c_vars_branching_priorities(modelptr::Ptr{Cvoid}, p::Vector{Tuple{Symbol,Symbol,Int,Float64}})
    sp_bcid = Array{Cint,1}(undef, 8)
    for (varname, sp_type, sp_id, priority) in p
        from_index_to_BaPCodindex(sp_id, sp_bcid)
        sp_bctype = sptype_to_int(sp_type)
        status = set_var_priority_in_master!(modelptr, varname, sp_bctype, sp_bcid, Cdouble(priority))
        (status != 1) && error("Cannot set branching priority on variables named $varname.")
    end
end

function wbcr_new(
    c_model::Ptr{Cvoid},
    sp_bctype::Integer,
    sp_bcid::Vector{Cint},
    nb_nodes::Integer,
    nb_es::Integer,
    nb_ps::Integer,
    nb_cs::Integer,
)
    ptr = @bcr_ccall("new", Ptr{Cvoid}, (Ptr{Cvoid}, Cint, Ptr{Int}, Cint, Cint, Cint, Cint),
        c_model, Cint(sp_bctype), sp_bcid, Cint(nb_nodes), Cint(nb_es), Cint(nb_ps), Cint(nb_cs))
end

function new_network!(
    c_model::Ptr{Cvoid},
    sp_id::Int,
    sp_type::Symbol,
    nb_nodes::Int,
    nb_es::Int,
    nb_ps::Int,
    nb_cs::Int,
)
    sp_bcid = Array{Cint,1}(undef, 8)
    from_index_to_BaPCodindex(sp_id, sp_bcid)
    sp_bctype = sptype_to_int(sp_type)
    return wbcr_new(c_model, sp_bctype, sp_bcid, nb_nodes, nb_es, nb_ps, nb_cs)
end

function wbcr_new_resource(c_net::Ptr{Cvoid}, res_id::Integer)
    status = @bcr_ccall("newResource", Cint, (Ptr{Cvoid}, Cint),
        c_net, Cint(res_id))
    (status != 1) && error("Cannot create resource $res_id.")
end

function wbcr_set_as_main_resource(c_net::Ptr{Cvoid}, res_id::Integer, stepvalue::Float64)
    @bcr_ccall("setAsMainResource", Cvoid, (Ptr{Cvoid}, Cint, Cdouble), c_net, Cint(res_id), Cdouble(stepvalue))
end

function wbcr_set_vertex_consumption_lb(c_net::Ptr{Cvoid}, n_id::Integer, res_id::Integer, lb::Float64)
    @bcr_ccall("setVertexConsumptionLB", Cint, (Ptr{Cvoid}, Cint, Cint, Cdouble),
        c_net, Cint(n_id), Cint(res_id), Cdouble(lb))
end

function wbcr_set_vertex_consumption_ub(c_net::Ptr{Cvoid}, n_id::Integer, res_id::Integer, ub::Float64)
    @bcr_ccall("setVertexConsumptionUB", Cint, (Ptr{Cvoid}, Cint, Cint, Cdouble),
        c_net, Cint(n_id), Cint(res_id), Cdouble(ub))
end

function wbcr_set_arc_consumption_lb(c_net::Ptr{Cvoid}, arc_id::Integer, res_id::Integer, lb::Float64)
    @bcr_ccall("setArcConsumptionLB", Cint, (Ptr{Cvoid}, Cint, Cint, Cdouble),
        c_net, Cint(arc_id), Cint(res_id), Cdouble(lb))
end

function wbcr_set_arc_consumption_ub(c_net::Ptr{Cvoid}, arc_id::Integer, res_id::Integer, ub::Float64)
    @bcr_ccall("setArcConsumptionUB", Cint, (Ptr{Cvoid}, Cint, Cint, Cdouble),
        c_net, Cint(arc_id), Cint(res_id), Cdouble(ub))
end

function wbcr_attach_elementarity_set_to_node(c_net::Ptr{Cvoid}, n_id::Integer, es_id::Integer)
    @bcr_ccall("attachElementaritySetToNode", Cint, (Ptr{Cvoid}, Cint, Cint),
        c_net, Cint(n_id), Cint(es_id))
end

function wbcr_attach_elementarity_set_to_edge(c_net::Ptr{Cvoid}, edge_id::Integer, es_id::Integer)
    @bcr_ccall("attachElementaritySetToEdge", Cint, (Ptr{Cvoid}, Cint, Cint),
        c_net, Cint(edge_id), Cint(es_id))
end

function wbcr_add_arc_to_mem_of_elementarity_set(c_net::Ptr{Cvoid}, edge_id::Integer, es_id::Integer)
    status = @bcr_ccall("addEdgeToMemOfElementaritySet", Cint, (Ptr{Cvoid}, Cint, Cint),
        c_net, Cint(edge_id), Cint(es_id))
    (status != 1) && error("Cannot add arc $edge_id to memory of elementarity set $es_id.")
end

function wbcr_add_vertex_to_mem_of_elementarity_set(c_net::Ptr{Cvoid}, n_id::Integer, es_id::Integer)
    status = @bcr_ccall("addVertexToMemOfElementaritySet", Cint, (Ptr{Cvoid}, Cint, Cint),
        c_net, Cint(n_id), Cint(es_id))
    (status != 1) && error("Cannot add vertex $n_id to memory of elementarity set $es_id.")
end

function wbcr_add_vertex_to_packing_set(c_net::Ptr{Cvoid}, n_id::Integer, ps_id::Integer)
    status = @bcr_ccall("addVertexToPackingSet", Cint, (Ptr{Cvoid}, Cint, Cint),
        c_net, Cint(n_id), Cint(ps_id))
end

function wbcr_add_edge_to_packing_set(c_net::Ptr{Cvoid}, edge_id::Integer, ps_id::Integer)
    status = @bcr_ccall("addEdgeToPackingSet", Cint, (Ptr{Cvoid}, Cint, Cint),
        c_net, Cint(edge_id), Cint(ps_id))
end

function wbcr_add_vertex_to_covering_set(c_net::Ptr{Cvoid}, n_id::Integer, cs_id::Integer)
    status = @bcr_ccall("addVertexToCoveringSet", Cint, (Ptr{Cvoid}, Cint, Cint),
        c_net, Cint(n_id), Cint(cs_id))
end

function wbc_add_generic_lim_mem_one_cut(c_model::Ptr{Cvoid})
    status = @bcr_ccall("addGenericLimMemOneCut", Cint, (Ptr{Cvoid},),
        c_model)
    (status != 1) && error("Cannot add the generic lim-mem-one-rank cut.")
end

function wbcr_add_generic_capacity_cut(
    c_model::Ptr{Cvoid},
    max_cap::Int,
    dem::Vector{Int},
)
    cint_dem = [Cint(d) for d in dem]
    status = @bcr_ccall("addGenericCapacityCut", Cint,
        (Ptr{Cvoid}, Cint, Ptr{Cint}, Cint, Cint, Cdouble, Cdouble, Cint),
        c_model, Cint(max_cap), cint_dem, Cint(length(dem)), Cint(0), 3.0, 1.0, Cint(-1))
    (status != 1) && error("Cannot add the generic capacity cut generator.")
end

function wbcr_set_source(c_net::Ptr{Cvoid}, n_id::Integer)
    @bcr_ccall("setSource", Cint, (Ptr{Cvoid}, Cint), c_net, Cint(n_id))
end

function wbcr_set_sink(c_net::Ptr{Cvoid}, n_id::Integer)
    @bcr_ccall("setSink", Cint, (Ptr{Cvoid}, Cint), c_net, Cint(n_id))
end

function wbcr_new_arc(c_net::Ptr{Cvoid}, src::Integer, dst::Integer, cost::Float64)
    edge_id = @bcr_ccall("newArc", Cint, (Ptr{Cvoid}, Cint, Cint, Cdouble),
        c_net, Cint(src), Cint(dst), Cdouble(cost))
    return edge_id
end

function wbcr_attach_bcvar_to_arc(c_net::Ptr{Cvoid}, edge_id::Integer, c_model::Ptr{Cvoid}, var_col::Int, coeff::Float64)
    @bcr_ccall("attachBcVarToArc", Cint, (Ptr{Cvoid}, Cint, Ptr{Cvoid}, Cint, Cdouble),
        c_net, Cint(edge_id), c_model, Cint(var_col), Cdouble(coeff))
end

function wbcr_set_edge_consumption_value(c_net::Ptr{Cvoid}, edge_id::Integer, res_id::Integer, value::Float64)
    @bcr_ccall("setEdgeConsumptionValue", Cint, (Ptr{Cvoid}, Cint, Cint, Cdouble),
        c_net, Cint(edge_id), Cint(res_id), Cdouble(value))
end

function wbcr_set_as_nondisposable_resource(c_net::Ptr{Cvoid}, res_id::Integer)
    @bcr_ccall("setAsNonDisposableResource", Cvoid, (Ptr{Cvoid}, Cint),
        c_net, Cint(res_id))
end

function wbcr_set_special_as_nondisposable_resource(c_net::Ptr{Cvoid}, res_id::Integer)
    @bcr_ccall("setSpecialResourceAsNonDisposable", Cvoid, (Ptr{Cvoid}, Cint),
        c_net, Cint(res_id))
end

function wbcr_set_vertex_special_consumption_lb(c_net::Ptr{Cvoid}, n_id::Integer, res_id::Integer, lb::Float64)
    @bcr_ccall("setVertexSpecialConsumptionLB", Cint, (Ptr{Cvoid}, Cint, Cint, Cdouble),
        c_net, Cint(n_id), Cint(res_id), Cdouble(lb))
end

function wbcr_set_vertex_special_consumption_ub(c_net::Ptr{Cvoid}, n_id::Integer, res_id::Integer, ub::Float64)
    @bcr_ccall("setVertexSpecialConsumptionUB", Cint, (Ptr{Cvoid}, Cint, Cint, Cdouble),
        c_net, Cint(n_id), Cint(res_id), Cdouble(ub))
end

function wbcr_set_edge_special_consumption_value(
    c_net::Ptr{Cvoid},
    edge_id::Integer,
    res_id::Integer,
    value::Float64,
)
    @bcr_ccall("setEdgeSpecialConsumptionValue", Cint, (Ptr{Cvoid}, Cint, Cint, Cdouble),
        c_net, Cint(edge_id), Cint(res_id), Cdouble(value))
end

function wbcr_set_elementarity_sets_distance_matrix(c_net::Ptr{Cvoid}, matrix::Array{Array{Float64,1},1}, nb_es::Integer)
    if length(matrix) != nb_es
        error("distance matrix size does not equal to the number of elementarity sets")
        for row in matrix
            if length(row) != nb_es
                error("distance matrix size does not equal to the number of elementarity sets")
            end
        end
    end
    cdouble_matrix = [[Cdouble(i) for i in row] for row in matrix]
    @bcr_ccall("setElemSetsDistanceMatrix", Cint, (Ptr{Cvoid}, Ptr{Ptr{Cdouble}}, Cint),
        c_net, cdouble_matrix, Cint(nb_es))
end

function wbcr_create_oracle(
    c_net::Ptr{Cvoid},
    c_model::Ptr{Cvoid},
    sp_type::Integer,
    sp_id::Vector{Cint},
    save_standalone::Bool,
    standalone_filename::String,
)
    @bcr_ccall("createOracle", Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Cint, Ptr{Cint}, UInt8, Ptr{UInt8}),
        c_net, c_model, sp_type, sp_id, save_standalone, standalone_filename)
end

function new_oracle!(c_net::Ptr{Cvoid}, c_model::Ptr{Cvoid}, sp_type::Symbol, sp_id::Int)
    sp_bcid = Array{Cint,1}(undef, 8)
    from_index_to_BaPCodindex(sp_id, sp_bcid)
    sp_bctype = sptype_to_int(sp_type)
    wbcr_create_oracle(c_net, c_model, sp_bctype, sp_bcid, false, "")
end

function c_optimize(modelptr::Ptr{Cvoid}, solution::Ptr{Cvoid})
    @bcs_ccall("optimize", Cvoid, (Ptr{Cvoid}, Ptr{Cvoid}), modelptr, solution)
end

function new_sol!()
    @bcsol_ccall("new", Ptr{Cvoid}, ())
end

function c_getValueOfVar(mptr::Ptr{Cvoid}, solution::Ptr{Cvoid}, colid::Int)
    value = Ref{Cdouble}(0.0)
    @bcsol_ccall("getValueOfVar", Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Cint, Ref{Cdouble}),
        mptr, solution, Cint(colid), value)
    return Float64(value[])
end

function c_getMultiplicity(solution::Ptr{Cvoid})
    mult = Ref{Cint}(0)
    @bcsol_ccall("getMultiplicity", Cint, (Ptr{Cvoid}, Ref{Cint}),
        solution, mult)
    return Int(mult[])
end

function c_getProblemFirstId(solution::Ptr{Cvoid})
    fid = @bcsol_ccall("getProblemFirstId", Cint, (Ptr{Cvoid},), solution)
    return Int(fid)
end

function c_start(solution::Ptr{Cvoid}, mptr)
    @bcsol_ccall("start", Cint, (Ptr{Cvoid}, Ptr{Cvoid}), solution, mptr)
end

function c_next(solution::Ptr{Cvoid})
    @bcsol_ccall("next", Cint, (Ptr{Cvoid},), solution)
end

function c_getArcs(solution::Ptr{Cvoid})
    nb_arcs = @bcsol_ccall("getNbNodes", Cint, (Ptr{Cvoid},), solution) - 1
    arcs_ids = [Cint(0) for _ in 1:nb_arcs]
    @bcsol_ccall("getArcsIds", Cint, (Ptr{Cvoid}, Ptr{Cint}, Cint), solution, arcs_ids, nb_arcs)
    return [Int(aid) for aid in arcs_ids]
end
