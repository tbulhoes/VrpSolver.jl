function build_model(data::DataCVRP, app::Dict{String,Any})

   E = edges(data) # set of edges of the input graph G′
   n = nb_customers(data)
   V = [i for i in 0:n] # set of vertices of the graphs G′ and G
   V⁺ = [i for i in 1:n] # set of customers of the input graph G′
   Q = veh_capacity(data)

   # Formulation
   cvrp = VrpModel()
   @variable(cvrp.formulation, x[e in E], Int)
   @objective(cvrp.formulation, Min, sum(c(data,e) * x[e] for e in E))
   @constraint(cvrp.formulation, deg[i in V⁺], sum(x[e] for e in δ(data, i)) == 2.0)

   #println(cvrp.formulation)

   # Build the model directed graph G=(V,A)
   function build_graph()

      v_source = v_sink = 0

      L = max(app["minr"], lowerBoundNbVehicles(data))
      U = min(app["maxr"], n)

      # node ids of G from 0 to n
      G = VrpGraph(cvrp, V, v_source, v_sink, (L, U))
      cap_res_id = add_resource!(G, main = true) # R = R_M = {cap_res_id}

      for i in V
         l_i, u_i = 0.0, Float64(Q) # accumulated resource consumption interval [l_i, u_i] for the vertex i
         set_resource_bounds!(G, i, cap_res_id, l_i, u_i)
      end

      # Build set of arcs A from E (two arcs for each edge (i,j))
      for (i,j) in E
         arc_id = add_arc!(G, i, j)
         add_arc_var_mapping!(G, arc_id, x[(i,j)])
         #option for adding an arc and already mapping it to a set of vars
         #arc_id = add_arc!(G, i, j, [x[(i,j)]])
         set_arc_consumption!(G, arc_id, cap_res_id, (d(data, i) + d(data, j))/2)
         arc_id = add_arc!(G, j, i)
         add_arc_var_mapping!(G, arc_id, x[(i,j)])
         set_arc_consumption!(G, arc_id, cap_res_id, (d(data, i) + d(data, j))/2)
      end

      return G
   end

   G = build_graph()
   add_graph!(cvrp, G)
   #println(G)

   set_vertex_packing_sets!(cvrp, [[(G,i)] for i in V⁺])

   define_elementarity_sets_distance_matrix!(cvrp, G, [[distance(data, (i, j)) for j in V⁺] for i in V⁺])
 
   add_capacity_cut_separator!(cvrp, [ ( [(G,i)], d(data, i) ) for i in V⁺], Q)
 
   set_branching_priority!(cvrp, "x", 1)

   function edge_ub_callback()
      for (i,j) in E
        e = (i,j)
         if i != 0 && get_value(cvrp.optimizer, x[e]) > 1.001
            println("Adding edge ub cut for e = ", e)
            add_dynamic_constr!(cvrp.optimizer, [x[e]], [1.0], <=, 1.0, "edge_ub")
         end
      end
   end
   add_cut_callback!(cvrp, edge_ub_callback, "edge_ub")

   return (cvrp, x)
end
