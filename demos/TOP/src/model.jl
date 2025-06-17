function build_model(data::DataTOP)

   A = arcs(data) # set of arcs of the input graph G
   n = num_customers(data)
   V1 = [i for i in 0:(n+1)]
   V⁺ = [i for i in 1:n]
   F = fleet_size(data)
   T = max_duration(data)

   # Formulation
   top = VrpModel()
   @variable(top.formulation, 0 <= y[i in V⁺] <= 1, Int)
   @variable(top.formulation, x[a in A], Int)
   @objective(top.formulation, Min, sum(-p(data, i) * y[i] for i in V⁺))
   @constraint(top.formulation, deg[i in V⁺], sum(x[a] for a in A if a[2] == i) == y[i])
   #println(top.formulation)

   function build_graph()

      v_source = 0
      v_sink = n + 1
      L, U = 0, F

      # node ids of G from 0 to n
      G = VrpGraph(top, V1, v_source, v_sink, (L, U))
      time_res_id = add_resource!(G, main=true)
      for i in V1
         l_i, u_i = 0.0, T # accumulated resource consumption interval [l_i, u_i] for the vertex i
         set_resource_bounds!(G, i, time_res_id, l_i, u_i)
      end

      # Build set of arcs 
      for (i, j) in A
         arc_id = add_arc!(G, i, j)
         add_arc_var_mapping!(G, arc_id, x[(i, j)])
         set_arc_consumption!(G, arc_id, time_res_id, t(data, (i, j)))
      end

      return G
   end

   G = build_graph()
   add_graph!(top, G)
   #println(G)

   set_vertex_packing_sets!(top, [[(G, i)] for i in V⁺])
   define_elementarity_sets_distance_matrix!(top, G, [[t(data, (i, j)) for j in V⁺] for i in V⁺])

   set_branching_priority!(top, "y", 1)
   set_branching_priority!(top, "x", 1)

   return (top, x, y)
end
