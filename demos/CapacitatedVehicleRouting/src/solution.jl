mutable struct Solution
   cost::Union{Int,Float64}
   routes::Array{Array{Int}}
end

# build Solution from the variables x
function getsolution(data::DataCVRP, optimizer::VrpOptimizer, x, objval, app::Dict{String,Any})
   E, dim = edges(data), dimension(data)
   adj_list = [[] for i in 1:dim] 
   for e in E
      val = get_value(optimizer, x[e])
      if val > 0.5
         push!(adj_list[e[1]+1], e[2])
         push!(adj_list[e[2]+1], e[1])
         if val > 1.5
            push!(adj_list[e[1]+1], e[2])
            push!(adj_list[e[2]+1], e[1])
         end
      end
   end
   visited, routes = [false for i in 2:dim], []
   for i in adj_list[1]
      if !visited[i]
         r, prev = [], 0
         push!(r, i)
         visited[i] = true
         length(adj_list[i+1]) != 2 && error("Problem trying to recover the route from the x values. "*
                                             "Customer $i has $(length(adj_list[i+1])) incident edges.")
         next, prev = (adj_list[i+1][1] == prev) ? adj_list[i+1][2] : adj_list[i+1][1], i
         maxit, it = dim, 0
         while next != 0 && it < maxit
            length(adj_list[next+1]) != 2 && error("Problem trying to recover the route from the x values. "* 
                                                   "Customer $next has $(length(adj_list[next+1])) incident edges.")
            push!(r, next)
            visited[next] = true
            aux = next
            next, prev = (adj_list[next+1][1] == prev) ? adj_list[next+1][2] : adj_list[next+1][1], aux
            it += 1
         end
         (it == maxit) && error("Problem trying to recover the route from the x values. "*
                                "Some route can not be recovered because the return to depot is never reached")
         push!(routes, r)
      end
   end 
   !isempty(filter(a->a==false,visited)) && error("Problem trying to recover the route from the x values. "*
                              "At least one customer was not visited or there are subtours in the solution x.")
   if !app["noround"]
      objval = trunc(Int, round(objval))
   end
   return Solution(objval, routes)
end

function print_routes(solution)
   for (i,r) in enumerate(solution.routes)
      print("Route #$i: ") 
      for j in r
         print("$j ")
      end
      println()
   end
end

# checks the feasiblity of a solution
function checksolution(data::DataCVRP, solution)
   dim, Q = dimension(data), veh_capacity(data)
   visits = [0 for i in 2:dim]
   sum_cost = 0.0
   for (i,r) in enumerate(solution.routes)
      sum_demand, prev = 0.0, 0
      for j in r
         visits[j] += 1
         (visits[j] == 2) && error("Customer $j was visited more than once")
         sum_cost += distance(data, (prev,j))
         sum_demand += d(data, j)
         prev = j
      end
      sum_cost += distance(data, (prev,0))
      (sum_demand > Q) && error("Route #$i is violating the capacity constraint. Sum of the demands is $(sum_demand) and Q is $Q")
   end
   !isempty(filter(a->a==0,visits)) && error("The following customers were not visited: $(filter(a->a==0,visits))")
   (abs(solution.cost-sum_cost) > 0.001) && error("Cost calculated from the routes ($sum_cost) is different from that passed as"*
                                                                                                  " argument ($(solution.cost)).") 
end

# read solution from file (CVRPLIB format)
function readsolution(app::Dict{String,Any})
   str = read(app["sol"], String)
   breaks_in = [' '; ':'; '\n';'\t';'\r']
   aux = split(str, breaks_in; limit=0, keepempty=false) 
   sol = Solution(0, [])
   j = 3
   while j <= length(aux)
      r = []
      while j <= length(aux)
         push!(r, parse(Int, aux[j]))
         j += 1
         if contains(lowercase(aux[j]), "cost") || contains(lowercase(aux[j]), "route")
            break
         end
      end
      push!(sol.routes, r)
      if contains(lowercase(aux[j]), "cost")
         if app["noround"]
            sol.cost = parse(Float64, aux[j+1])
         else 
            sol.cost = parse(Int, aux[j+1])
         end
         return sol
      end
      j += 2 # skip "Route" and "#j:" elements
   end
   error("The solution file was not read successfully. The file must be in the CVRPLIB format.")
   return sol
end

# write solution in a file
function writesolution(solpath, solution)
   open(solpath, "w") do f
      for (i,r) in enumerate(solution.routes)
         write(f, "Route #$i: ")
         for j in r
            write(f, "$j ") 
         end
         write(f, "\n")
      end
      write(f, "Cost $(solution.cost)\n")
   end
end

# write solution as TikZ figure (.tex) 
function drawsolution(tikzpath, data, solution)
   open(tikzpath, "w") do f
      write(f,"\\documentclass[crop,tikz]{standalone}\n\\begin{document}\n")
      # get limits to draw
      pos_x_vals = [i.pos_x for i in data.G′.V′]
      pos_y_vals = [i.pos_y for i in data.G′.V′]
      scale_fac = 1/(max(maximum(pos_x_vals),maximum(pos_y_vals))/10)
      write(f,"\\begin{tikzpicture}[thick, scale=1, every node/.style={scale=0.3}]\n")
      for i in data.G′.V′
         x_plot = scale_fac*i.pos_x
         y_plot = scale_fac*i.pos_y
         if i.id_vertex == 0 # plot depot
            write(f, "\t\\node[draw, line width=0.1mm, rectangle, fill=yellow, inner sep=0.05cm, scale=1.4] (v$(i.id_vertex)) at ($(x_plot),$(y_plot)) {\\footnotesize $(i.id_vertex)};\n")
            # Uncomment to plot without vertex id
            #write(f, "\t\\node[draw, rectangle, fill=yellow, scale=1.4] (v$(i.id_vertex)) at ($(x_plot),$(y_plot)) {};\n")
         else
            write(f, "\t\\node[draw, line width=0.1mm, circle, fill=white, inner sep=0.05cm] (v$(i.id_vertex)) at ($(x_plot),$(y_plot)) {\\footnotesize $(i.id_vertex)};\n")
            # Uncomment to plot without vertex id
            #write(f, "\t\\node[draw, circle, fill=white] (v$(i.id_vertex)) at ($(x_plot),$(y_plot)) {};\n")
         end
      end
      for r in solution.routes
         #=prev = r[1] # Uncomment (and comment below) to hide edges with the depot
         for i in r[2:end]
            e = (prev,i)
            write(f, "\t\\draw[-,line width=0.8pt] (v$(e[1])) -- (v$(e[2]));\n")
            prev = i
         end=#
         prev = 0
         for i in r
            e = (prev,i)
            edge_style = (prev == 0) ? "dashed,-,line width=0.2pt,opacity=.2" : "-,line width=0.8pt"
            write(f, "\t\\draw[$(edge_style)] (v$(e[1])) -- (v$(e[2]));\n")
            prev = i
         end
         write(f, "\t\\draw[dashed,-,line width=0.2pt,opacity=.2] (v0) -- (v$(prev));\n") 
      end
      write(f, "\\end{tikzpicture}\n")
      write(f, "\\end{document}\n")
   end   
end
