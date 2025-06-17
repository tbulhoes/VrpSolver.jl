import Unicode

mutable struct Vertex
   id_vertex::Int
   pos_x::Float64
   pos_y::Float64
   demand::Float64
end

# Undirected graph
mutable struct InputGraph
   V′::Array{Vertex} # set of vertices (access with id_vertex + 1)
   E::Array{Tuple{Int64,Int64}} # set of edges
   cost::Dict{Tuple{Int64,Int64},Float64} # cost for each edge
end

mutable struct DataCVRP
   G′::InputGraph
   Q::Float64 # vehicle capacity
   depot_id::Int
   coord::Bool # instance with NODE_COORD_SECTION
   round::Bool # Is the distance matrix rounded?
end

customers(data::DataCVRP) = [i.id_vertex for i in data.G′.V′[2:end]] # return set of customers

# Euclidian distance
function distance(data::DataCVRP, arc::Tuple{Int64, Int64})
   e = (arc[1] < arc[2]) ? arc : (arc[2],arc[1])
   if haskey(data.G′.cost, e) # use already calculated value
      return data.G′.cost[e]
   elseif data.coord 
      u, v = arc
      vertices = data.G′.V′ 
      # array <vertices> is indexed from 1 (depot is vertices[1], customer 1 is vertices[2], and so on) 
      x_sq = (vertices[v+1].pos_x - vertices[u+1].pos_x)^2
      y_sq = (vertices[v+1].pos_y - vertices[u+1].pos_y)^2
      if data.round
         return floor(sqrt(x_sq + y_sq) + 0.5)
      end
      return sqrt(x_sq + y_sq)
   else
      return 0.0
   end
end

contains(p, s) = findnext(s, p, 1) != nothing

function readCVRPData(app::Dict{String,Any})

   str = Unicode.normalize(read(app["instance"], String); stripcc=true)
   breaks_in = [' '; ':'; '\n']
   aux = split(str, breaks_in; limit=0, keepempty=false)

   G′ = InputGraph([],[],Dict())
   data = DataCVRP(G′, 0, 0, false, !app["noround"])
    
   dim = 0
   for i in 1:length(aux)
      if contains(aux[i], "DIMENSION")
         dim = parse(Int, aux[i+1])
      elseif contains(aux[i], "CAPACITY")
         data.Q = parse(Float64, aux[i+1])  # the method parse() convert the string to Int64
      elseif contains(aux[i], "NODE_COORD_SECTION")
         data.coord = true
         j = i+1
         while aux[j] != "DEMAND_SECTION" 
            v = Vertex(0, 0, 0, 0)
            v.id_vertex = parse(Int, aux[j])-1 # depot is forced to be 0, fisrt customer to be 1, and so on
            v.pos_x = parse(Float64, aux[j+1])
            v.pos_y = parse(Float64, aux[j+2])
            push!(G′.V′, v) # add v in the vertex array
            j+=3
         end
      elseif contains(aux[i], "DEMAND_SECTION")
         j = i+1
         while aux[j] != "DEPOT_SECTION"
            pos = parse(Int, aux[j])
            G′.V′[pos].demand = parse(Float64, aux[j+1])
            j += 2
         end
         data.depot_id = 0
         break
      end
   end
   
   if data.coord
      # E = {{i,j} : i,j ∈ V′, i < j}  
      for i in customers(data)
         e = (data.depot_id, i)
         push!(G′.E, e) # add edge between depot and customer
         data.G′.cost[e] = distance(data, e)  
         for j in customers(data) # add edges between customers
            if i < j
               e = (i,j)
               push!(G′.E, e) # add edge e
               data.G′.cost[e] = distance(data, e)  
            end
         end
      end
   end

   return data
end

edges(data::DataCVRP) = data.G′.E # return set of edges
c(data,e) = data.G′.cost[e] # cost of the edge e
dimension(data::DataCVRP) = length(data.G′.V′) # return number of vertices
d(data::DataCVRP, i) = data.G′.V′[i+1].demand # return demand of i
veh_capacity(data::DataCVRP) = data.Q
nb_customers(data::DataCVRP) = length(customers(data))

function lowerBoundNbVehicles(data::DataCVRP) 
   sum_demand = 0
   for i in data.G′.V′
      sum_demand += i.demand
   end
   return Int(ceil(sum_demand/data.Q))
end

# return incident edges of i
function δ(data::DataCVRP, i::Integer)
   incident_edges = Vector{Tuple}()
   for j in 0:i-1 push!(incident_edges, (j, i)) end
   for j in i+1:(length(data.G′.V′)-1) push!(incident_edges, (i, j)) end
   return incident_edges
end
