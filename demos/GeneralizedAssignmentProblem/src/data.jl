import Base.show, Base.print

mutable struct DataGap
   machines::UnitRange{Int}
   tasks::UnitRange{Int}
   weight::Array{Int,2}
   cost::Array{Int,2}
   capacity::Array{Int,1}
end

function readGAPData(path_file::AbstractString)
   # STEP 1 : pushing data in a vector.
   data = Array{Int,1}()
   open(path_file) do file
      for line in eachline(file)
         for peaceofdata in split(line)
            push!(data, parse(Int, peaceofdata))
         end
      end
   end

   nbmachines = data[1]
   nbtasks = data[2]
   offset = 2
   cost = reshape(data[offset+1 : offset+nbmachines*nbtasks], nbtasks, nbmachines)
   offset += nbmachines*nbtasks
   weight = reshape(data[offset+1 : offset+nbmachines*nbtasks], nbtasks, nbmachines)
   offset += nbmachines*nbtasks
   capacity = reshape(data[offset+1 : offset+nbmachines], nbmachines)
   return DataGap(1:nbmachines, 1:nbtasks, weight, cost, capacity)
end

tasks(data::DataGap) = data.tasks
machines(data::DataGap) = data.machines
capacity(data::DataGap) = data.capacity
c(data::DataGap, t::Int, k::Int) = data.cost[t,k]
w(data::DataGap, t::Int, k::Int) = Float64(data.weight[t,k])

function show(io::IO, d::DataGap)
   println(io, "Generalized Assignment dataset.")
   println(io, "nb machines = $(length(d.machines)) and nb tasks = $(length(d.tasks))")
   println(io, "Capacities of machines : ")
   for m in d.machines
      println(io, "\t machine $m, capacity = $(d.capacity[m])")
   end

   println(io, "Ressource consumption of tasks : ")
   for j in d.tasks
      println(io, "\t task $j")
      for m in d.machines
         print(io, "\t\t on machines $m : consumption = $(d.weight[j,m])")
         println(io, " and cost = $(d.cost[j,m])")
      end
   end
end
