using VrpSolver, JuMP, ArgParse

include("data.jl")
include("model.jl")

function parse_commandline(args_array::Array{String,1}, appfolder::String)
   s = ArgParseSettings(usage="##### VRPSolver #####\n\n"*
	   "  On interactive mode, call main([\"arg1\", ..., \"argn\"])", exit_after_help=false)
   @add_arg_table s begin
      "instance"
         help = "Instance file path"
      "--cfg", "-c"
         help = "Configuration file path"
         default = "$appfolder/../config/GAP_Classic.cfg"
      "--ub","-u"
         help = "Upper bound (primal bound)"
         arg_type = Float64
         default = 10000000.0
      "--batch","-b" 
         help = "batch file path" 
   end
   return parse_args(args_array, s)
end

function run_gap(app::Dict{String,Any})
   println("Application parameters:")
   for (arg,val) in app
      println("  $arg  =>  $(repr(val))")
   end
   flush(stdout)

   data = readGAPData(app["instance"])

   (model, x) = build_model(data)
   optimizer = VrpOptimizer(model, app["cfg"])
   set_cutoff!(optimizer, app["ub"])

   (status, solution_found) = optimize!(optimizer)

   println("########################################################")
   if solution_found # Is there a solution?
      println("Objective value : $(get_objective_value(optimizer))")
      for k in data.machines
         tasks_k = []
         for t in data.tasks
            if get_value(optimizer, x[t,k]) > 0.99
               push!(tasks_k, t)
            end
         end
         if !isempty(tasks_k)
            print("machine $(k):")
            for t in tasks_k
               print(t, " ") 
            end
            println("") 
         end
      end
   elseif status == :Optimal
      println("Problem infeasible")
   else
      println("Solution not found")
   end
end

function main(args)
   appfolder = dirname(@__FILE__)
   app = parse_commandline(args, appfolder)
   isnothing(app) && return
   if app["batch"] != nothing
      for line in readlines(app["batch"])
         if isempty(strip(line)) || strip(line)[1] == '#' 
            continue
         end
         args_array = [String(s) for s in split(line)]
         app_line = parse_commandline(args_array, appfolder)
         run_gap(app_line)
      end
   else
      run_gap(app)
   end
end

if isempty(ARGS)
   main(["--help"])
else
   main(ARGS)
end
