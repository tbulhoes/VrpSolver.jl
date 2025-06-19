using VrpSolver, JuMP, ArgParse

include("data.jl")
include("model.jl")
include("solution.jl")

function parse_commandline(args_array::Array{String,1}, appfolder::String)
   s = ArgParseSettings(usage="##### VRPSolver #####\n\n" *
                              "  On interactive mode, call main([\"arg1\", ..., \"argn\"])", exit_after_help=false)
   @add_arg_table s begin
      "instance"
      help = "Instance file path"
      "--cfg", "-c"
      help = "Configuration file path"
      default = "$appfolder/../config/TOP.cfg"
      "--ub", "-u"
      help = "Upper bound (primal bound)"
      arg_type = Float64
      default = 10000000.0
      "--batch", "-b"
      help = "batch file path"
   end
   return parse_args(args_array, s)
end

function run_top(app)
   println("Application parameters:")
   for (arg, val) in app
      println("  $arg  =>  $(repr(val))")
   end
   flush(stdout)

   instance_name = split(split(app["instance"], dirname(app["instance"]) * "/")[2], ".txt")[1]

   data = readTOPData(app["instance"])

   (model, x, y) = build_model(data)

   optimizer = VrpOptimizer(model, app["cfg"], instance_name)
   set_cutoff!(optimizer, app["ub"])

   (status, solution_found) = optimize!(optimizer)

   println("########################################################")
   if solution_found # Is there a solution?
      # for a in arcs(data)
      #    if get_value(optimizer, x[a]) > 0.0
      #       println("x[$a] = $(get_value(optimizer, x[a]))")
      #    end 
      # end
      sol = getsolution(data, x, get_objective_value(optimizer), optimizer)
      print_routes(sol)
      println("Cost $(sol.cost)")
   elseif status == :Optimal
      println("Problem infeasible")
   else
      println("Solution not found")
   end
   println("########################################################")
end

function main(args)
   appfolder = dirname(@__FILE__)
   app = parse_commandline(args, appfolder)
   isnothing(app) && return
   if !isnothing(app["batch"])
      for line in readlines(app["batch"])
         if isempty(strip(line)) || strip(line)[1] == '#'
            continue
         end
         args_array = [String(s) for s in split(line)]
         app_line = parse_commandline(args_array, appfolder)
         run_top(app_line)
      end
   else
      run_top(app)
   end
end

if isempty(ARGS)
   main(["--help"])
else
   main(ARGS)
end
