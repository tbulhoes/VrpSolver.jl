using Documenter, VrpSolver, JuMP

makedocs(
    modules = [VrpSolver],
    format = Documenter.HTML(prettyurls = false),
    sitename = "VRPSolver v2.0",
    pages    = Any[
        "Home"   => "index.md",
        #"Basic Model" => "basic_model.md",
        "Resources"   => "resources.md",
        "Methods" => "methods.md",
        "Parameters" => "parameters.md",
        "Disabled functions" => "disabled_functions.md"
    ]   

)

