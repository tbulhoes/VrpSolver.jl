using Documenter, VrpSolver, JuMP

makedocs(;
    modules = [VrpSolver],
    format = Documenter.HTML(; prettyurls = false),
    sitename = "VRPSolver v2.0",
    pages = Any[
        "Home" => "index.md",
        "Methods" => "methods.md",
        "Parameters" => "parameters.md",
        "Custom resources" => "custom_resources.md",
        "Disabled functions" => "disabled_functions.md",
        "Legacy installation (Docker)" => "legacy_installation.md",
    ],
)

