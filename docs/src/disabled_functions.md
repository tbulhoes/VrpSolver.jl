# Temporarily Disabled Functions

In the current release of **VrpSolver.jl**, some advanced branching and cut separation functionalities have been temporarily disabled to undergo maintenance and architectural improvements.

These features are planned to be reintroduced in future releases. If your current workflow strictly requires these methods, please refer to the [legacy installation page](legacy_installation.md).

## Affected Functions

The following functions are currently unavailable in the latest version:

* `enable_packset_ryanfoster_branching!`
* `enable_resource_consumption_branching!`
* `add_strongkpath_cut_separator!`
* `add_permanent_ryanfoster_constraint!`

---

## Legacy Access & Compatibility

If you need to use these specific functions for your research, you must use the older Docker-based installation of VrpSolver. See the [Legacy Installation (Docker)](legacy_installation.md) page for the legacy files and step-by-step setup instructions.
