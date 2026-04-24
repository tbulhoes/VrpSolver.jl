# Temporarily Disabled Functions

In the current release of **VrpSolver.jl**, some advanced branching and cut separation functionalities have been temporarily disabled to undergo maintenance and architectural improvements.

These features are planned to be reintroduced in future releases. If your current workflow strictly requires these methods, please refer to the legacy installation section below.

## Affected Functions

The following functions are currently unavailable in the latest version:

* `enable_packset_ryanfoster_branching!`
* `enable_resource_consumption_branching!`
* `add_strongkpath_cut_separator!`
* `add_permanent_ryanfoster_constraint!`

---

## Legacy Access & Compatibility

If you need to use these specific functions for your research, you must use the **older installation of VrpSolver.jl based on Docker**. This legacy environment contains the specific dependencies and core logic required for these routines to function.

### Instructions for Docker-based Installation

For step-by-step guidance on how to set up and run the legacy version using Docker, please visit the following link:

> [Insert link to the legacy installation instructions here]

---

