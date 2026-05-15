## Parameters

This section lists the parameters available to the user. These parameters should be put to the configuration file. If a parameter is missing in the configuration file, the solver will use its default value, shown below.

In addition, the user may provide the following information to improve the solver performance.
- Designation of the first main resource which will be used for the bidirectional labelling.
- Priorities for branching strategies.
- Distance matrix between elementarity sets for defining initial $ng$-sets and for local search heuristic separation of $l$-row rank-1 cuts with $l\geq 5$.
  
---

### Main parameters

```
GlobalTimeLimit = 21474836  # time limit in seconds to solve the model
```
User time limit for the execution in seconds.

```
MaxNbOfBBtreeNodeTreated = 100000
```
Limits the total number of explored nodes in the branch-and-bound tree.

```
MaxDepthInBBtree = 100000
```
Limits the maximum depth in the branch-and-bound tree.

```
treeSearchStrategy = 1
```
Defines the exploration strategy for the primary branch-and-bound tree. Value `0` corresponds to breadth-first exploration (the open node with the smallest lower bound is considered next); value `1` corresponds to depth-first exploration (the open node with the largest depth is considered first).

```
OpenNodesLimit = 1000  # max. number of open nodes in the primary tree
```
Sets the maximum number of open nodes in the primary branch-and-bound tree. When this limit is reached, newly created nodes are pushed to the secondary branch-and-bound tree, which is always explored in depth-first manner. The solver returns to the primary tree when the secondary one becomes empty.

```
DEFAULTPRINTLEVEL = 0  # verbosity of the BaPCod output
```
Possible values are `-2` (no output except errors and important warnings), `-1` (reduced output, one line per 10 column generation iterations), `0` (standard output, one line per one column generation iteration). Positive values are not recommended as the output quickly becomes overwhelming.

```
solverName = CPLEX_SOLVER  # underlying LP and MIP solver
```
Can be set to `CLP_SOLVER` if BaPCod was appropriately configured. Note that certain features cannot be used with CLP solver, as the latter is not a MIP solver, and the overall performance may be degraded.

---

### MIP solver parameters

```
MipSolverMaxBBNodes = 2000000  # max. number of nodes for the MIP solver
MipSolverMaxTime = 360000      # time limit in seconds for the MIP solver
MipSolverMultiThread = 1       # number of threads for the MIP solver
```
These options are valid for the underlying MIP solver, which is used to solve the pricing problems (if parameterized for that), the restricted master problem as a MIP in the corresponding heuristic, and the enumerated master (if the pricing functor supports subproblem solution enumeration). If the value of `MipSolverMultiThread` is equal to `0`, the number of threads is determined automatically by the underlying MIP solver. Setting `MipSolverMultiThread` to `1` makes the whole solution process use a single thread.

```
MasterMipSolverRightHandSideZeroTol = --
```
Zero tolerance for the right-hand side values in the master MIP solver.

```
MasterMipSolverReducedCostTolerance = 1e-7
```
Reduced cost tolerance for the master MIP solver.

---

### Column generation parameters

```
MaxNbOfCgIterations = 100000
```
Maximum number of column generation iterations.

```
MaxNbOfStagesInColGenProcedure = 3  # number of col. gen. phases
```
Column generation phases are used to specify several algorithms for solving the pricing problems (only in the case the pricing functor is defined). Usually, during phase zero, the pricing problems are solved exactly, and the larger is the phase number, the lighter the heuristic algorithm applied. The stages are solved successively, from phase `MaxNbOfStagesInColGenProcedure - 1` to phase zero. The column generation procedure passes to phase `k-1` once phase `k` has converged.

```
ColumnCleanupThreshold = 10000
ColumnCleanupRatio = 0.66
```
Once the number of columns in the restricted master LP exceeds `ColumnCleanupThreshold`, only `ColumnCleanupRatio` part of them (with smallest reduced cost) remain, and the others are removed. The columns participating in the basis of the restricted master LP are never removed.

```
ReducedCostFixingThreshold = 0.9
```
Determines how often the reduced cost fixing procedure of the pricing functor is called. It is called if the current integrality gap is less than `ReducedCostFixingThreshold` part relative to the integrality gap when the reduced cost fixing procedure was called the last time. If the value is equal to `0.0`, no reduced cost fixing is performed. If the value is equal to `1.0`, reduced cost fixing is called after each convergence of the column generation procedure.

```
TerminateCgWhenRoundedDbCannotImprove = 0  # off
```
If activated, terminates column generation when the rounded dual bound cannot improve the incumbent solution.

```
UseObjScalingFact = 0  # off
```
If activated, scales the objective function by a constant factor to improve numerical stability.

```
BapCodReducedCostTolerance = --
```
Tolerance for the reduced cost in the BapCod framework.

```
BapCodIntegralityTolerance = --
```
Tolerance for integrality in the BapCod framework.

```
optimalityGapTolerance = 1e-6
```
If the relative gap between lower and upper bound is below this value, the column generation procedure terminates. Also, the node is pruned in the branch-and-bound tree if the relative difference between the global upper bound and the lower bound of the node is below this tolerance.

```
relOptimalityGapTolerance = 1e-9
```
Relative optimality gap tolerance. If the relative gap between lower and upper bound is below this value, the solution is considered optimal.

---

### Cut generation parameters

```
CutTailingOffThreshold = 0.015
CutTailingOffCounterThreshold = 3
```
These parameters are used to control the tailing-off condition of cut separation. The tailing-off counter is initialized with zero at the beginning of each branch-and-bound node. After a cut generation round, if the relative decrease of the integrality gap is smaller than the value of `CutTailingOffThreshold`, the tailing-off counter is increased by one. When the tailing-off counter reaches the value of `CutTailingOffCounterThreshold`, the tailing-off condition is activated: the cut separators with smaller priority are called if they are defined, or branching is performed.

```
CutCleanupThreshold = 1
CutCleanupRatio = 0.66
```
If the number of cuts reaches `CutCleanupThreshold`, all non-active cuts are removed from the restricted master LP. The `CutCleanupRatio` part of cuts (with largest violation) remain, and the others are removed.

```
BapCodCutViolationTolerance = --
```
Tolerance for cut violation in the BapCod framework.

```
ColGenSpRelaxationImprovementPriority = 0
```
Priority for improving the column generation subproblem relaxation.

---

### Stabilization parameters

```
colGenDualPriceSmoothingAlphaFactor = 1.0
colGenDualPriceSmoothingBetaFactor = 0.0
```
These two parameters correspond to parameters $\alpha$ and $\beta$ for dual price smoothing. The first parameter corresponds to Wentges smoothing and the second to directional smoothing. Value `0.0` means the technique is not used; value `1.0` means it is used with automatic parameter setting. Any value in $(0, 1)$ fixes the corresponding parameter to this value.

```
colGenStabilizationFunctionType = none
colGenProximalStabilizationRule = undefined
StabilFuncKappa = 1.0
```
The first parameter sets the stabilization function type: `0` (penalty function stabilization is not used), `2` (3-piecewise linear function is used), `3` (5-piecewise linear function). The second parameter switches between curvature mode (value `0`) and explicit mode (value `1`). The third parameter sets the value for parameter $\kappa$. The penalty function stabilization should be used with caution as it may deteriorate the column generation performance.

```
StabilizationMinPhaseOfStage = 0
```
Limits the stabilization use only to column generation phases with this number and above.

---

### Primal heuristic parameters

```
MaxTimeForRestrictedMasterIpHeur = -1
CallFrequencyOfRestrictedMasterIpHeur = 5
MIPemphasisInRestrictedMasterIpHeur = 1
PolishingAfterTimeInRestrictedMasterIpHeur = -1
```
The first parameter sets the maximum time in seconds for the MIP solver called to solve the restricted master MIP. The second parameter sets the frequency of the heuristic; its value should be `1` to call it at every node of the branch-and-bound tree. The heuristic is called only at the root node if the value of the second parameter is not positive. The last two parameters correspond to `CPX_PARAM_MIPEMPHASIS` and `CPX_PARAM_POLISHAFTERTIME` of the Cplex MIP solver. The restricted master heuristic cannot be used with CLP solver.

```
MaxNumEnumSolsInRestrictedMasterIpHeur = 5000
```
Maximum number of enumerated solutions used in the restricted master heuristic.

```
DivingHeurUseDepthLimit = -1
CallFrequencyOfDivingHeur = 1
```
The first parameter sets the maximum depth in the branch-and-bound tree for using the diving heuristic. If its value is negative, the diving heuristic is not used. The second parameter is equivalent to `CallFrequencyOfRestrictedMasterIpHeur`.

```
RoundingColSelectionCriteria
```
Determines the criteria for column selection for rounding. This parameter should be initialized with a chain of integers separated by spaces. Each integer corresponds to a certain criterion, which is used only if all previous ones could not select the column for rounding. The criteria are:
- `2` — highest priority (a column from a higher priority subproblem is preferred)
- `4` — smallest distance to the closest non-zero integer
- `5` — distance to the closest non-zero integer weighted by the column cost
- `6` — closest value to its round-up
- `9` — least column cost

```
FixIntValBeforeRoundingHeur = true
```
If set to `true`, any column with integer value in the solution will be fixed before rounding a non-integer column. Otherwise, integer columns will be ignored (and thus may take different values later in the dive).

```
MaxNbOfCgIteDuringRh = 5000
```
Limits the number of column generation iterations in each node of the diving heuristic.

```
MaxLDSbreadth = 0
MaxLDSdepth = 0
```
These parameters correspond to `maxDiscrepancy` and `maxDepth` in the diving heuristic with Limited Discrepancy Search. If their values are positive, they serve to control the diving heuristic with LDS.

```
DivingHeurStopsWithFirstFeasSol = false
```
If set to `true`, the diving heuristic will stop as soon as it finds the first feasible solution (diving for feasibility).

```
DivingHeurPreprocessBeforeChoosingVar = false
```
If set to `true`, the preprocessing will be launched after rounding of each candidate column (thus the diving will be slower). If preprocessing determines infeasibility, the candidate will be discarded and the next one will be considered. When this parameter is `false`, a dive is stopped if the preprocessing determines infeasibility.

```
StrongDivingCandidatesNumber = 1
```
If the value of this parameter is greater than `1`, the strong diving heuristic will be activated. This parameter corresponds to `maxCandidates` in the paper.

```
EvalAlgParamsInDiving = max#cand.=1 max#cg.iters=10000 min.lvl.sp.restr.=1 min#cut.rounds=0 max#cut.rounds=0 red.cost.fix&enum.=false tree.size.ratio=1
```
An optional parameter sequence to set the behaviour of the column and cut generation procedure at every node of the diving heuristic. The instantiation is similar to the parameters for strong branching phases. If this parameter sequence is empty, the same parameters are used as for the column and cut generation in the main branch-and-bound tree.

```
LocalSearchHeurUseDepthLimit = -1
MaxFactorOfColFixedByLocalSearchHeur = 0.8
MaxLocalSearchIterationCounter = 3
```
The first parameter sets the maximum depth in the branch-and-bound tree for using the local search heuristic. If its value is negative, the heuristic is not used. The last two parameters correspond to `fixRatio` and `numIterations` in the paper.

```
UseInitialPrimalHeur = false
```
If set to `true`, an initial primal heuristic is applied before the branch-and-bound procedure.

---

### Strong branching parameters

```
SimplifiedStrongBranchingParameterisation = true
```
If set to `true`, the strong branching parameters are set automatically using the four parameters below.

```
StrongBranchingPhaseOneCandidatesNumber = 100   # <p1>
StrongBranchingPhaseOneTreeSizeEstimRatio = 0.2 # <p2>
StrongBranchingPhaseTwoCandidatesNumber = 3     # <p3>
StrongBranchingPhaseTwoTreeSizeEstimRatio = 0.02 # <p4>
```
These parameters define the simplified strong branching setting. When `SimplifiedStrongBranchingParameterisation = true`, they set the number of candidates and tree size estimation ratios for phases one and two of strong branching.

```
SafeDualBoundScaleFactor = -1  # off
```
Numerically safe dual bound multiplier. When active, scales the dual bound to improve numerical stability.

```
StrongBranchingPhaseOne = max#cand.=100 max#cg.iters=0 min.lvl.sp.restr.=0 min#cut.rounds=0 max#cut.rounds=0 red.cost.fix&enum.=false tree.size.ratio=0.2
StrongBranchingPhaseTwo = max#cand.=3 max#cg.iters=10000 min.lvl.sp.restr.=1 min#cut.rounds=0 max#cut.rounds=0 red.cost.fix&enum.=false tree.size.ratio=0.02
StrongBranchingPhaseThree = max#cand.=1 max#cg.iters=10000 min.lvl.sp.restr.=0 min#cut.rounds=0 max#cut.rounds=10000 red.cost.fix&enum.=true tree.size.ratio=1
StrongBranchingPhaseFour = not active
```
Strong branching phase parameters. The parameter sequence is empty if the corresponding phase is not active. For each active phase, the parameters specify: whether the phase is exact, the maximum number of candidates evaluated, tree size ratio to stop, maximum number of column generation iterations (non-exact phases only), minimum column generation phase, minimum and maximum number of cut generation rounds, whether reduced cost fixing is performed, and the frequency of column generation output.

```
StrongBranchingUseHistory = 1  # on
```
If activated, uses the history of strong branching evaluations to avoid re-evaluating the same candidates.

---

### Debug output parameters

```
printMasterPrimalSols = 0
```
Controls printing of master primal solutions. Higher values produce more detailed output.

---

### VRPSolver parameters

```
RCSPstopCutGenTimeThresholdInPricing = 10  # seconds
```
Time threshold for the labeling algorithm to stop adding non-robust cuts.

```
RCSPhardTimeThresholdInPricing = 25  # seconds
```
Time threshold for the labeling algorithm to rollback to the state before the last cut round.

```
RCSPredCostFixingTimeThreshold = 100
```
Time threshold for the reduced cost fixing procedure in pricing.

```
RCSPnumberOfBucketsPerVertex = 25
```
Number of buckets per vertex.

```
RCSPdynamicBucketSteps = 1
```
Dynamic bucket step size adjustment: `0` = fixed (uses `RCSPnumberOfBucketsPerVertex`), `1` = adjusted dynamically but aggregated (same for all vertices), `2` = adjusted independently per vertex.

```
RCSPuseBidirectionalSearch = 2
```
Bi-directional search mode.

```
RCSPapplyReducedCostFixing = 1  # on
```
Bucket arc elimination via reduced cost fixing.

```
RCSPmaxNumOfColsPerIteration = 30
```
Maximum number of generated columns per iteration.

```
RCSPmaxNumOfColsPerExactIteration = 150
```
Maximum number of generated columns per exact iteration.

```
RCSPmaxNumOfLabelsInEnumeration = 500000
```
Maximum number of labels in the enumeration.

```
RCSPmaxNumOfLabelsInHeurEnumeration = 0
```
Maximum number of labels in the heuristic enumeration.

```
RCSPmaxNumOfEnumeratedSolutions = 5000000
```
Maximum number of paths in the enumeration.

```
RCSPmaxNumOfEnumSolutionsForMIP = 10000
```
Maximum total number of enumerated paths for MIP.

```
RCSPmaxNumOfEnumSolsForEndOfNodeMIP = 10000
```
Maximum total number of enumerated paths for MIP at the end of a node.

```
RCSPinitNGneighbourhoodSize = 8
```
Initial size of ng-sets.

```
RCSPmaxNGneighbourhoodSize = 8
```
Maximum size of ng-sets.

```
RCSPrankOneCutsMaxNumRows = 5
```
Limited-memory rank-1 cuts parameter: maximum number of rows.

```
RCSPrankOneCutsMaxNumPerRound = 100
```
Limited-memory rank-1 cuts parameter: maximum number of cuts per round.

```
RCSPrankOneCutsMemoryType = 0
```
Limited-memory rank-1 cuts parameter: memory type.

```
RCSPrankOneCutsLSnumIterations = 1000
```
Number of local search iterations for rank-1 cuts separation.

```
RCSPallowRoutesWithSameVerticesSet = true  # on
```
Avoids inserting routes with the same set of vertices from the same pricing into the Master LP.

```
RCSPredCostFixingFalseGap = 0  # off
```
Edge and path elimination based on a false gap $FG = (UB - LB) / \texttt{RCSPredCostFixingFalseGap}$.

```
LocArtVarInConvexityConstr = 0  # off
```
Add local artificial variables to convexity constraints: `0` = no, `1` = yes.
