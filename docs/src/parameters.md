## Parameters

Table below list the VRPSolver parameters available to the user. In addition the user may provide the following information to improve the
solver performance.
- Designation of the first main resource which will be used for the bi-directional labelling.
- Priorities for branching strategies.
- Distance matrix between elementarity sets for defining initial $ng$-sets and for local search heuristic separation of $l$-row rank-1 cuts with $l\geq 5$.

|                 Parameter                 |                Description                |   Default value  |
|:-----------------------------------------:|:-----------------------------------------:|:----------------:|
|              GlobalTimeLimit              |     User time limit for the execution     | 21474836 seconds |
|          MaxNbOfBBtreeNodeTreated         |    Max number of nodes in the b&b tree    |      100000      |
|    RCSPstopCutGenTimeThresholdInPricing   | Time threshold for the labeling algorithm |    10 seconds    |
|       RCSPhardTimeThresholdInPricing      | Time threshold for the labeling algorithm |    20 seconds    |
|        RCSPnumberOfBucketsPerVertex       |    Calculation of step size for buckets   |        25        |
|           RCSPdynamicBucketSteps          |    Calculation of step size for buckets   |      1 (on)      |
|         RCSPuseBidirectionalSearch        |           Bi-directional search           |         2        |
|         RCSPapplyReducedCostFixing        |           Bucket arc elimination          |      1 (on)      |
|     RCSPmaxNumOfColsPerExactIteration     | Max. # of generated columns per iteration |        150       |
|        RCSPmaxNumOfColsPerIteration       | Max. # of generated columns per iteration |        30        |
|        StabilizationMinPhaseOfStage       |      Minimum stage for stabilization      |  0 (everywhere)  |
|      RCSPmaxNumOfLabelsInEnumeration      |    Max. # of labels in the enumeration    |       1000000       |
|      RCSPmaxNumOfEnumeratedSolutions      |     Max. # of paths in the enumeration    |       1000000       |
|      RCSPmaxNumOfEnumSolutionsForMIP      |  Max. total # of enumerated paths for MIP |       10000       |
|   RCSPmaxNumOfEnumSolsForEndOfNodeMIP     |   Max. total # of enumerated paths for MIP at the end of a node  |      10000        |
|        RCSPinitNGneighbourhoodSize        |          Initial size of ng-sets          |         8        |
|         RCSPmaxNGneighbourhoodSize        |          Maximum size of ng-sets          |         8        |
|       RCSPrankOneCutsMaxNumPerRound       |   Limited-memory rank-1 cuts parameter 1  |        100       |
|         RCSPrankOneCutsMaxNumRows         |   Limited-memory rank-1 cuts parameter 2  |         5        |
|         RCSPrankOneCutsMemoryType         |   Limited-memory rank-1 cuts parameter 3  |         2        |
|   RCSPallowRoutesWithSameVerticesSet      |   Avoid inserting routes with the same set of vertices from the same pricing into the Master LP.  |         true (on)        |
|           RCSPredCostFixingFalseGap              |   Edge and path elimination based on a false gap FG=(UB-LB)/RCSPredCostFixingFalseGap  |         0 (off)        |
|   RCSPmaxNumOfLabelsInHeurEnumeration     |        Max. number of labels in the heuristic enumeration.   |         -- (MaxTimeForRestrictedMasterIpHeur must be active)        |
|           CutTailingOffThreshold          |   Cut generation tailing-off parameter 1  |        0.02 (2%)        |
|       CutTailingOffCounterThreshold       |   Cut generation tailing-off parameter 2  |         3        |
|          SafeDualBoundScaleFactor         |   Numerically safe dual bound multiplier  |     -1 (off)     |
|  StrongBranchingPhaseOneCandidatesNumber  |   Strong branching parameter for phase 1  |        100       |
| StrongBranchingPhaseOneTreeSizeEstimRatio |   Strong branching parameter for phase 1  |        0.3       |
|  StrongBranchingPhaseTwoCandidatesNumber  |   Strong branching parameter for phase 2  |         3        |
| StrongBranchingPhaseTwoTreeSizeEstimRatio |   Strong branching parameter for phase 2  |        0.1       |
|      MaxTimeForRestrictedMasterIpHeur     |        Restricted master heuristic        |     -1 (off)     |
|          DivingHeurUseDepthLimit          |        Diving heuristic (with LDS)        |     -1 (off)     |
|                MaxLDSdepth                |        Diving heuristic (with LDS)        |         0        |
|               MaxLDSbreadth               |        Diving heuristic (with LDS)        |         0        |
|   CallFrequencyOfRestrictedMasterIpHeur   |  Frequency (in terms of nodes) in which the restricted master heuristic is invoked (it is first called at the root node)  |   -- (MaxTimeForRestrictedMasterIpHeur must be active) |
| MIPemphasisInRestrictedMasterIpHeur       |        Set the MIP solver to concentrate on improving primal bound at the expense of improving lower bound.        |         -- (MaxTimeForRestrictedMasterIpHeur must be active)        |
| MaxNumEnumSolsInRestrictedMasterIpHeur    |  Max. number of enumerated solutions in the restricted master heuristic.   |         -- (MaxTimeForRestrictedMasterIpHeur must be active)        |

