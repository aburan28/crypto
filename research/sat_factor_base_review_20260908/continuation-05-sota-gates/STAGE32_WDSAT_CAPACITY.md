# Stage 32: WDSat capacity correction

The completed Stage 26 n=59 cell exposed two source-specific WDSat assertions. Input ordinal 27 requires `__MAX_BUFFER_SIZE__` 32,588 and ordinal 34 requires 32,804, while the frozen Stage 26 build provided 32,264. The requirement calculation includes the strict-inequality margin used by WDSat.

Stage 32 changes only that macro to 32,804 and retains every other WDSat limit, source commit, compiler recipe, solver command, one-CPU policy, 120-second watchdog, source verifier, model parser, and exact point-witness validator. The corrected config SHA-256 is `4a73c3b5a14ded98f749d282b597594ae3ac05355a39cadcb9345a90bf735352`.

Only the two failed blind IDs are rerun. Their original `solver_error` rows remain in Stage 26. A corrected timeout is an inconclusive result; SAT is admitted only with an exact lifted point witness. Any repeated assertion, malformed terminal, or backend error fails Stage 32 rather than being substituted into the panel.

Source acquisition, the corrected WDSat build, packet verification, source verification, both solver runs, optional witness validation, one-CPU elapsed time, total CPU, and sampled process-tree memory are retained separately. The correction closes only the two-row capacity defect and does not provide Magma, external review, or SOTA evidence.
