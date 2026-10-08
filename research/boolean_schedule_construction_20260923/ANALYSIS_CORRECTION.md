# Additive correction of an analysis-only naming mismatch

All 24 fixture processes finished with zero exit codes and completed verification.
The original producer then exited with status 1 because its analyzer generated
`constructiondispatch_*`, while the executed and declared dispatch arms are
`construction_dispatch_*`. The assertion failed before a scientific result was
written. The original source, protocol, raw observations, receipts and failed
analyzer are retained unchanged in `run_01`.

`run_01/EXECUTION_STATUS.json` records the observed failure and leaves its uncaptured
analysis execution timestamp null. A manifest seals that complete original set.
`corrected_analysis.py` repairs the dispatch name and accepts a separate output
directory. The formulas, thresholds, guards, input cohorts and measurements are
unchanged. `analysis_01` binds the sealed input manifest and both analyzer hashes,
and contains the corrected results and its own manifest. The correction is not a
new timing run and does not change the experiment's classification or eligibility.

The top-level `run.py` is a replay convenience updated to invoke the corrected
analyzer into a nested `analysis` directory on a fresh execution. Its historical
version remains in the original execution bundle. `reanalyze.py --out <NEW_PATH>`
replays the corrected verifier against the sealed original execution without
modifying it. Both entry points refuse an existing output directory.

Python regression checks retain the original failure, verify both artifact sets,
replay the corrected analysis byte-for-byte, and exercise null ceilings, observer
overhead, missing work, changed sources, censored results and reference omission.
These checks are producer validation, not external review or a speedup claim.
