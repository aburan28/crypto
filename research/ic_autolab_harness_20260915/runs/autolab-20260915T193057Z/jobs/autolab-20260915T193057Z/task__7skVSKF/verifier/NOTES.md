# Engineering Notes — Index Calculus Cost Reduction

## Session 1: batch_trials sweep

### Hypothesis

The incumbent uses `batch_trials=16` (relation batch size = 16). Each batch
generates 16 target/decomposition pairs in parallel, then runs a rank check
on the accumulated relation matrix. The phase share for
relation-verification/filtering/LA is ~29% of total cost.

A smaller batch (8) reduces surplus relation processing when the rank check
stops collection early, but increases the number of rank checks (more
repeated echelon work). A larger batch (24) amortizes rank checks better
but processes more surplus relations that never raise rank.

### Expected change

- `batch_trials=8`: ~3% reduction in relation phase cost from fewer surplus
  relations, partially offset by more rank checks. Net ~1-2% total cost
  improvement.
- `batch_trials=24`: ~2% cost increase from more surplus relations processed.

### Target

Probe both batch_trials=8 and batch_trials=24 against the incumbent at 16.
The 20% threshold is unlikely from this alone; this is a calibration probe
to map the batch-size cost surface for the next iteration.

### Plan

1. Try `batch_trials=24` first (larger batch, tests whether surplus processing
   dominates). If it wins, the relation phase is dominated by surplus work.
2. Then try `batch_trials=8` (smaller batch, tests whether rank-check
   overhead dominates). If it wins, the relation phase is dominated by
   repeated rank checks.
3. If neither wins, try `collection_window` with a window of 4 or 8 to
   reduce per-probe scan cost instead.
