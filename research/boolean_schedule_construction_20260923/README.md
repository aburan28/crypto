# Discovery cost bound for direct schedule construction

This standalone mathematical diagnostic follows PR #678. It measures whether
construction-only optimization has enough potential to deliver the requested
dramatic improvement. It does not implement that optimization or use holdouts,
curve inputs, imported targets, scalar recovery, production solvers or rho runs.

The diagnostic is complete. Under its unchanged-scan assumptions, six relevant
16-point groups and three dispatcher groups have guarded ceilings below 2.0.
The remaining groups are inconclusive, including every 64-point group. No
constructor optimization or performance promotion follows. See
[CONCLUSION.md](CONCLUSION.md) for the full accounting and limitations.

## Shared-kernel contract

The retained compiled 16-point and 64-point schedules are split into builders and
non-inlined scan functions. Profiled and unprofiled rebound arms use the exact same
builders and scan function. The diagnostic separates input encoding, schedule
construction, scanning and release. Their sum is bounded by the measured solve
cost. Complete totals also include result validation and all wrapper overhead.

The original 49 methods remain in the same binary. Original-policy verification
and the search-status reference are prepared outside arm timers but inside worker
receipts. Every profiled/rebound result must match the original policy's outcome,
first model, complete semantic work and trace. Independent SAT checks evaluate
the original equations; UNSAT requires the completed search reference. Caps remain
censored and block complete-result claims.

The dispatcher uses the same frozen width rule: 16-point schedules through n20,
64-point schedules above n20. It is not claimed to do less work than its selected
component. Extraction may alter generated code or data layout, so rebound timing
is measured against the original code, not assumed equal.

## Conditional cost ceiling

Under a fixed unchanged scan, even free construction leaves `T_new >= T_scan`.
The conditional maximum ratio to the strongest complete reference is therefore
`T_reference/T_scan`. This is an accounting identity under that assumption;
observed wall-clock values are not universal runtime lower bounds.

The diagnostic gives construction every measured advantage. From the scan timer
it subtracts the entire positive paired profile-minus-rebound total difference
and the maximum of 4,096 observed back-to-back clock pairs. A zero or negative
remaining scan estimate produces a null ceiling. This conservative adjustment
does not certify an absolute timing-error bound; the result remains conditional
on the recorded measurement and unchanged-scan model.

For a representative group, both 95% paired timing intervals—profile/rebound and
rebound/original—must fit within [0.8,1.25]. Otherwise it is inconclusive. A guarded
optimistic-ceiling upper bound below 2.0 rules out the construction-only target
for that measured group. A ceiling above 2.0 is merely an opportunity, not a gain.
No guard or threshold may be changed after seeing the measurements.

## Fixed discovery experiment

The protocol uses n12/16/20/24, discovery seeds 17 and 937, and all three retained
families. Fifty-five methods—49 retained and six diagnostic arms—run eight times
per fixture, giving 10,560 observations over 24 systems. Every deterministic random
order is followed by its reverse. The reference contains the 49 retained methods
and all three rebound arms. Its pointwise minimum is charged in complete totals.

Only n16/n20/n24 enter the prospective construction gate. Any guarded failure
blocks a universal constructor-only claim on this discovery grid. This does not
rule out changing the scan, data representation, inference algorithm or hardware;
those would be different mechanisms requiring new controls and measurements.

## Evidence and scope

Sixty retained Rust checks are joined by checks of extracted scans at every small
cap, exact original models/work on the discovery grid, and initial blocks against
direct equation values. Raw samples, per-worker receipts, compiler and source
hashes, clock observations and analysis are preserved in each fresh output folder.
Large temporary files use the checkout volume. Frozen artifacts are never edited.

The original execution completed all fixture processes but its analyzer failed
on a dispatch-arm naming mismatch. [ANALYSIS_CORRECTION.md](ANALYSIS_CORRECTION.md)
documents the unchanged execution and additive corrected analysis. New replays
use the corrected verifier and write their results under `analysis/`.

```sh
python3 research/boolean_schedule_construction_20260923/run.py --out /tmp/construction-replay
```

Choose an unused output path on a volume with enough space. The runner's compilation
scratch stays on the checkout volume. The WDSat/full-curve suite does not apply to
this general generated-system diagnostic. No speedup, exponent, generic-group
advance, production cost or rho crossover is inferred. All unmeasured broader
costs stay null; the original dramatic-gain objective remains open.
