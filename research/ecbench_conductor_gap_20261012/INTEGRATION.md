# Integration and current validation

Snapshot 7b1dbaef preserves the measured evidence and reports before merging
upstream c2963dcff359814bf8203b946b88e9bb62431795. The library-test observations
in VALIDATION.md refer to the original measured baseline, not a green result
for every later upstream revision. Upstream's only change in the measured
ecbench/IC implementation paths is in a record-deserializer unit test; the
solver, accounting and query implementations used here remain unchanged.

The N83 validation issue is tracked at
[the named problem record](../../problems/n83-recorded-query-sampler-contract-mismatch.md).
It remains open and does not invalidate the two successful prime-field
session audits. Original logs, session files and snapshot hashes are retained.

The native analysis regenerates byte-identically. Both analyzer controls pass;
static site validation passes all 48 tests. The leaderboard, lab-browser and
panel-index render checks pass. All six PDF pages were visually inspected,
including the regenerated report with the validation limitation.

No new mathematical theorem is asserted, so the added upstream Lean-proof
rule has no new theorem to formalize in this observation-only change.
