# Corrected source-bound static SAT on the fresh paired point

The [paired allocation](../paired-fresh-n17a1/target-panel.json) fixed
`[52411,72106]` before any outcome of the earlier SAT correctness pilot.
This [panel](panel.json) binds one corrected static-SAT arm on that public
point. It does not use a target scalar. The exact n17a1 source curve has
63 geometric factor-base points, 62 usable cofactor images and 29 signed
Frobenius columns. The full relation matrix and column logs must be built
from this arm's own natural queries, never imported from another arm.

The fresh StdRng08 ordinary seed `2026092955` produces 256 distinct
nonzero subgroup scalars with no overlap with the earlier 104 F5 queries,
32 static-SAT natural queries or 256 planned SAT correctness-pilot queries.
No query was chosen by a solver outcome. Collection stops only after full
rank and scalar replay of every column log, or after the frozen 256-attempt
cap. The descent seed is `2026092956`, with at most 64 target-dependent
`aG+bQ` queries. The source-receipted static CryptoMiniSat remains one
thread, one model, 1,000,000 conflicts and 120 seconds per query; the Rust
exporter has 60 seconds. All attempts, including capped and failed ones,
remain raw rows. No global wall or hard RSS cap is enforced by this local
runner; those limits are explicitly `null` in the workload envelope.

This revision classifies the solver's exit code 15 with
`s INDETERMINATE` at the declared conflict cap as
`CONFLICT_BUDGET_INCONCLUSIVE`, never UNSAT or an intrinsic solver crash.
It times source/model checking within each query separately from the
exporter/SAT work. The one-target online clock starts after all reusable
preparation is ready and stops **immediately** after independent scalar
replay, before final progress or summary writes. Its five exclusive phases
must add exactly to that interval. The run retains monotonic start/stop
timestamps, stop event, target-query attempts, source/model verification,
matrix rank trajectory, final LA and scalar certificate. A missing phase,
failed rank, unverified target or missing same-point rho leaves the speedup
unknown.

The exact candidate/workload/method/source manifests and run ID are sealed
before execution. A complete run still needs independent raw-evidence
replay. All four paired arms from the target allocation remain scheduled
irrespective of this arm's outcome; this isolated run alone establishes no
global winner or speedup.
