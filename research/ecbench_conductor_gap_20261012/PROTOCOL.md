# Large-prime conductor-gap IC follow-up

Frozen before execution, 2026-10-12 UTC. Classification: **boundary control**;
operation-count diagnostic. This extends the archived D = -19 CM-pair study,
not its unrelated small-prime-power fiber experiments.

## Requested scope and checkpoint

Resume measuring ECDLP hardness across large conductor gaps. The archived
`research/conductor_gap_cm_pairs_20261010` has verified IC measurements at
ell = 263 and 8209, but only generic solvers at ell = 131101. Extend the
implemented IC comparison to that existing 37-bit pair, and replicate the
29-bit comparison with fresh targets. Broader hardness and cryptographic-scale
conclusions remain open; this bounded experiment does not discharge them.

## Inputs, algorithms and accounting

Reuse exact archived models, subgroup orders and generators. Select the crater
and **first listed** floor in each archived spec, before measuring. These
are two fixed pairs, not a random sample of isogeny classes or floor curves.
The prior Sage point-count certificates remain historical evidence; Sage is
unavailable here. Native workload validation and exact log checks are required.

At each ell use eight independent planted one-target workloads per curve,
target seed 2026101201; one measured round and one warm-up, order `alternate`,
algorithm seed 20261012. Target generation hashes the curve model, so equal
target indices on different curves do **not** mean equal scalars or transported
targets. Cross-curve inference uses independent target samples. No map is
computed and no transfer cost is assumed to be zero.

Arms: `rho.negation`, an identical `rho.negation` A/A control, and
`ic.pipeline` with the two-summand negation-folded MITM oracle, walked queries,
and `structured-gauss`. Factor-base sizes are fixed before execution: 256
abscissae at 29 bits, 2048 at 37 bits. Actual points, columns and table sizes
come from records. This is a new frozen configuration, not a timing comparison
against the historical 64-column configuration. No parameter tuning after
reading results.

Each execution is cold: base selection, table construction, queries, relation
handling, linear algebra, recovery and verification are retained. Primary unit
is the harness's counted group-addition-equivalent S = GAE / sqrt(r), with
its recorded unpriced-work/lower-bound flags. Phase counters supplement S;
they do not replace it. Ratios of lower bounds are diagnostics, not full-cost
bounds or speedups. The generic reference is matched one-target rho; its floor
is sqrt(pi/4) in the harness unit. No asymptotic fit from two sizes.

## Prediction and frozen analysis

Prediction: IC floor/crater mean counted S is near one at both gaps.
Compute mean S separately for each curve and arm; floor/crater is the ratio
of those means. Bootstrap each curve's eight workloads independently (20,000
resamples, seed 20261012), preserving all rounds of a sampled workload.
Report 95% and 99% percentile intervals; use the 99% interval for the two-pair
lead decision. Also report a two-sided label-permutation test of the absolute
log ratio (20,000 permutations, plus-one p-value). A/A is an exact per-workload
count equality check. A finite-sample rho crater/floor difference is a control
observation, not automatically a harness defect or a conductor effect.

A repeated IC lead requires both 99% intervals to exclude one in the same
direction, an effect exceeding 10% in both pairs, and exact A/A equality. A
lead still requires a new independent pair, representation controls and
independent validation before a structural conclusion. Otherwise record the
uncertainty; an interval containing one is not an equivalence proof.

Retain all failed, timed-out, exhausted and interrupted runs. Any measured
failure makes that arm's inferential comparison inadmissible; do not silently
analyze the successful subset. Warm-ups remain in raw evidence but not in
the primary estimates. No survival analysis or failure-cost substitution is
introduced after observing outcomes.

## Resource limits and execution

One child at a time, RAYON_NUM_THREADS=1; 1200-second per-execution watchdog,
100,000,000 relation-query cap. No spending cap: watchdogs protect the host.
Attempt native ecbench CPU reservation. Below-L2 runs may supply operation
counts only; wall time is descriptive and is not used for an advantage claim.
Do not run builds or tests during a session. Audit all deterministic measured
runs by exact replay; keep the audit receipts and all refused starts.

## Deliverables and limits

Committed protocol/specs; source and binary hashes; exact commands and host;
sealed sessions and replay receipts; native analysis JSON/Markdown; report,
vector graph and inspected PDF; dashboard coverage; review-ready PR.
No hypothesis-state promotion or claim of intrinsic/cryptographic hardness.
The 57-bit ECC2K-130 conductor-gap endpoints remain unconstructed.

The repository's `conductor check` was attempted before editing; the CLI and
`.conductor/WORKFLOW.md` are not installed in this cloud checkout. Changes
are confined to a new study directory and a native analysis example, with
dashboard publication paths added explicitly after results. This is a tool
availability deviation, not a scientific exemption or a competing-file claim.
