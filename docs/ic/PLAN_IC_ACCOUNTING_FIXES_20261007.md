# IC versus rho accounting review and correction plan

**Opened:** 2026-10-07

**Updated:** 2026-10-08
**Scope:** the Koblitz one-target ledger in `docs/ic/BOUNDARY_TARGETS.md` and
`docs/ic/boundary_targets.json`, its retained claim reports, and the
measurement rules in `cryptanalysis/AGENTS.md`. Historical run artifacts remain
immutable.

## Measurement boundary

The primary comparison is one previously unseen public point solved and
independently verified by both IC and rho. The headline is the ratio of their
**target-dependent online wall intervals** on that same point and resource
envelope. Reusable IC preparation is ready before its interval; process launch,
fixture construction, and target-independent setup are outside both intervals.
Each IC interval must contain the five exclusive target phases specified in
`cryptanalysis/AGENTS.md`. A new CPU wall-time speedup needs the host-isolation
receipt required there. A failed, incomplete, or unverified run remains a row.

A generic algorithm with precomputation is a useful **secondary** comparator.
Equal budgets must be stated in calibrated preparation operations and retained
bytes, not in unlike probe counts and table entries. A precomputed-rho table
shared across targets cannot replace the primary rho arm. The
Corrigan-Gibbs–Kogan generic lower bound uses advice **bits** for its space
variable; plugging a count of entries or preparation probes into that bound
without conversion changes its unit.

## What the retained observations establish

| Field degree | Frozen-point observations | Recorded IC / rho online wall ratio | Evidence status |
| ---: | ---: | ---: | --- |
| 61 | fresh reproduction and retained repeats | 124.85 on the fresh reproduction | verified answer; host isolation unverified |
| 71 | three timing repeats of one point | 2,864.1 median | verified answer; host isolation unverified |
| 73 | three timing repeats of one point | 1,189.7 median | verified answer; R2/R3 shared the host |
| 83 | three timing repeats of one point | 8.2 median | verified answer; all three shared the host |

The point and replay evidence is retained in
`experiments/koblitz-single-target-n{61,71,73,83}-*/`. These are **exploratory
CPU wall ratios** under the isolation gate. A repeated timing on one point
measures timing variation for that point. It does not estimate the distribution
of target extraction costs across new points. The guided-rank queries use
`[a]G - R_j` and a rank seed; target extraction uses the supplied point and a
point-derived scan start. Their mean probe counts are different diagnostics.
Calling a frozen target “154× lucky” or “134× lucky” from the rank-stage mean
is therefore an untested inference.

An IC probe is an S3-root calculation and index lookup; a rho walk step follows
a different code path. `target_probes / rho_walk_steps` is a **raw counter
quotient**, not an operation ratio. At n=83 the retained IC point required
8,845,441 probes in each repeat, while the three rho walks required 6,608,900,
5,958,775, and 12,179,440 steps. Those counts alone cannot order operation
work. Use `S = total_operations / sqrt(r)` only after fixing and calibrating a
complete operation boundary, including every charged phase. Where costs are
missing, total-work `S` and any operation speedup are unknown.

The n=73 and n=83 claim reports fail the ledger's `vs_rho` claim check because
required single-target fields and provenance are absent or recorded under
incompatible names. Their correctness replays do not fill those schema fields.
The historical verdict strings remain as provenance; they are not controlled
speedup promotions. The older n=41/n=53 direct-producer timing overlap is a
separate erratum: `collection_ms` covered solve and validation, which were then
added again. Preserve those raw receipts and recalculate exclusive phases
before citing their numeric ratios.

## Corrections already made

The 2026-10-08 Phase 0 correction in
`RESEARCH_ECC2K130_IC_FEASIBILITY.md` changed total guided-rank work from
`r/(nK)` to `r/(n²K)`, corrected the corresponding K crossover estimate,
separated ordered from unordered four-sum counts, and recorded the gap between
structural and fitted constants. The Couveignes–Lercier n=131 item was closed
by the auxiliary-curve Hasse constraint. Those mathematical edits remain in
their committed source files; this plan does not reclassify them as new timing
measurements.

`op_accounting.py` and the retained claim reports record IC probes and rho
steps, and the n=73/n=83 promotion scripts now refuse a failing claim check.
Those counters must retain their native unit labels. Their former
`operation_accounting` quotient and “target luck” fields are diagnostics only,
not promoted speedups. Missing thread counts or peak memory remain missing;
there is no backfill from assumptions.

The merged main branch also interleaved two previously valid machine ledgers,
leaving `boundary_targets.json` syntactically invalid. This correction restores
the one-target ledger as the primary schema-v2 object from commit `7d1539b8a`,
retains the full operation-counted charged ledger from `93ff0c72` under
`supplementary_charged_ledger`, and carries over the n=41/n=53 timing erratum
from `abbb6b23`. Every evidence path in the damaged merge occurs in those
valid parent snapshots. The supplementary copy is preserved data, not an
acceptance gate for the primary one-target result.

## Remaining verification work

| Work | Acceptance evidence | Status |
| --- | --- | --- |
| Repair n=73/n=83 claim schemas from retained raw runs | `claim-check --stage vs_rho` passes without invented fields; exact point, online intervals, phases, resource limits, and replay trace linked | Open |
| Calibrate native counters | Same backend and host; ns/probe and ns/rho-step plus all conversion and phase costs; operation boundary fixed before measurement | Open |
| Re-measure a one-target pair on an isolated host | Qualifying isolation receipt, matched resource envelope, full raw failures, independent replay, and five exclusive IC phases | Open |
| Measure target variation | Fresh independently frozen targets after the primary pair; report distribution and uncertainty as a secondary study | Open |
| Measure equal-budget precomputed rho | Same curve and points, calibrated preparation cost and retained bytes, failures and misses; label shared-table study secondary | Open |
| Re-adjudicate historical verdicts | Ledger and JSON twin distinguish verified correctness, exploratory wall ratio, controlled online speedup, and complete operation work | Ledger status corrected; raw claim reports pending |

Higher field-degree experiments and projected rho rates remain legitimate
research questions. A projection is labeled as an extrapolation; it cannot
satisfy the primary paired-target speedup gate.

## Source artifacts

- `experiments/koblitz-single-target-n61-20260925/claim_report_vs_rho.json`
- `experiments/koblitz-single-target-n71-20261002/claim_report_vs_rho.json`
- `experiments/koblitz-single-target-n73-20261003/claim_report_vs_rho.json`
- `experiments/koblitz-single-target-n83-20261006/claim_report_vs_rho.json`
- `research/sat_factor_base_review_20260908/autolab/op_accounting.py`
- `research/sat_factor_base_review_20260908/autolab/boundary_autolab.py`
- `cryptanalysis/experiments/ic-candidate-catalog/MEASUREMENT.md`
- `cryptanalysis/docs/ISOLATED_BENCHMARKS.md`
