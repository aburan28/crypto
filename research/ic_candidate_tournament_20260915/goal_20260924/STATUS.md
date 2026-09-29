# September 24 bounded IC goal

Status: the archived reference panel has completed under the frozen protocol:
1,290/1,290 native/profile pairs, three IC sources and eighteen rho configurations.
[PR 761](https://github.com/aburan28/crypto/pull/761) merged the implementation at
`e9b540eef8775f8d0d65e24c319174d8fefa7460`. The
[qualification report and durable evidence](reference-qualification/README.md)
select `pairinv` for IC cold instructions and online time, `rho_incumbent_4` for
rho cold instructions and `rho_pairinv_4` for rho online time. PR 765 merged the evidence at `c2ca50318d4b9a3dc4c259e59106f2717553da31`,
including exact Linux replay of 22 retained evidence sets and every exported table.
Independent macOS receipt replay passed all 1,290 jobs, with two one-ULP
derived-summary differences documented in the report.
No reference selection is a promotion. The three-attempt pair-table improvement
budget is complete without a promoted challenger: see
[BOUNDED-GOAL-RESULT.md](BOUNDED-GOAL-RESULT.md).
[Round one](improvement/round1/README.md) (workflow `36140265516`),
[round two](improvement/round2/README.md) (seed `2026092552`) and
[round three](improvement/round3/README.md) (workflow `36463687634`, seed
`2026092553`) each retained the qualified `pairinv` incumbent after full
verification. Selected challengers `stop6`, `stop5_word` and `stop7_word` all
failed the 0.8 cold gates and/or familywise online rules. Altogether
10,203/10,203 native/profile pairs were verified. Do not redispatch any of the
three registered attempts or retune on their confirmation/replay points. The
next registered comparison is
[generic-backend-qualification](generic-backend-qualification/README.md)
(seed `2026092901`); its measurement is pending.
Canonical admission is merged in both drivers; the earlier
[driver controls](driver-admission/README.md) preserve their fixed-vector scope.
Public-point input and single-target native intervals are implemented in the
[follow-up controls](public-inputs/README.md); PR 748 merged at
`4a898fcd3464b71439bd9451bd466a6ec217ddc9`, with 39/39 IC and 13/13 rho
native/profile pairs passing final Linux validation and independent replay.
Foundation: [PR 704](https://github.com/aburan28/crypto/pull/704), merged at
`c04879dbb305423d4e7bad9b9e9dd2188f000a35`.
Canonical record/field/accounting work: [PR 718](https://github.com/aburan28/crypto/pull/718),
merged at `e3b0a1b6bc823ed2205c8333bc43474fe5f1a0b4`.
Optimized producer work: [PR 733](https://github.com/aburan28/crypto/pull/733).
The tested implementation head `6b5a402558c2c883baa52c2ef585f99886108167` passed
all 39 Linux native/profile pairs. The [admission report and durable evidence](producer-admission/README.md)
retain independent transport replay and the earlier failed controls. This
qualifies the instrumentation on fixed vectors; it does not select an incumbent.

## Objective and stopping rule

Complete at most three improvement rounds. Target at least 20% lower complete
cold instruction cost **and** native cold time than the strongest qualified
archived IC incumbent, on the declared workload only. Finish with a qualifying
winner or an audited finding that none qualified in those rounds. Scaffolding,
an incomplete campaign or exhausted implementation time is not completion.

Before measurement, freeze at least five exact curve/subgroup cells, at least
60 fresh single-target confirmation cases, three process repetitions per target,
resources, accounting unit, floor/reference and stopping rules. Up to 16 complete
pipelines per round; retain six diverse candidates including an exploration
slot. Measure actual combinations, including mechanisms with individually slower
parents. One challenger per round enters confirmation; never retune on it.

Predeclare a familywise 95% rule across at most three confirmation attempts.
Promotion requires both 20% point improvements, uncertainty excluding regression,
no per-cell regression above 10%, independently certified answers for every
target, and passing confirmation and replay. The precise allocation, estimator,
familywise rule and panel must be sealed before the first new measured round;
this status note does not substitute for that protocol.

The current user measurement contract makes **single-target online wall time**
the headline metric, excluding reusable preparation and fixture generation. The
complete cold instruction/time goals above are additional acceptance gates.
Report both boundaries explicitly; neither batch amortization nor a cheap solver
stage substitutes for the one-target result.

## Baseline inventory and archived provenance

| Source | Durable identity | Why retained |
|---|---|---|
| Round 0020 `both` | Source manifest `563eb460f29d9ef09a2adbde4770a466d3dc16566f5f08f1e6d8238325b104d4` | Last formally promoted single-target winner |
| Round 0023 `scaled` | Source manifest `55154f73c35b1f55240b35fd1a5e8e2488c6df41444b5114949e41741c0c3db1` | Improved the archived incumbent; nonpromotion under the old rho objective does not disqualify it as a stronger IC reference |
| Round 0024 `pairinv` proposal | `campaign_20260916/round24-pairinv.patch` over `scaled` | Retained source change and public equivalence check; inspect and qualify before choosing a baseline |

The archive manifest and restore tool are in `../evidence/`. Archive bytes:

- Round 0020 SHA-256 `944a4bd2d6e54fc95fc7bb108ef09a5a1319da45882cf5565c30a8d30f653c2f`.
- Round 0023 SHA-256 `e2aa7111bb03ae606a3af6f249c9bd18beb3df1cc3ca4f813812010336dd47da`.

Both archives restored successfully during this audit. Restoration checks hashes;
it is not a fresh correctness/performance qualification. Preserve the original
sealed sources. Instrument derived copies and measure observer overhead and
equivalence against originals on public development fixtures. Also inspect the
current `icx` engine before claiming the strongest compatible incumbent.

Correction from the materializer's source-hash check: the historical winner
summary's `source_root` names round 0020 `both`, but its top-level source hash
`69de47e3...` belongs to that round's incumbent. The archived candidate registry
and independently hashed `both/source-manifest.json` identify `both` as
`563eb460...`. Use the verified candidate source, not the inherited summary hash;
the historical summary remains unchanged as evidence of the discrepancy.

Source review found two qualification hazards in archived `scaled`: the tiny
fast path dispatches before the configured `linear_algebra` selection and
actually uses incremental Gaussian elimination, and `solve_target` can emit an
empty direct-collision witness. The current independent oracle rejects such a
witness as IC. Resolve/report the actual dispatched method and qualify against
the stricter admission rule. A nominal LA configuration sweep is not evidence
that multiple LA solvers ran. The current `icx` runner constructs planted target
scalars; a public-target adapter is required before that path joins this panel.

## Next gates

1. Canonical record and base census regression/CI: implemented in `identity.py`,
   `measurement.py`, `test_records.py` and `ci_smoke.py`; merged in PR 718.
2. Optimized archived producers now export eleven exclusive phases, ordinary-query
   outcomes/rank and matrix diagnostics; 39/39 pairs independently replay. Combined
   old labels remain unknown under the new schema. Admission is now wired into
   both drivers; the [driver control protocol](driver-admission/PROTOCOL.md)
   covers canonical records, online timing, retained failures and frozen audit. The public-point
   and single-target native timing follow-up passes local controls, Linux
   integration and transported evidence replay. Rho reusable arithmetic/Frobenius
   preparation is excluded from its online interval and retained in cold cost.
3. The bounded reference-quality checks and five-cell development qualification
   have executed; see the report above. Bind its selected complete sources and
   rho settings to the new protocol before any held-out data.
4. Seal the panel and familywise protocol, then run bounded rounds. Keep all
   failures and source/fixture/profiler artifacts; archive them durably, update
   the existing scoreboard, and merge implementation/evidence PRs.

The local development machine is macOS arm64. Calibrated measurements require
the existing Linux amd64 / Valgrind 3.22.0 workflow; native local timings cannot
substitute for that accounting model. The new measurements are development reference selection, not an improvement-round claim.

## Reference-quality checks before comparison

Review of the prepared `scaled` source identified two concrete issues. PR 761
completed the allocator correction and independent row-arithmetic controls; it
did not optimize the row kernel:

- In `examples/ic_tournament_worker.rs`, the arena region guarantees 4096-byte
  alignment but allocation rounds only the offset for arbitrary requested
  alignments. Harden requests above that guarantee with a system-allocator
  fallback or correct absolute-address alignment, and exercise the release
  allocator controls. Keep old measured source snapshots immutable.
- In `src/cryptanalysis/koblitz_tiny_ic.rs`, `Echelon::push` performs modular
  arithmetic across every column, including already-zero prefixes. Review
  trailing-column reduction and bounded-modulus arithmetic against independent
  scalar-field rank/solution checks before treating this kernel as the strongest
  compatible reference. A faster kernel must be measured; source inspection
  alone does not establish a gain. These small matrices do not by themselves
  justify replacing Gaussian elimination with a large sparse solver.

Inspect the current compatible `icx`/rho paths as well as restored sources.
Freeze a cross-campaign target exclusion set and the familywise confirmation rule
before enabling promotion. The bounded runner now binds both selected rho sources/settings and the IC
incumbent to the accepted qualification digest. The next deliverable is its
executed bounded rounds; more integration controls alone will not complete this goal.


Reference qualification now has measured development evidence. `pairinv`'s online
ratio to the old incumbent is 0.9781 [0.9249, 1.0302], its cold instruction ratio
is 0.9845 and its cold native ratio is 1.0064. This is not a 20% gain. Larger
requested rho widths clip to the same effective width on several cells; the
report preserves those counts. The existing exact evaluator remains unchanged,
and Linux archive replay checks both raw receipts and derived selection. The first-round runner validates those bindings and exclusions before preparing
fresh targets. Execute the registered diversified pipeline budget after its
implementation PR passes; do not infer an improvement from these controls.


## Generic query admission checkpoint

The [generic query accounting protocol](generic-query-accounting/PROTOCOL.md)
preserves every attempted collection/descent query, typed frontend outcomes and
solver counters, plus terminal failed descents. The worker exports these records
and retains actual attempted matrix solves. Independent group and bounded
negative-answer replay checks accounting only; complete generic scientific
admission still needs exact source/base/matrix binding and
exclusive public-target timing. This is not an improvement round. The incumbent
remains selected and two rounds remain under the frozen goal protocol.

The [generic supplied-point follow-on](generic-public-inputs/RESULTS.md) separates
fixture creation from measured jobs and places the outer online interval after
reusable IC/rho preparation through independent scalar replay. Its 37 final local
worker controls pass, including seven intended preparation failures with null
online intervals. Combined legacy phase dumps remain unqualified for scientific
cost comparison; full generic admission and the remaining two rounds are open.

PR 803 merged at `62ef21ec1e083f197593edbe1309c5ddf60b7789` after all applicable
checks passed, including Linux integration and strict archived-round replay.
The [independent query-law controls](generic-query-law/RESULTS.md) now replay
7,436 pinned Rust RNG/probe values, all 35 archived IC reports, and 47 fresh
controls (40 complete, seven intentionally incomplete). Wrong seeds, batch
partitions and collection/descent rules are rejected even when group equations
remain valid. This is accounting admission, not a new measured improvement round.

PR 805 merged the query-law checks at
`c7c2922c116b2ec3ca84a2066a8b9c63a782a39d`, with all applicable checks passing.
The [exclusive generic phase follow-on](generic-exclusive-phases/RESULTS.md)
now passes two retained local 147-pair panels (126 complete and 21 deliberately
incomplete pairs per panel). It separates query/PDP/checking/matrix/LA/descent
work and independently checks native clock closure. Strict sessions reject
phase changes on another thread. Linux instruction closure is exercised by
the PR integration checks. Shared-host mode ratios remain too unstable to
qualify overhead or comparative performance. This accounting work consumes no
round; exact generic admission and optimized-reference qualification precede
the remaining two rounds.


## Generic five-cell readiness

The [registered readiness check](generic-reference-readiness/README.md) passes
15/15 jobs: ten pair-table IC solves with dense/sparse relation-LA policies and
five rho solves on the five development cells. All raw reports and canonical
records independently replay after fresh archive extraction. The two n23a1 IC
runs retain unresolved queries; no larger proved-negative PDP claim is admitted.
This is admission readiness only, pending the implementation/evidence PR and
parent PR 822. Calibrated generic reference and instrumentation qualification
still precedes comparative ranking. Preserve all five exposed points in
`generic-reference-readiness/fixtures.json` as exclusions for the next registered
confirmation panel; keep the sealed round-one history unchanged. One round is
closed, no challenger qualified, and two remain.

## Three-round closeout and next registration

Rounds two and three are archived
([PR 893](https://github.com/aburan28/crypto/pull/893),
[PR 916](https://github.com/aburan28/crypto/pull/916)). The audited negative
result of the pair-table campaign is
[BOUNDED-GOAL-RESULT.md](BOUNDED-GOAL-RESULT.md). The next preregistered
comparison is
[generic-backend-qualification](generic-backend-qualification/PROTOCOL.md):
qualify generic F4/F5/SAT and sparse relation-LA complete pipelines against the
optimized incumbent and strong rho on fresh points that exclude every exposure
from all three sealed rounds. Panel byte SHA-256
`83c640a03b4239b918851f6f1b8450e2fe27dc710fb99f306e3481a99d8875cf`, seed
`2026092901`. Its one measured dispatch was canceled at the six-hour job cap;
artifact upload also failed, leaving completion and comparative costs unknown.
The [censored result](generic-backend-qualification/RESULT.md) is not a family
qualification or solver-performance verdict. A fresh protocol must exclude all
potentially exposed points and retain partial results before a job timeout.
The [second registration](generic-backend-qualification-v2/PROTOCOL.md) freezes
seed `2026092902`, panel SHA-256
`d283a869b0412228d1c66260fdfd8f387d7243bd15456c7febf3c46ee5da27a8`,
one process on each of 25 distinct points and 250 trial slots. It excludes the
25 reconstructed first-run points and keeps the same source-bound F4/F5, SAT,
incumbent and matched-rho arms. Its one permitted dispatch is
[Actions run 36580669479](https://github.com/aburan28/crypto/actions/runs/36580669479)
(attempt one); never redispatch this seed. On 2026-09-29 the measured step
recorded a 300-minute timeout, `pack_partial_campaign.py` failed, and the
upload step found no campaign archive (only the earlier pack-smoke artifact).
The workflow object was still `in_progress` when those step conclusions were
read, so this note does not replace a terminal closeout. Qualification and
cost rows stay unknown. A timed-out run with no retained bundle is operationally
censored for family verdicts.

Post-registration source audit (PR
[#952](https://github.com/aburan28/crypto/pull/952)): the pinned `m=3` Semaev
template needs `4n` Boolean variables on the registered ambient
`subgroup_orbits` bases, which exceeds `MAX_VARS=64` on every cell, so all
twenty v2 F4/F5-family layouts are statically `unsupported` before solving.
That finding does not rewrite live receipts or decide SAT arms. The disclosed-point
[standard-subspace dimension-6 F4/F5 pilot](generic-f4-subspace-pilot/run-20260929/RESULT.md)
has executed once on the five inventory-control points (`max_trials=1`,
`node_budget=4096`). All ten `f4`/`f5` jobs dispatched into MatrixF4/MatrixF5
with `unsupported: false`, exhausted the node budget, and accepted zero
relations. It is a factor-base-policy diagnostic with `promotion_eligible=false`.
Do not rerun that budget. Do not register a fresh competitive F4/F5 panel from
it. Actions run 36580669479 is the sole empirical record for seed `2026092902`;
do not redispatch that seed.
