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
No reference selection is a promotion. Zero of the three improvement rounds have
run. The [bounded protocol](improvement/PROTOCOL.md) now fixes the nominal
familywise rule, 5,133-point historical exclusions and the first 16-pipeline
registry. The implementation PR must pass candidate controls before dispatch.
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
