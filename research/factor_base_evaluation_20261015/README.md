# Controlled elliptic-curve factor-base comparison

This experiment compares multiplicative-coset, low-x, low-digit-weight and random
point factor bases on 48 admitted prime-order curves over 24 fields at 16, 20 and
24 bits. Each comparison holds the target curve fixed and matches the actual
number of canonical ± point classes.

The 24-bit exact-join comparison favored low-x in this cohort:

| Factor base | Independent rows per total CPU second |
| --- | ---: |
| Low-x | 1,352 |
| Coset | 1,167 |
| Random | 1,144 |
| Low digit weight | 859 |

Coset/random throughput was 1.020×, with a descriptive 95% whole-field bootstrap
interval of 0.959–1.083. This does not establish a coset throughput advantage.

A separate quadratic-algebra norm solver passed the 24-bit positive-case gate
(14/16 planted targets and 16/16 known-positive uniform targets). On fresh
24-bit inputs it produced 90 independent rows in 229.902 CPU seconds, versus
115 rows in 0.268 CPU seconds for the exact join. Setup, unsuccessful searches,
verification, rank work, supervisor CPU and measured budget overruns are charged.
The norm solver is a Sage prototype; the exact join uses native C. These are
bounded prototype measurements and do not establish NIST-size complexity.

## Evidence

- [Report (PDF)](results/report.pdf) and [HTML version](results/report.html).
- [Per-case measurements](results/cases.csv) and [aggregate CSV](results/generic_summary.csv).
- [Machine-readable summary](results/summary.json).
- [Mathematical validation receipt](results/validation.json): 26,497 independent
  point-relation identities and 402 matrix ranks checked in Sage, in addition to
  3,024 native point-operation checks and complete tiny-group solution sets.
- [Complete original run](artifacts/factor_base_evaluation.zip): source, exact
  curves and bases, targets, decompositions, independent rows, logs, original
  build artifact, discarded timing attempt and a SHA-256 manifest.

`protocol.json` records the fixed design. `gate_sampling_addendum.json` records
the separate conditional gate input stream. `timing_wrapper_addendum.json`
records the replacement of deferred SIGPROF timing with supervised workers.
The first timing attempt remains in the evidence archive and is excluded from
the final comparisons. Worker budgets have measured monitoring/scheduling
overruns; every actual overrun is recorded and charged.

## Historical snapshot and repository status

The complete experiment predates selecting this repository as its destination.
Its Sage/Python/C sources and original compiled artifact are preserved byte for
byte inside the archive. This PR imports historical evidence, not a new native
research implementation. Historical Sage/Python CPU timings remain exploratory
stage diagnostics. There is no full-DLP solve, matched rho reference, counted
end-to-end cost, or repository isolation/A/A receipt; those values remain unset.
The reference host was an Intel Xeon Platinum 8573C, Linux x86_64, pinned CPU 0.
Pinning alone does not satisfy this repository's isolation requirements.

Repository policy requires a native implementation and matched native reruns
before extending this experiment or using its timings as new performance claims.
No new experiment or statistical calculation was run while packaging this PR.
No new mathematical theorem is claimed; there is no Lean proof claim.

## Inspect and verify

Run `sh verify-artifacts.sh` with standard `unzip`, `sha256sum` and `cmp` tools.
It checks the original ZIP CRCs, all 32 archived manifest hashes, and the
byte-identical review copies. It does not rerun mathematical verification.
Use `unzip artifacts/factor_base_evaluation.zip -d reference-run` to inspect
source snapshots, raw outputs and original execution commands. Historical
execution commands in the archive are provenance, not a repository-native
reproduction workflow.

The receipt reports 26,497 saved relation identities and 402 matrix ranks,
3,024 native point-operation checks, and 905 complete tiny-group solution sets
for each of two- and three-point decompositions. These Sage checks were completed
before import; PR packaging repeats integrity checks only.

## Requirements and validation status

| Requested step | Status | Evidence |
| --- | --- | --- |
| 1. Admit prime-order curves | Verified in historical run | Archived curves and Sage receipt |
| 2. Compare four matched bases | Verified in historical run | Protocol, per-case CSV, raw bases and targets |
| 3. Check decomposition and fixed-column rank | Verified in historical run | Saved rows, rank receipt and norm-solver timeouts |
| Preserve reproducibility artifacts in a PR | Complete in this change | Full archive, review copies and hash checker |
| Native repository implementation and rerun | Not performed; required before extension | Historical Sage sources retained in archive |

`conductor check` and `cargo test --release --lib` were attempted during packaging;
both exited 127 because the commands are unavailable in this workspace.
The referenced `.conductor/WORKFLOW.md` returned 404 at the pinned base commit.
No library or example source is changed. Repository CI and review remain gates
before merge; unavailable local checks are not reported as passes.

## Visual records and canonical graphs

The original run/results report remains [PDF](results/report.pdf) and
[HTML](results/report.html). [Rendered figures](reports/render_figures.py) read
only the committed summary, producing SVG and 150 dpi PNG from the same source;
the original figures stay frozen under `results/`. Matching browsing copies live
in [benchmarks](../../benchmarks/20261015-factor-base-evaluation/).

Canonical graph applicability was checked against root AGENTS.md sections 7a–7c:
`docs/index-calculus-scoreboard.html`, `docs/ic/progress-timeline.json`,
`docs/ic/leaderboard.json` and `docs/browser/` remain unchanged because this import
adds no IC/rho ratio, whole-pipeline measurement, registered curve, ecbench
session or tournament identity. Stage-only CPU diagnostics must not populate
those end-to-end axes. No existing value or verdict is superseded.

## Construction reference

Petit, Kosters and Messeng, *Algebraic approaches for the Elliptic Curve Discrete
Logarithm Problem over prime fields*, PKC 2016, §§3.3–3.4:
<https://people.maths.ox.ac.uk/petit/files/16PKC_primeECDLP.pdf>.
