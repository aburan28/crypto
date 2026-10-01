# The `ic` tool programme: rounds and baselines

The plan, its rules and its backlogs are in
[`research/notes/index-calculus/IC_TOOL_PROGRAM.md`](../notes/index-calculus/IC_TOOL_PROGRAM.md).
This directory holds what the plan freezes and what each round measures.
It adds each part with the round that creates it:

| path | what it holds | added by |
|:--|:--|:--|
| `suite/v1/` | the frozen suite: parameter files, `SUITE.json` with every file's SHA-256, and the script that wrote them | R01 |
| `harness/` | the shared runner (ABAB order, isolation, resumable, never overwrites) and the analysis | R01 |
| `fuzz/` | B6's generator of v2 documents, pinned by `SHA256SUMS` | B6 |
| `harness/bround.py` | Track B's steps for every B round: conformance on both arms, the 90-row pin, the v1 to v2 translation check, and the ABAB timing check against R01's A/A bands | Track B, before B0 runs |
| `rounds/R<k>-<slug>/` | one round: `PROTOCOL.md` committed before its candidate code, then its runs (as `runs.tar.xz` with its SHA-256), analysis, decision and README | each round |
| `baselines.json` | the ledger's rows as data: baseline, commit, binary hash, host and per-size figures | R01, then every accepted round |
| `conformance/v1/` | the conformance (C) suite: one case per defect or refusal, and its runner | B0 |
| `design/` | Track B's designs: the schema, the checks, the routing, and the cases each B step is judged against ([`schema-v2.md`](design/schema-v2.md)); the F1 sampled level ([`f1-sampled.md`](design/f1-sampled.md)) | B1, B7 |
| `conformance/v2/` | B1's cases and their parameter files, the script that writes and checks them, and the runner (v1's cases, then v2's) | B1 |
| `conformance/v2-b2/` | B2's cases (C032–C051) beside B1's frozen ones | B2 |
| `conformance/v2-b2b/` | B2b's cases (C059–C070) beside the earlier steps' frozen ones | B2b |
| `conformance/v2-b7a/` | B7a's cases (C071–C077), each an earlier step's frozen document at `fidelity: F1` | B7a |
| `conformance/v2-b3/` | B3's cases (C052–C058) beside B1's frozen ones | B3 |
| `conformance/v2-b3b/` | B3b's cases (C078–C087) and measurement 5's public targets | B3b |
| `conformance/run.py` | the runner from B3 on: every step's cases in order, with the `until` rule applied, and from B3b the `supersedes` rule | B3, B3b |
| `track-b/` | Track B's implementations on record before their measurements: a git bundle of every local branch, B0 to B3b, with its SHA-256 and how to restore it | B3b's amendments |

## Status

| round | what | state | PR |
|:--|:--|:--|:--|
| plan | goals, suite, measurement, loop, backlogs | merged | #1103 |
| R01 | baseline v0: freeze suite v1, pin against §23, profile every phase, A/A | complete ([`rounds/R01-baseline-v0/`](rounds/R01-baseline-v0/README.md)) | #1104 (declaration), #1115 (results) |
| R02 | the AVX-512 batched addition for fields with `n + deg t = 66` | **rejected** ([`rounds/R02-wide-tail-kernel/`](rounds/R02-wide-tail-kernel/README.md)): collection 1.29–1.31× and cold 1.18–1.27× faster at the two wide sizes, but the holdouts at `icv1-f2m59-tm943548413-98844ecc` read 1.161 [1.089, 1.239], whose lower end is not above 1.10; the kernel is kept as `candidate.patch` | #1117 (declaration), #1152 (results) |
| R02b | R02's wide-tail kernel re-tested: fresh holdouts with forty pairs a size, the power fixed before running, R02's callgrind control restated by kernel; the kernel's last test | declared ([`rounds/R02b-wide-tail-retest/PROTOCOL.md`](rounds/R02b-wide-tail-retest/PROTOCOL.md)); runs after R03's decision, before B0 | #1157 (declaration) |
| R03 | the curve's construction at composite degrees: `factorise_u64` tests primality only when a division changes what is left | declared ([`rounds/R03-curve-construction/PROTOCOL.md`](rounds/R03-curve-construction/PROTOCOL.md)); runs after R02's decision and R04 | #1119 (declaration) |
| B0 | refusals and fixes for the survey's defects, with conformance suite v1 | declared ([`rounds/B0-refusals/PROTOCOL.md`](rounds/B0-refusals/PROTOCOL.md), [`conformance/v1/`](conformance/v1/cases.json)); runs after R03 | #1119 (declaration) |
| B1 | schema v2, the checks for binary instances, and the Koblitz pipelines on imported instances | declared ([`rounds/B1-schema-v2/PROTOCOL.md`](rounds/B1-schema-v2/PROTOCOL.md); design [`design/schema-v2.md`](design/schema-v2.md); cases [`conformance/v2/`](conformance/v2/cases.json)); runs after B0 | #1125 (declaration) |
| B2 | prime and extension fields, the other curve forms and importers, estimates and budgets | declared ([`rounds/B2-fields-forms-estimates/PROTOCOL.md`](rounds/B2-fields-forms-estimates/PROTOCOL.md); cases [`conformance/v2-b2/`](conformance/v2-b2/cases.json)), amended before measuring (amendment 1: B2b split off, the estimates' models and [`estimates.json`](rounds/B2-fields-forms-estimates/estimates.json)); runs after B1, and after B3 if B3 is accepted first | #1140 (declaration), #1141 (amendment 1) |
| B2b | `kic` on subfield curves with `k > 1` (paired with `rho-negation`), `kic` alone, order certificates | declared ([`rounds/B2b-subfield-kic-certificates/PROTOCOL.md`](rounds/B2b-subfield-kic-certificates/PROTOCOL.md); cases [`conformance/v2-b2b/`](conformance/v2-b2b/cases.json); the subfield sweep [`sweep.py`](rounds/B2b-subfield-kic-certificates/sweep.py)), amended before measuring (amendment 1: v2 target rules, v1's generator rule and translation for `k > 1`, the report keys); runs after B2 | #1144 (declaration), #1145 (amendment 1) |
| R04 | where a scanned summand's time goes: per-stage counters inside the scan, in a probe build (a stage diagnostic) | complete ([`rounds/R04-scan-probes/`](rounds/R04-scan-probes/README.md)), a stage diagnostic of class accounting: the batched subtraction leads the scan at three of the four largest sizes (37–41%), the admitted keys at `icv1-f2m53-tm56619371-dac20a85` (42%); the probes' overhead interval lies inside R01's A/A band at no size, so the shares are reported, not trusted; a check after the run finds nearly every admitted key a false positive of the presence filter | #1128 (declaration), #1146 (amendment 1), #1150 (amendment 2), #1158 (results) |
| B3 | two-word binary fields (`n ≤ 126`), `rho-koblitz` on them, `solve: rho`, and the gate's rho at F0 at `n = 83` | declared ([`rounds/B3-two-word-rho/PROTOCOL.md`](rounds/B3-two-word-rho/PROTOCOL.md); cases [`conformance/v2-b3/`](conformance/v2-b3/cases.json)); runs after B1 | #1139 (declaration) |
| B3b | the index calculus on two-word fields (`63 ≤ n ≤ 126`) at F0: two-word kernels for every phase, a second pipeline over the same algorithm, the same as the one-word pipeline counter for counter where both run | declared ([`rounds/B3b-two-word-kic/PROTOCOL.md`](rounds/B3b-two-word-kic/PROTOCOL.md); design [`design/two-word-kic.md`](design/two-word-kic.md); cases [`conformance/v2-b3b/`](conformance/v2-b3b/cases.json)), amended before any measurement (amendment 1: `r ≤ h` is refused at one word only; amendment 2: C086 and C087 supersede C031 and C053, whose `kic` width gate B3b moves, and C082's document name is shortened to the schema's limit); implemented, on record in [`track-b/`](track-b/README.md); runs after B3 | #1151 (declaration), #1161 (amendments 1–2) |
| B6 | fuzzing and differential checks: a seeded generator of valid and corrupted documents, every answer replayed in independent Python arithmetic, in CI and as a 20,000-document campaign | declared ([`rounds/B6-fuzzing/PROTOCOL.md`](rounds/B6-fuzzing/PROTOCOL.md); generator [`fuzz/fuzz_v2.py`](fuzz/fuzz_v2.py)); runs after B2 | #1143 (declaration) |
| B7a | the F1 sampled level at one-word sizes: `kic` and every rho extrapolated from samples on the instance, the yield counted, with predictive intervals; validated against F0 | declared ([`rounds/B7a-f1-sampled/PROTOCOL.md`](rounds/B7a-f1-sampled/PROTOCOL.md); design [`design/f1-sampled.md`](design/f1-sampled.md); cases [`conformance/v2-b7a/`](conformance/v2-b7a/cases.json); carried constants [`carried.json`](rounds/B7a-f1-sampled/carried.json)); amendment 1, before any run: the count from distinct keys, the trial granule, the synthetic linear algebra solved in full, the first-repetition context `χ`, partial tables; runs after B2b | #1148 (declaration), #1149 (amendment 1) |
