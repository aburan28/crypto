# The `ic` tool programme: rounds and baselines

The plan, its rules and its backlogs are in
[`research/notes/index-calculus/IC_TOOL_PROGRAM.md`](../notes/index-calculus/IC_TOOL_PROGRAM.md).
This directory holds what the plan freezes and what each round measures.
It adds each part with the round that creates it:

| path | what it holds | added by |
|:--|:--|:--|
| `suite/v1/` | the frozen suite: parameter files, `SUITE.json` with every file's SHA-256, and the script that wrote them | R01 |
| `harness/` | the shared runner (ABAB order, isolation, resumable, never overwrites) and the analysis | R01 |
| `rounds/R<k>-<slug>/` | one round: `PROTOCOL.md` committed before its candidate code, then its runs (as `runs.tar.xz` with its SHA-256), analysis, decision and README | each round |
| `baselines.json` | the ledger's rows as data: baseline, commit, binary hash, host and per-size figures | R01, then every accepted round |
| `conformance/v1/` | the conformance (C) suite: one case per defect or refusal, and its runner | B0 |
| `design/` | Track B's designs: the schema, the checks, the routing, and the cases each B step is judged against | B1 |
| `conformance/v2/` | B1's cases and their parameter files, the script that writes and checks them, and the runner (v1's cases, then v2's) | B1 |
| `conformance/v2-b2/` | B2's cases (C032–C051) beside B1's frozen ones | B2 |
| `conformance/v2-b3/` | B3's cases (C052–C058) beside B1's frozen ones | B3 |
| `conformance/run.py` | the runner from B3 on: every step's cases in order, with the `until` rule applied | B3 |

## Status

| round | what | state | PR |
|:--|:--|:--|:--|
| plan | goals, suite, measurement, loop, backlogs | merged | #1103 |
| R01 | baseline v0: freeze suite v1, pin against §23, profile every phase, A/A | complete ([`rounds/R01-baseline-v0/`](rounds/R01-baseline-v0/README.md)) | #1104 (declaration), #1115 (results) |
| R02 | the AVX-512 batched addition for fields with `n + deg t = 66` | declared with amendments 1–3 ([`rounds/R02-wide-tail-kernel/PROTOCOL.md`](rounds/R02-wide-tail-kernel/PROTOCOL.md)); running | #1117 (declaration) |
| R03 | the curve's construction at composite degrees: `factorise_u64` tests primality only when a division changes what is left | declared ([`rounds/R03-curve-construction/PROTOCOL.md`](rounds/R03-curve-construction/PROTOCOL.md)); runs after R02's decision | #1119 (declaration) |
| B0 | refusals and fixes for the survey's defects, with conformance suite v1 | declared ([`rounds/B0-refusals/PROTOCOL.md`](rounds/B0-refusals/PROTOCOL.md), [`conformance/v1/`](conformance/v1/cases.json)); runs after R03 | #1119 (declaration) |
| B1 | schema v2, the checks for binary instances, and the Koblitz pipelines on imported instances | declared ([`rounds/B1-schema-v2/PROTOCOL.md`](rounds/B1-schema-v2/PROTOCOL.md); design [`design/schema-v2.md`](design/schema-v2.md); cases [`conformance/v2/`](conformance/v2/cases.json)); runs after B0 | #1125 (declaration) |
| B2 | prime and extension fields, the other curve forms and importers, estimates and budgets | declared ([`rounds/B2-fields-forms-estimates/PROTOCOL.md`](rounds/B2-fields-forms-estimates/PROTOCOL.md); cases [`conformance/v2-b2/`](conformance/v2-b2/cases.json)); runs after B1, and after B3 if B3 is accepted first | this PR (declaration) |
| R04 | where a scanned summand's time goes: per-stage counters inside the scan, in a probe build (a stage diagnostic) | declared ([`rounds/R04-scan-probes/PROTOCOL.md`](rounds/R04-scan-probes/PROTOCOL.md)); runs after B1 | #1128 (declaration) |
| B3 | two-word binary fields (`n ≤ 126`), `rho-koblitz` on them, `solve: rho`, and the gate's rho at F0 at `n = 83` | declared ([`rounds/B3-two-word-rho/PROTOCOL.md`](rounds/B3-two-word-rho/PROTOCOL.md); cases [`conformance/v2-b3/`](conformance/v2-b3/cases.json)); runs after B1 | #1139 (declaration) |
| B3b | the index calculus on two-word fields, at F1 | to be declared | — |
