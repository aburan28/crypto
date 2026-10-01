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

## Status

| round | what | state | PR |
|:--|:--|:--|:--|
| plan | goals, suite, measurement, loop, backlogs | merged | #1103 |
| R01 | baseline v0: freeze suite v1, pin against §23, profile every phase, A/A | complete ([`rounds/R01-baseline-v0/`](rounds/R01-baseline-v0/README.md)) | #1104 (declaration), #1115 (results) |
| R02 | the AVX-512 batched addition for fields with `n + deg t = 66` | declared with amendments 1–3 ([`rounds/R02-wide-tail-kernel/PROTOCOL.md`](rounds/R02-wide-tail-kernel/PROTOCOL.md)); running | #1117 (declaration) |
| R03 | the curve's construction at composite degrees: `factorise_u64` tests primality only when a division changes what is left | declared ([`rounds/R03-curve-construction/PROTOCOL.md`](rounds/R03-curve-construction/PROTOCOL.md)); runs after R02's decision | this PR (declaration) |
| B0 | refusals and fixes for the survey's defects, with conformance suite v1 | declared ([`rounds/B0-refusals/PROTOCOL.md`](rounds/B0-refusals/PROTOCOL.md), [`conformance/v1/`](conformance/v1/cases.json)); runs after R03 | this PR (declaration) |
