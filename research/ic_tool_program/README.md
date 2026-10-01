# The `ic` tool programme: rounds and baselines

The plan, its rules and its backlogs are in
[`research/notes/index-calculus/IC_TOOL_PROGRAM.md`](../notes/index-calculus/IC_TOOL_PROGRAM.md).
This directory holds what the plan freezes and what each round measures.
It adds each part with the round that creates it:

| path | what it holds | added by |
|:--|:--|:--|
| `suite/v1/` | the frozen suite: parameter files, `SUITE.json` with every file's SHA-256, and the script that wrote them | R01 |
| `harness/` | the shared runner (ABAB order, isolation, resumable, never overwrites) and the analysis | R01 |
| `rounds/R<k>-<slug>/` | one round: `PROTOCOL.md` committed before its candidate code, then its runs, analysis, decision and README | each round |
| `baselines.json` | the ledger's rows as data: baseline, commit, binary hash, host and per-size figures | R01, then every accepted round |

## Status

| round | what | state | PR |
|:--|:--|:--|:--|
| plan | goals, suite, measurement, loop, backlogs | merged | #1103 |
| R01 | baseline v0: freeze suite v1, pin against §23, profile every phase, A/A | complete ([`rounds/R01-baseline-v0/`](rounds/R01-baseline-v0/README.md)) | #1104 (declaration), results PR |
| R02 | the AVX-512 batched addition for fields with `n + deg t = 66` | declared, next | — |
| B0 | refusals and fixes for the survey's defects, with conformance suite v1 | written, after R02 | — |
