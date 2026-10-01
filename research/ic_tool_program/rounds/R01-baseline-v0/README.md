# R01: baseline v0

The `ic` tool programme's first round
(`research/notes/index-calculus/IC_TOOL_PROGRAM.md` §12). It freezes
suite v1, pins v0 against ledger §23, profiles every phase, and measures
the host's noise floor. It has no candidate and claims no speedup.

**Status: declared, not yet run.** The results, the analysis and v0's
ledger row are added with the runs.

| file | what it is |
|:--|:--|
| `PROTOCOL.md` | the declaration: v0, the frozen inputs, the steps in order, the figures, the stop rules |
| `run.py` | the steps, resumable, never overwriting (`IC=<v0> IC_COMMIT=<commit> python3 run.py all`) |
| `analyse.py` | every figure (`python3 analyse.py > analysis.json`) |

The suite is `../../suite/v1/`, and the shared runner and statistics are
in `../../harness/`.
