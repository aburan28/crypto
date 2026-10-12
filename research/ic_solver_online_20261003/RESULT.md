# First public-target development panel

This panel used the registered point `[36389,93850]` (seed `20261003017`),
the exact frozen worker binary, prepared 62-point usable base and 29-column
certified log table, and three fixed process orders. All nine worker reports
said `complete` and independently replayed the recovered scalar `31760`.
Each IC solve needed two target attempts: first `proved_unsat`, then a witness.
The five exclusive IC online phases sum exactly to each reported online total.

| Repetition | F5 online ms | Inherited F4 online ms | Rho online ms | F5 / F4 | Rho / F4 |
| --- | ---: | ---: | ---: | ---: | ---: |
| R1 | 10154.197 | 45.360 | 0.629 | 223.86 | 0.01386 |
| R2 | 8757.058 | 42.003 | 0.097 | 208.49 | 0.00232 |
| R3 | 8700.804 | 47.995 | 0.153 | 181.29 | 0.00318 |

Median paired F5 / inherited-F4 online ratio: **208.49×** (range
181.29–223.86×). Median paired rho / inherited-F4 ratio: **0.00318×**;
inherited F4 remained about 315× slower than this rho reference on the
median paired ratio. F5's and inherited F4's online intervals were almost
entirely `target_pdp` (over 99.9% for both). F5 used 2,546 splits in the
first unsuccessful target PDP attempt and 168 in the successful attempt;
inherited F4 used 68 and 13, respectively. These counts explain why a
standalone matrix-kernel gain cannot be projected onto this path.

**Capture defect:** The registered shell wrapper ran `/usr/bin/time -l` around
each worker. On this sandboxed macOS host, `time` emitted the worker's timing
then failed its own `sysctl kern.clockrate` read with exit status 1. Thus
each process receipt records wrapper exit 1 even though worker stdout reports
complete, with verified scalar replay and complete online ledgers. The frozen
gate requires a clean process exit, so this panel is **diagnostic, not a
passing end-to-end speedup claim**. The raw stdout, stderr, dates, exit
statuses, host and compiler records are retained under `runs/`. No run was
discarded, replaced or rerun within this registered panel.

The n17 result also has a scaling limit: rho solved this point in tens of
walk additions, whereas this IC path still performs a full Boolean PDP for
the target. This result neither establishes an ECC2K-130 advantage nor a
general F5-versus-F4 ranking. A separate fresh-target confirmation is
registered in `confirmation/PROTOCOL.md` before its execution.
