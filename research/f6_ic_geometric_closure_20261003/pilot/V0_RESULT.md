# F6-IC v0 pilot: geometric closure did not engage

This is the frozen, six-process, one-public-target comparison registered in
[`PROTOCOL.md`](../PROTOCOL.md) and [`PILOT_AMENDMENT.md`](../PILOT_AMENDMENT.md).
The source commit was `5eacbd2b42294a10c56c7360f258e7d45c71f1d3`, the worker binary SHA-256
was `d8d840a9feb98da45a51b62a9487dde47884a33e1026f8b021a8515352c1cc74`,
and the workload ID was `4de3483deb46`. The previously unseen public point
was `[86087,73652]` on the registered n17 Koblitz model; fixture generation
constructed no scalar. Both candidates used 62 actual usable points and 29
folded relation columns. The candidate records correct inherited F4's actual
`highest-free` split rule; the older immutable record's error is documented
in [`BASELINE_ERRATUM.md`](../BASELINE_ERRATUM.md).

| Repetition | Inherited F4 online ms | F6-IC v0 online ms | F4 / F6 reductions | F4 / F6 splits | F6 residual lookups |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 3.789 | 5.572 | 16 / 16 | 9 / 9 | 0 |
| 2 | 3.858 | 4.114 | 16 / 16 | 9 / 9 | 0 |
| 3 | 3.842 | 3.858 | 16 / 16 | 9 / 9 | 0 |

All six processes exited zero, reported `complete`, recovered scalar `25280`,
independently replayed it, and closed all five exclusive online phase costs
exactly to their reported online interval. Every solve used one target query
and one verified PDP witness. F6 checked 11 fully fixed summand coordinates
per run, but reached no residual closure or geometric refutation before F4
found the witness. The identical reductions and splits falsify v0's intended
mechanism for this workload. No 2x gain occurred.

The Mac host has no host-level CPU isolation receipt. These wall times are
**exploratory**; the controlled wall-time speedup and any IC-versus-rho
speedup remain **unknown**. This pilot made no cross-method run. Raw outputs,
stderr, exit codes, timestamps, source/input hashes, per-attempt stats and
five-phase ledgers are preserved in [`runs/`](runs/) and
[`v0-measurements.jsonl`](v0-measurements.jsonl). The latter has one row per
`(candidate_id, workload_id, run_id)`, including each raw-output SHA-256.

Decision: preserve v0 as a negative receipt. Version 1 will attempt exact
residual closure when **one** of three summand x-codes is fixed, while the
second summand ranges over the small exact base. This moves the group test
early enough to change the search tree, at a measured cost in point additions.
