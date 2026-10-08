# F6-IC v3 pilot: batched point work, 2x gate missed

The preregistered v3 kernel batched two groups of independent packed point
additions per fixed-coordinate choice. Its source commit was
`f56aa9bd1851c4cd792ba7175d45c1b31f5a25a7`, worker binary SHA-256
`fe315508346da22d60f26a1f278d1fde5dea86f73e23f95072c2df8700e05165`,
and workload ID `dbffd5dbfc8c`. The fresh public target was
`[47109,42247]`, constructed without a scalar. Both candidate manifests,
inputs and exact hashes are in [`freeze.txt`](freeze.txt). No rho arm ran.

| Repetition | Inherited F4 online ms | F6-IC v3 online ms | F4 / F6 reductions | F4 / F6 splits | F6 packed lookups | F6 batch groups |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 87.461 | 58.102 | 696 / 295 | 292 / 160 | 16,801 | 534 |
| 2 | 88.125 | 58.892 | 696 / 295 | 292 / 160 | 16,801 | 534 |
| 3 | 89.361 | 59.037 | 696 / 295 | 292 / 160 | 16,801 | 534 |

All six processes exited zero, retained all five attempts (four exact
no-relation outcomes and one group-verified witness), recovered scalar
`41891`, and independently replayed it. The five exclusive online phase
costs sum exactly to each online interval. F6 counted 33,645 logical point
additions per call, including general witness replay; the batch kernel
amortised their inversions but did not omit the work. Exact n9 controls
compared batched, scalar packed and general branch decisions for both
partial-bit values and positive/negative cases.

The complete online F6 calls were descriptively about two thirds of F4's
time on this target, short of the preregistered 2x engineering target. These
macOS CPU times are **exploratory** under the host-isolation rule. The
controlled wall-time speedup, IC-versus-rho ratio, and any end-to-end attack
improvement remain **unknown**. Raw attempts, phase ledgers, stderr, exit
codes and timestamps are in [`runs/`](runs/) and the keyed
[`v3-measurements.jsonl`](v3-measurements.jsonl).

Decision: stop this branch of F6-IC. The exact geometric oracle is a useful
IC-specific algorithm experiment, but in the prepared n17 pilot its search
reduction did not yield a 2x complete-call gain. It is not promoted as an
F4/F5 replacement, and no larger-field or asymptotic claim follows. The
next optimization work returns to the shared F4/F5 arithmetic and matrix
paths, where the users' original profiles showed larger exclusive costs.
