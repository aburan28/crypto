# F6-IC v2 pilot: packed geometry helps, 2x target unmet

The preregistered v2 path used the exact single-word `FastCurve` for the
one-fixed-summand residual enumeration and general group replay for its
witness. Source commit was `9d9f5de881f3f548681fd00a84e5b9e16b6cad04`,
binary SHA-256 was `b9046be2973771cebd33a2efe55552f3f13ad4c890daccbd4597ed0c19bea909`,
and the fresh one-target workload ID was `63154d0a34c9`. The public point was
`[73951,104451]`; the fixture constructed no scalar. Candidate and input
hashes are in [`freeze.txt`](freeze.txt). This was an IC-variant pilot only.

| Repetition | Inherited F4 online ms | F6-IC v2 online ms | F4 / F6 reductions | F4 / F6 splits | F6 packed lookups | F6 point additions |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 119.205 | 82.581 | 911 / 383 | 384 / 209 | 22,223 | 44,449 |
| 2 | 118.584 | 83.067 | 911 / 383 | 384 / 209 | 22,223 | 44,449 |
| 3 | 118.850 | 82.966 | 911 / 383 | 384 / 209 | 22,223 | 44,449 |

All six processes exited zero and recovered scalar `28062` with independent
replay. They retained all six target-query attempts: five exact no-relation
outcomes, then a verified witness. All five exclusive phase costs summed to
each online interval. All 22,223 F6 residual lookups used the packed path.
Total F4 reductions fell by 58%, but complete online time did not reach the
preregistered exploratory 2x target. The fast point kernel improved the
engineering cost; it did not change the exact search decisions.

The Mac host is unisolated, so these wall times are **exploratory** and the
controlled wall-time speedup remains **unknown**. No new rho run or IC/rho
claim was made. [`runs/`](runs/) and [`v2-measurements.jsonl`](v2-measurements.jsonl)
retain every raw output, exit code, attempt, phase ledger and output digest.

Decision: preserve v2. The remaining 44,449 packed additions each pay for
an inversion. The next bounded variant batches the two groups of independent
additions per fixed-coordinate branch, using the existing exact batch
inversion kernel, and charges the batch operations in the same PDP interval.
