# Current-source cold CPU for the four cells the instruction ledger left pending

The [preregistered protocol](PROTOCOL.md) was committed before this panel ran and before the Callgrind ledger was inspected. [Actions run 36694337420](https://github.com/aburan28/crypto/actions/runs/36694337420) is its one dispatch, on merged main `fe75eea69dc66660c3920342d83f17e3f4eb1151`. It materialized the blocked-prefilter twenty-file freeze, built the pinned W64 compact and 32-walk normal-basis rho v3 binaries, and ran `off_a`, `rho`, `off_b` on the frozen public Q. Do not dispatch the panel again.

Each arm is a fresh full process on one reserved Linux CPU with `RAYON_NUM_THREADS=1`. The charged cost is wait4 user+system CPU from startup through rank, linear algebra, every target recovery and the group checks. There are 20 blocks at n37/L1 and five blocks at each other cell. A second Linux host replayed every rank trace and scalar with `verify_cold.py --relocated`. The hosted receipts and that replay agree on status, A/A, isolation and the paired intervals. `method_crossover` stays null: the protocol requires both the instruction unit and this CPU unit to win, plus a separately frozen disjoint-Q confirmation, and neither extra condition is met here.

| n / L | Blocks | A/A median (95%) | Cold CPU / rho median (95%) | Eligible | Verified targets per arm |
| --- | ---: | --- | --- | --- | ---: |
| 37 / 1 | 20 | 0.9860 (0.9730–1.0078) | 3.7362 (3.6901–3.7803) | yes | 1 |
| 37 / 1,024 | 5 | 1.0022 (0.9912–1.0092) | 1.7064 (1.6933–1.7146) | yes | 1,024 |
| 41 / 1 | 5 | 0.9922 (0.9356–1.0480) | 17.0499 (16.0731–18.1962) | yes | 1 |
| 53 / 1 | 5 | 0.9983 (0.9950–1.0025) | 27.4483 (27.1951–27.6616) | yes | 1 |

All four A/A medians lie in [0.9, 1.1], all four A/A intervals contain 1, and every isolation log has zero contended samples. Every CPU interval lies above 1. Wall ratios are retained in the receipts and are not the protocol ratio. The n37 wall medians are about 1.0004 and 0.9998; the n37/L1024 wall interval is 0.8824–1.3035 and is too wide to carry a wall claim. The n41/L1 and n53/L1 wall medians are 2.9817 and 26.7823.

**Decision.** These four current-source cells have no matched-rho cold CPU crossover. Together with the already eligible n41/n53 L=1,024 rows (blocked/rho 1.3282 and 1.4727), every latest-source cold CPU cell on this Q set is above one. The batch Ir leads of 0.5165, 0.6221 and 0.7773 remain an instruction-unit observation. Class: accounting. The missing CPU column is filled; the algorithm and the instruction counts are unchanged. No calibrated group-addition S, disjoint-Q confirmation, n83 result or GF(2^131) transfer follows.

The raw cell tarballs, smoke archive and second-host receipts are under [evidence/run_36694337420](evidence/run_36694337420/MANIFEST.json). The manifest SHA-256 of `n37_L1.tar.gz` is `53f9b0e3413a375a5897025cf6b4f4c20f4318489e5b1b6ddf7f414fdcf2eeb2`.
