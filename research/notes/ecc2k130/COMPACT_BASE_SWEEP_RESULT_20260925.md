# Compact four-sum sparse-base sweep at n=37 and n=41

**Decision:** The current materialized-root, exhaustive-S3 compact producer has no demonstrated end-to-end advantage over signed-Frobenius rho at these sparse bases. This is an occupancy and rank boundary diagnostic, not a shared-base crossover measurement. Root occupancy stays high, but sparse target relations require expensive complete scans on misses. Per-miss optimizations could reduce that cost; they cannot enlarge the fixed four-sum support ceiling. The n=131 transfer remains a no-go for this architecture under its stated materialized-root assumption, now with the sharper [unordered-count bound](../../sat_factor_base_review_20260908/autolab_orbit_extract_20260924/COMPACT_ORBIT_N131_UNORDERED_ADDENDUM_20260925.md). No n=131 attack, successful challenge log, exponent improvement, or full-DLP speedup is claimed.

## Frozen question and accounting

The [preregistered protocol](COMPACT_BASE_SWEEP_PROTOCOL_20260925.md) fixed two curves, three nonsaturated candidate bases per curve, a nested SHA-derived scalar stream, 64/512/4096 stage rules, three already fixed point-only Q holdouts per curve, and same-Q signed-Frobenius rho L=1 processes before producer measurement. The input specification SHA-256 is `700ab214498bdd18cc49530e0233ed2f135628124fdb3eea48cbda9e1ef6666b`. The instrumented compact source SHA-256 is `c2bc8b05087df69bef9593363e9d7c112e843ef16da122da50eb29ab22115f09`; the release binary SHA-256 is `c668b19ff89d429d8e47f06ba9696bb0d97f4fc92f58f717098687c5cabbb305`. Receipts retain the actual Git commit, command, environment, input and output hashes, process wall, child CPU, peak RSS, stderr, and even setup failures. All times below are cold, one thread, on the same host; rho arms are separate L=1 processes and their sum is not a batched-rho reference.

Independent replay reconstructs the signed-Frobenius factor-base orbits, validates each S3 witness and group lift, solves relation rank modulo the subgroup order, checks every solved representative log, and checks all scalar predictions after first full rank. The point-only holdouts and rho scalars are checked against the previously merged clean Q controls. The archive includes compact raw stdout, exact source snapshots, hashed inputs, failures, independent validation, and a fail-closed replay script.

## Training measurements

Every 512-query arm below finished with `nR²` regular pair states and zero exceptional states. `U/(2nR²)` is the unique-root occupancy, not the probability that a target has a four-sum relation. The exact last column is the universal *necessary* support ceiling `min(1,C(F+3,4)/(q-1))` for a uniform nonzero target; it is neither a fitted hit-rate model nor a prediction for this extractor. A 95% Wilson interval for hit fraction treats the frozen hash stream as an independent uniform sample; it is descriptive and does not make the deterministic stream random by assertion.

| n | R | F | hits / 512 | 95% Wilson interval | rank / R | unique roots / candidates | cold wall | peak RSS | exact support ceiling |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 37 | 1 | 74 | 4 | 0.0030–0.0199 | 1 / 1 | 71 / 74 | 2.148 s | 14.88 MB | 0.005868 |
| 37 | 2 | 148 | 45 | 0.0663–0.1156 | 2 / 2 | 288 / 296 | 8.208 s | 15.20 MB | 0.090248 |
| 37 | 3 | 222 | 175 | 0.3020–0.3839 | 3 / 3 | 651 / 666 | 15.532 s | 15.60 MB | 0.450829 |
| 41 | 4 | 328 | 1 | 0.0003–0.0110 | 1 / 4 | 1,288 / 1,312 | 50.289 s | 17.33 MB | 0.000893 |
| 41 | 8 | 656 | 9 | 0.0093–0.0331 | 8 / 8 | 5,168 / 5,248 | 206.061 s | 18.50 MB | 0.014164 |
| 41 | 12 | 984 | 40 | 0.0579–0.1046 | 12 / 12 | 11,640 / 11,808 | 496.248 s | 19.51 MB | 0.071490 |

In all six arms, each failed query exhausted exactly `2n²R²=F²/2` S3 calls: 2,738, 10,952, 24,642, 53,792, 215,168, and 484,128 in table order. The total 512-query S3 counts were 1,392,392; 5,198,767; 9,029,730; 27,492,228; 108,497,946; and 231,639,023. Base selection, regular-state scan, and index construction each took at most 11 ms inside the producer; the target loop dominated wall time. Thus this is a **failed-query cost** issue even where RSS is modest at toy size. Root collisions removed 1.4–4.1% of root candidates, so a memory-saving representation may help RSS but cannot alone supply missing four-sum coverage.

Across the six 512-query arms, the observed funnel was `2,784,781/4/4`, `10,397,513/45/45`, `18,059,375/175/175`, `54,984,456/1/1`, `216,995,888/9/9`, and `463,278,026/40/40` for partner roots / indexed partner hits / group-lift attempts. Every indexed hit lifted and produced a relation in these streams. These counts identify the target-conditioned partner-root/index step as the measured bottleneck, but do **not** show that all misses lie outside the factor-base four-sum support: an exact group-law oracle is needed to separate tuple-sum collisions and extractor incompleteness. In particular, n=37 R=3's 95% hit interval lies below its exact combinatorial ceiling; the ceiling itself can be loose because distinct multisets can share a sum.

The frozen stage rule extended only n=37 R=1 and n=41 R=4,8 to 4,096 training queries: their 512-query hit counts were below `R+4` and projected cold wall stayed below 1,800 s. Each extension was a new cold process with the same nested scalar prefix; the 512-query run was not resumed or charged as an incremental tail. All 64-query calibrations (hits in table order: 0, 7, 23, 0, 0, 5), all 512-query arms, and both early wrapper failures remain in the archive. The failures were a runner path error before producer launch and a validation-header error after a successful producer exit; neither is silently counted as a producer observation.

| n | R | hits / 4,096 | 95% Wilson interval | rank / R | total S3 calls | cold wall | peak RSS | exact support ceiling |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 37 | 1 | 14 | 0.00204–0.00573 | 1 / 1 | 11,180,684 | 17.072 s | 24.05 MB | 0.005868 |
| 41 | 4 | 2 | 0.00013–0.00178 | 2 / 4 | 220,232,615 | 404.120 s | 26.74 MB | 0.000893 |
| 41 | 8 | 56 | 0.01054–0.01771 | 8 / 8 | 870,910,598 | 1,630.414 s | 26.38 MB | 0.014164 |

The n=41 R=4 extension gained only one additional relation in the last 3,584 target queries and never attained rank four. n=41 R=8 attained full rank, but its roughly 871 million S3 calls were charged across successful and unsuccessful training queries. The 95% interval at n=41 R=8 overlaps its exact counting ceiling; this does not establish saturation or prove completeness.

## Fixed point-only targets and rho control

Each base faced the same three previously frozen Q points for its curve in a new point-only compact process. No holdout scalar label reached the compact producer. Across all six bases, **zero of three** point queries produced a relation, so even the full-rank training bases recovered **zero of three** holdout logs. Independent replay verified the failed-query S3 counts, point inputs and complete factor-base logs where rank existed. Three Q values per base are too few for a stable coverage estimate; the concrete full-DLP outcome on these fixed Q values is nevertheless zero completions.

| n | R | training rank | point-only relations/logs out of 3 | compact cold wall for 3 Q | same-Q rho L=1 logs out of 3 |
|---:|---:|---:|---:|---:|---:|
| 37 | 1 | 1 / 1 | 0 / 0 | 0.0617 s | 3 |
| 37 | 2 | 2 / 2 | 0 / 0 | 0.0602 s | 3 |
| 37 | 3 | 3 / 3 | 0 / 0 | 0.1742 s | 3 |
| 41 | 4 | 2 / 4 | 0 / 0 | 0.3506 s | 3 |
| 41 | 8 | 8 / 8 | 0 / 0 | 1.2712 s | 3 |
| 41 | 12 | 12 / 12 | 0 / 0 | 3.1489 s | 3 |

The independent cold rho L=1 process walls on the exact n=37 Q values were 0.2758, 0.0537 and 0.0611 s; at n=41 they were 0.4089, 0.2758 and 0.2853 s. All six recovered scalars satisfy `[d]G=Q`. These are **six separate single-log controls**. Adding them would not measure shared-work batched rho, and a ratio using only the failed compact point-query wall would exclude relation collection and base-log recovery. The completed-workload IC/rho speed ratio remains unset.

The [committed evidence archive](compact_base_sweep_20260925/evidence/) contains 27 raw run directories (six 64-query calibrations, six 512-query arms, three 4,096-query arms, six point-only holdouts, six rho controls), two wrapper-failure directories, exact source snapshots, input and stdout hashes, and independent reports. Its `SHA256SUMS` file has SHA-256 `d20e66aad75e2f68ed7143ce7df7989666db5d370ef8592f7c6646271cebf1c4`. Reproduce all hash, group-law, rank, point-log and rho checks with `python3 research/notes/ecc2k130/compact_base_sweep_20260925/replay_all.py`; gzip-compressed stdout is decoded by the verifier. The archive is 1.3 MB in the repository, so no external evidence location is needed.

## Next separately gated experiment

A fresh protocol should freeze a new, disjoint Q stream and run an exact group-law four-sum oracle on selected small bases before measuring another compact variant. Construct all pair sums of the signed base points with a hashed witness index; for each fresh Q, exhaustively test whether `Q - pair_sum` appears, including identity and repeated-point cases. Independently verify recovered four-point witnesses. Run the compact extractor on the **same point-only Q** and classify each oracle-positive/extractor-negative target. This separates base-support misses from S3/index misses, and permits a meaningful proposal for pair-index changes. Charge oracle construction, RSS and query wall. The existing sweep must not be retrofitted to this oracle after seeing its outcomes.

A five- or six-summand path is a separate solver/index architecture hypothesis. Exact unordered counting reaches a 1% *necessary* threshold at `F=60,591,280` or `4,121,293` respectively, versus `3,574,951,633` for four summands at n=131. No cost, rank, or memory conclusion transfers with those numbers. Any such variant needs its own preregistration, group-law correctness oracle, memory accounting, complete descent and linear algebra, and same-Q full-DLP comparison.
