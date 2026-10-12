# Held-out K sweep: no eager-index timing crossover against strong rho

The preregistered n41/n53 L=1,024 K grid completed on Linux and the two downloaded raw artifacts passed independent macOS replay. Each arm saw exactly the same **point-only** Q file in each of five balanced cold blocks. Both files contain 1,024 distinct Q and have zero overlap with the prior point-panel corpora. Every one of the 40 compact processes reached full rank and recovered all 1,024 logs; the 10 normal-basis rho processes recovered the same logs. The independent verifier replayed all 40 rank transitions, base logs, every four-sum and every recovered scalar, plus every rho answer. That is 51,200 verified arm-target outputs. There were no timeouts, failed targets, rank misses or OOMs. Five repeats measure runner variation on the same Q, not target-distribution uncertainty.

The four K values were frozen before these Q were generated: n41 K=128/192/255/320 and n53 K=220/330/440/550. K=255/440 are the prior selected controls. Every compact child constructs its own signed-orbit factor base and complete eager S3 index, reaches rank, solves all targets and verifies them; the matched reference is PR #943's 32-walk batched-inversion normal-basis signed-Frobenius rho on the same host and Q. Scalar labels stayed in the verifier's separate fixture file. Source/input hashes, order, limits and acceptance rule are in [the freeze](FROZEN.json) and [preregistered protocol](PROTOCOL.md).

| n | Arm | Cold CPU median, s | Wall median, s | RSS median, MiB | Paired CPU / rho | 95% paired CPU interval | Verified | S=operations/√r |
|---:|---|---:|---:|---:|---:|---|---:|---|
| 41 | K=128 | 4.4844 | 4.4899 | 87.8 | 3.491 | 3.461–3.509 | 5,120/5,120 | unset |
| 41 | K=192 | 3.1510 | 3.1551 | 173.5 | 2.443 | 2.422–2.478 | 5,120/5,120 | unset |
| 41 | K=255 | 3.2144 | 3.2191 | 329.4 | 2.493 | 2.462–2.530 | 5,120/5,120 | unset |
| 41 | K=320 | 3.8953 | 3.8992 | 621.4 | 3.029 | 3.009–3.049 | 5,120/5,120 | unset |
| 41 | normal-basis rho | 1.2891 | 1.2923 | 39.6 | 1.000 | reference | 5,120/5,120 | unset |
| 53 | K=220 | 34.9846 | 34.9917 | 327.9 | 5.553 | 5.527–5.574 | 5,120/5,120 | unset |
| 53 | K=330 | 19.4660 | 19.4730 | 661.7 | 3.087 | 3.077–3.102 | 5,120/5,120 | unset |
| 53 | K=440 | 16.4827 | 16.4891 | 1281.3 | 2.616 | 2.565–2.727 | 5,120/5,120 | unset |
| 53 | K=550 | 17.4305 | 17.4376 | 1420.0 | 2.766 | 2.749–2.786 | 5,120/5,120 | unset |
| 53 | normal-basis rho | 6.3025 | 6.3057 | 150.2 | 1.000 | reference | 5,120/5,120 | unset |

All eight compact/rho paired CPU intervals exclude parity on the losing side. The n41 grid's lowest complete CPU cost is K=192: 3.1510 s versus 3.2144 s at prior K=255. Their **direct paired** CPU ratio is 0.979 (95% interval 0.974–0.990), a small exploratory K-policy improvement with roughly half the RSS; selecting it from this evaluation grid requires a new Q holdout for a confirmatory claim. Even at K=192, compact costs **2.443×** the matched rho (2.422–2.478). At n53, prior K=440 remains the grid minimum at 16.4827 s and **2.616×** rho (2.565–2.727). The n41 and n53 cells ran on different EPYC models, so only within-cell paired ratios are interpreted; cross-n exponents are not fitted. Rho recovered all 1,024 Q in each block, with median 3,932,363/19,617,899 walk steps and 131,052/653,871 batch inversions at n41/n53.

The measured phase tradeoff explains the loss. Timings below are medians of each compact child's in-process wall phases, while probe columns are rank-probe means and target-probe means from the retained traces. They are not a calibrated common operation unit.

| n | K | S3 states | Index build, s | Rank, s | Targets, s | Rank probes/attempt | Target probes/Q |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 41 | 128 | 671,744 | 0.410 | 0.476 | 3.560 | 10,247 | 10,235 |
| 41 | 192 | 1,511,424 | 0.941 | 0.370 | 1.796 | 4,084 | 4,737 |
| 41 | 255 | 2,666,025 | 1.678 | 0.408 | 1.079 | 2,481 | 2,534 |
| 41 | 320 | 4,198,400 | 2.620 | 0.459 | 0.767 | 1,494 | 1,594 |
| 53 | 220 | 2,565,200 | 1.566 | 6.071 | 27.243 | 85,007 | 82,865 |
| 53 | 330 | 5,771,700 | 3.534 | 4.111 | 11.710 | 36,440 | 35,111 |
| 53 | 440 | 10,260,800 | 6.287 | 3.282 | 6.779 | 19,946 | 20,130 |
| 53 | 550 | 16,032,500 | 9.863 | 2.984 | 4.464 | 12,636 | 12,872 |

For all eight K choices, the eager index materialized exactly **K²·n** regular states, and each root table retained almost as many distinct canonical entries. Reducing n53 K from 440 to 220 cut index construction from 6.287 to 1.566 s, but raised rank from 3.282 to 6.071 s and 1,024-target work from 6.779 to 27.243 s; mean target probes rose from 20,130 to 82,865 per Q. Increasing K to 550 cut target work to 4.464 s but took 9.863 s to build the index. The rank stage reached K independent rows in exactly K attempts at every K: the **attempt-count** ratio to the K+L floor is 1.000, while the probe burden remains large. This is a construction-versus-query cost tradeoff, not a relation-yield shortfall.

The closest n41 arm, K=192, spends 1.796 s in its target phase alone against rho's entire 1.292 s wall median; at n53 K=440, the target phase is 6.779 s against rho's 6.306 s. Conversely, K=255/550 reduce target cost but build a larger index. Improving only setup or only rank cannot cross at the measured grid minima. The next candidate must jointly reduce S3 index work and relation-query cost at equal useful workload. One concrete untested lead is to quotient the ordered `(left,right,relative)` state enumeration by its swap/Frobenius involution before building the index; this requires a proof that all canonical roots and four-sum witnesses remain reachable, an exceptional-case test, and a new held-out full-log comparison. The present panel does **not** measure that policy or any degree-263 descendant.

The [sealed archive](evidence/k_grid_0d6d/MANIFEST.json) holds both complete raw tarballs, process streams, base/rank/target traces, Linux receipts, second-machine receipts, host/binary/source hashes, and SHA-256 checksums. All discrete replay results agree exactly; the two endpoints of one n53 derived wall interval differ only at floating-point roundoff between Linux and macOS, within the archive's 10⁻¹² equivalence check. The passing [GitHub run](https://github.com/aburan28/crypto/actions/runs/36570106399) compiled PR head `0d6db9d48b2ba8353ddeaaa3100417463a37cdda` in synthetic checkout `32063e0f673876ca86d421c62ba3afaebb19972e`, whose main parent was `0adc50ecd67b868eced5e0de713d70bfd317b9fc`. The archived Q and source SHA-256 values are those frozen before this run.

**Decision.** This is a quantitative no-go for the tested **eager-index K grid on these held-out n41/n53 batches** against the strongest matched rho. It does not prove a no-go for untested K, a new index or PDP oracle, oriented degree-263 descendants, n83, or ECC2K-130. Common fully charged `S`, an operation-counted attack-speed ratio, generic-group-floor comparison and n131 transfer assumptions remain unset. The next evidence-ranked step is a separate, frozen joint index/query policy test, plus a common calibrated operation ledger for both IC and rho. The descendant-native/transported/pullback comparison still needs equal-useful-size n≥3 PDP workloads on held-out points; the earlier n131 m=2 zero-hit smoke was too sparse to decide that question.
