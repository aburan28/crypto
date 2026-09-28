# n53 true shared-log batch: six verified same-Q pairs

The [pre-outcome protocol](PROTOCOL.md) froze one cold 220-column base-log training solve, three disjoint point-only Q blocks, L=32 prefixes, a resource-gated L=128 extension, and matched fresh-table signed-Frobenius Kuhn–Struik rho. The six producer pairs completed on one GitHub Actions host in [run 36122257965](https://github.com/aburan28/crypto/actions/runs/36122257965). Every training relation, orbit label, point relation, and recovered scalar passed the in-run independent group-law replay. The durable [119-file raw archive](evidence/archive_manifest.json) has SHA-256 `6100f9b47da4a5fb6fa3ee152f4ec7376368fe09749db1aa196e4df204d4dfeb`.

For the preregistered **per-block** comparison, even the incomplete IC child-process lower cost exceeds same-Q batch rho in all six pairs. The lower IC/rho wall ratio is 2.688–3.034 at L=32 and 1.541–1.730 at L=128; the conservative verified upper ratio is 4.321–4.867 and 2.472–2.788, respectively. This is a no-crossover observation for this implementation and these fixed streams. It is not a confidence interval, a common-operation-unit `S`, or a result on ECC2K-130.

## Frozen input, rank, and stopping gates

The clean checkout descends from main merge `fc27150df3238b6863ed5618c721e7fd8b6ce403`, with exact compact Rust source SHA-256 `c2bc8b05087df69bef9593363e9d7c112e843ef16da122da50eb29ab22115f09`. The point/validator manifest SHA-256 is `f1843670a169d65645bf83886aeffee2e60e25627faa0002763eb1e7c13d8362`; all six point-file digests are in the frozen protocol. The 220-orbit/23,320-point certified base has point-set hash `d859319015ea405fd18aee41b51396ce4edcab64ef66265d8edcdeb5e040eb71` and gzip SHA-256 `23397af2ef668aed0775bcb409e1ae19555357ded635452818c9a3812f679d08`. The panel records all other Rust/Python source, executable and per-input hashes, commands, environment, stderr, CPU, peak RSS, host and load.

The single cold training driver took **188.449 s** wall and 191.087 s user+system CPU; its Rust child took 124.513 s. Sequential training peak process-group RSS was 448,335,872 bytes (427.6 MiB). All 512/512 frozen training relations and 23,320 orbit-point labels passed independent replay; rank 220 first occurred at relation 461, and the next 51 labels were predicted from the solved base logs. The producer scanned 2,565,200 regular states, built 5,081,560 index entries, and made 23,689,479 S3 calls across the 512 targets. No training target failed.

All three L=32 pairs completed before the extension gate. Four times the slowest L32 IC point-process wall was 47.809 s, below the frozen 300 s threshold; its maximum point-process RSS was 415,072,256 bytes, below 1.75 GiB. L=128 then completed in all three blocks. There were zero failed, censored or replaced targets in the final run: 480/480 IC point relations and 480/480 rho scalars were replayed across six pairs (384 distinct Q values, since each L32 is a prefix of L128). Neither child received a scalar label.

## Cold matched wall and native work

Each IC lower cost is the cold Rust training child plus the cold point child, omitting necessary Python base-log solving/recovery. Each conservative upper cost is the complete training driver, cold point child, and Python point recovery/replay driver. Training is charged once **in each block comparison**. Rho is its cold batch child with a fresh table. Seconds include each child process's setup/index construction; the external independent rho audit is recorded separately. Order alternated IC–rho, rho–IC, IC–rho in each rung.

| L | Block | IC point child s | IC lower s | IC verified upper s | Rho child s | Lower / rho | Upper / rho |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 32 | 0 | 9.543 | 134.056 | 215.025 | 44.181 | 3.034 | 4.867 |
| 32 | 1 | 11.952 | 136.465 | 218.080 | 45.053 | 3.029 | 4.841 |
| 32 | 2 | 8.969 | 133.482 | 214.618 | 49.666 | 2.688 | 4.321 |
| 128 | 0 | 30.133 | 154.646 | 248.083 | 100.342 | 1.541 | 2.472 |
| 128 | 1 | 29.061 | 153.574 | 247.505 | 88.760 | 1.730 | 2.788 |
| 128 | 2 | 28.217 | 152.730 | 246.561 | 90.851 | 1.681 | 2.714 |

| L | Block | IC query S3 calls | Rho walk steps | Rho group additions | IC sequential peak MiB | Rho peak MiB |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 32 | 0 | 1,455,598 | 3,178,704 | 3,390,964 | 427.6 | 22.8 |
| 32 | 1 | 1,927,642 | 3,243,975 | 3,461,047 | 427.6 | 22.8 |
| 32 | 2 | 1,186,859 | 3,572,547 | 3,812,458 | 427.6 | 39.7 |
| 128 | 0 | 4,896,732 | 7,226,561 | 7,711,219 | 427.6 | 76.4 |
| 128 | 1 | 5,172,524 | 6,389,577 | 6,817,438 | 427.6 | 39.6 |
| 128 | 2 | 4,861,980 | 6,534,920 | 6,972,861 | 427.6 | 39.6 |

The exact `query_observations` in the archive also contain partner roots, indexed hits, group-lift attempts, per-query hit/miss and time; rho summaries contain table entries, cross-target solves, collisions, canonicalizations and other `charges`. A successful-partner `trials` count is **not** all S3 calls. S3 calls and rho group additions are different units; no measured conversion to a shared group-addition-equivalent `S` was frozen, so `S` and any attack-speed crossover remain unset.

## Shared training across three disjoint blocks

The secondary three-block portfolio charges the training solve once per L rung and sums three fresh IC point processes/recoveries against three fresh same-Q rho tables. At L32, the IC child-only lower/verified upper are 154.978/270.825 s versus rho 138.899 s (ratios 1.116/1.950). At L128, they are 211.925/365.250 s versus rho 279.952 s (ratios **0.757/1.305**). The L128 portfolio therefore spans rho under the preregistered incomplete lower and conservative verified upper bounds; it does not establish an algorithmic break-even. This accounting does not model reuse beyond the frozen three blocks or any L>128 regime. The verified-upper CPU totals are 281.455 s (L32) and 375.847 s (L128), versus rho 138.860 s and 279.894 s.

The certified base was frozen before this panel. Its earlier construction is not included in these training-plus-query bounds; charging it would only increase IC cost. Compilation and public-Q creation are outside both process timings. Host load is recorded for each arm, but three deterministic streams do not support a distributional speed claim. In particular, these n53 measurements do not transfer to n131 or the Certicom ECC2K-130 challenge, where the current materialized four-sum index already fails the separate coverage/memory admission bound.

## Retained failure and next decision

The [first run](https://github.com/aburan28/crypto/actions/runs/36120881174) stopped after both L32 block-0 children because an added checker incorrectly asserted `s3_calls >= partner_roots`; S3 may return two partner roots per call. Its unmodified [23-file archive](failed_attempt_1/archive_manifest.json), [original checker](failed_attempt_1/verify_pair_original.py), and [failure explanation](failed_attempt_1/FAILED_ATTEMPT_1.md) remain in Git. The corrected checker alone changed; frozen Q streams, producers, costs and gates did not. Focused CI rehashes and retroactively replays that attempt without counting it toward the six-pair timing decision, then separately replays every relation and scalar in the accepted 119-file archive.

The next bounded measurement should be **preregistered L=384** on the exact union of the three disjoint L128 Q blocks: one fresh IC training solve, one cold point-query/index process, and one matched cold rho process with a single shared table on those identical 384 Q values. The driver should separately time the mandatory base-log linear solve and scalar recovery from independent audit replay, keeping all stages and failures in the primary complete cost. This directly resolves the L128 portfolio interval without reusing three separate rho tables or extrapolating from it. The target order, seeds, source hashes, ceilings and stop rules must be frozen before that follow-up's outcomes. Separately, any n131 design must pass a lower-memory, low-exhaustive-miss supported-target-mass and useful-rank admission bound before solver tuning; a conductor change alone is not such a gate.
