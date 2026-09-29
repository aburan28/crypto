# Prime F4 lazy AVX2 experiment

[CI run 36458164446](https://github.com/aburan28/crypto/actions/runs/36458164446) completed the frozen one-thread experiment on an AMD EPYC 7763 with AVX2. The AVX2 boundary test and prime F4 tests passed. All 132 benchmark processes across four seed workloads returned successfully, and every wide, current and candidate call matched its workload's basis and solve fingerprints, step count, matrix shape and basis size. The runner pinned one CPU and used one Rayon thread. Both release binaries remain in the CI artifact; their SHA-256 digests, source digests, CPU features, load, full output and per-call status are in the [compressed raw receipt](runs/36458164446-t1.json.gz). The compressed file's SHA-256 is `a0059f6e0460ec1607b8d9dd888b6492f563b2d3d08d02bcf829770ab10da237`; the uncompressed JSON SHA-256 is `38d6d08f5ef1c5628f97427db794440bd214c4c748696996bcab4c4b8e895b72`.

The primary case is `quad_n8_p65521`. Ratios below are medians of five paired reference/candidate ratios, with exact five-pair bootstrap 95% intervals. `Current` is the existing deferred 32-bit default; `wide` is the original 64-bit arithmetic. The primary marginal median complete-call times on the frozen seed were 2303.42 ms wide, 1391.04 ms current and 1313.13 ms candidate. Pair ratios use only calls from the same host, seed and phase, so they need not equal ratios of these marginal medians.

| Seed workload | Current/candidate | 95% interval | Wide/candidate | 95% interval | A/A range |
| --- | ---: | ---: | ---: | ---: | ---: |
| Frozen | 1.059× | 1.055–1.062× | 1.756× | 1.745–1.766× | 0.994–1.009× |
| Holdout A | 1.059× | 1.054–1.067× | 1.761× | 1.753–1.778× | 0.995–1.023× |
| Holdout B | 1.058× | 1.047–1.068× | 1.759× | 1.752–1.773× | 0.995–1.001× |
| Fresh holdout C | 1.058× | 1.055–1.072× | 1.748× | 1.743–1.759× | 0.991–1.006× |

The primary current/candidate gain is below the protocol's 1.15× promotion threshold, and the direct wide/candidate result is below the user's 2× complete-call target. The `quart_n3_p31` holdout A cell (0.9912× versus A/A minimum 0.9979×) and `cubic_n3_p29` fresh holdout C cell (0.9869× versus A/A minimum 0.9939×) also miss the no-smaller-regression gate. The original deferred scalar kernel remains the default. `F4_FP_LAZY_AVX2=1` stays available as an explicit opt-in for further research; no four-thread control or timing-driven repeat is run. The one-off CI workflow was removed after this receipt was archived.

These are complete prime-field F4 solver-call timings only. They do not establish a 2× complete-call gain, one-target IC online time or DLP speedup.
