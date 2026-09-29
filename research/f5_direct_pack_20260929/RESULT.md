# Direct packed F5 rows: a measured build gain, 2× complete-call target still open

The [frozen CI run 36527288378](https://github.com/aburan28/crypto/actions/runs/36527288378)
completed all 112 one-thread processes on an AMD EPYC 7763 with AVX2 and
BMI2. The direct path activated in all seven cases on all four seeds. Every
process succeeded. Default, prior and new arms matched canonical row-space
fingerprint, rank, F4/selected row counts, column count, pruning and
criterion word operations. Prior and new also matched raw row fingerprint,
output term count and reduction word operations. The F5 exactness tests and
the local [ARM64 preflight](preflight/README.md) passed.

The complete 1,041,815-byte JSON receipt is committed as
[runs/36527288378-t1.json.gz](runs/36527288378-t1.json.gz) (82,886 bytes;
compressed SHA-256
`3f8cfd3a802d657773c4456b4d39a1dd71ecd3b12c922c663504a22f9f77ac5f`,
uncompressed SHA-256
`ffb302a0d59664e1652ed2738c61e3afc9fa53b2534c1e503c8ad517afdb5f88`).
It retains every process output/status, source and release-binary digests,
host features and load, affinity, phase timings and signatures. The CI
artifact also retains the release binary.

The table gives medians of five paired **complete-call** ratios on the
primary `f5_n24_m24_d4` case. Intervals are exact five-pair bootstrap 95%
intervals. The default returns reduced rows; prior and new both use
selective echelon output, fused row counting and the AVX2 table XOR option.
New adds direct packed-row construction.

| Seed | Default / prior | Default / new (95% interval) | Prior / new (95% interval) | Default A/A range |
| --- | ---: | ---: | ---: | ---: |
| Frozen | 1.619× | **1.955× (1.932–1.964×)** | **1.210× (1.193–1.216×)** | 1.000–1.011× |
| Holdout A | 1.601× | 1.945× (1.553–1.976×) | 1.220× (0.955–1.230×) | 0.996–1.003× |
| Holdout B | 1.655× | 1.993× (1.972–2.014×) | 1.206× (1.198–1.210×) | 0.995–1.006× |
| Holdout C | 1.649× | 2.003× (1.975–2.034×) | 1.214× (1.198–1.225×) | 0.996–1.018× |

On the frozen primary, marginal median call times were 251.74 ms default,
155.65 ms prior and 128.55 ms new. Prior/new build medians were 29.08/7.78
ms; reduction was 63.86/63.98 ms and unpacking 62.20/56.12 ms. These
phase medians do not sum to the marginal call medians because they come from
different samples. The build-only ratio was 3.753×. No smaller case's
default/new complete-call median fell below its own A/A minimum. Holdout A
had a timing outlier, reflected in its wider interval; its median still
clears the predeclared holdout check.

The frozen **incremental** gate passes: prior/new is above 1.10 including
its lower interval bound, all three holdout medians exceed their A/A maxima,
and there is no smaller-case regression by the stated rule. The direct 2×
gate does not pass because the frozen default/new median and lower bound are
below 2.00, although one holdout median exceeds it. Direct packing remains
an explicit opt-in; no default change or four-thread claim follows. The
one-off timing workflow was removed after archiving this receipt. This is
a matrix-F5 solver-call stage diagnostic, not an IC online-time or DLP
speedup.
