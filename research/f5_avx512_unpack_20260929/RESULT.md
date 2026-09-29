# AVX-512 F5 row unpacking: exact output, no measured gain

The [frozen CI run 36525322701](https://github.com/aburan28/crypto/actions/runs/36525322701)
completed all 112 processes on one pinned thread of an AMD EPYC 9V45 with
AVX-512F, AVX2 and BMI2. The AVX-512 boundary and F5 row-space tests passed.
Every process succeeded; the default, prior combined, and new combined arms
matched canonical row-space fingerprint, rank, pruning and criterion work on
all seven cases and all four seeds. Prior and new arms also matched raw row
fingerprint and output term count exactly. The full 887,191-byte JSON receipt
is committed as [runs/36525322701-t1.json.gz](runs/36525322701-t1.json.gz)
(80,307 bytes, compressed SHA-256
`68c8255e75c6911e34a251734b937b418aa25af0f54e48069f7972108aed1ac8`,
uncompressed SHA-256
`54408a082f900016cbb2fe7ab762f70637bc04405b2a2440601a2ceb9b7d0d38`).
It retains every process's output/status, the release binary digest, source
digests, host features, affinity, load, and per-case timing. The CI artifact
also retains the release binary.

The table reports medians of five paired complete-call ratios on the primary
`f5_n24_m24_d4` case, with exact five-pair bootstrap 95% intervals. The
default is the current reduced output; prior and new both use selective
echelon output, fused row building and the AVX2 XOR option. New alone adds
AVX-512 unpacking.

| Seed | Default / prior | Default / new (95% interval) | Prior / new | Default A/A range |
| --- | ---: | ---: | ---: | ---: |
| Frozen | 1.466× | **1.478× (1.421–1.557×)** | 1.002× | 0.997–1.047× |
| Holdout A | 1.455× | 1.477× (1.444–1.518×) | 1.016× | 0.986–1.207× |
| Holdout B | 1.476× | 1.477× (1.460–1.507×) | 1.004× | 0.839–1.002× |
| Holdout C | 1.466× | 1.460× (1.412–1.495×) | 0.996× | 0.992–1.015× |

On the frozen primary, marginal median full-call times were 167.15 ms
default, 112.63 ms prior and 113.26 ms new. Prior/new median phase times
were 18.87/18.73 ms for build, 44.13/44.50 ms for reduction, and
48.22/48.57 ms for unpacking. The prior and new forms each returned
13,734,979 monomials across 6,924 pivot rows; the default reduced form
returned 13,127,926. No smaller case's default/new complete-call median
fell below its own A/A minimum. Some A/A ranges are wide, so the small
prior/new differences are not a reliable gain; the direct primary ratio
is far below the predeclared 2× gate.

The hypothesis is rejected for this implementation. AVX-512 compress-store
does not materially accelerate the observed unpack cost, so it stays an
explicit opt-in and the default is unchanged. The one-off timing workflow
was removed after the receipt was archived. Per the frozen stop rule, no
four-thread control or timing-driven repeat is run. This is a complete
matrix-F5 solver-call stage diagnostic, not an IC online-time or DLP speedup.
