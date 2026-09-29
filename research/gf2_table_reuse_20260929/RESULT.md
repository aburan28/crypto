# Table reuse reaches the frozen one-thread 2× F5 complete-call gate

The [frozen one-thread CI run 36528934150](https://github.com/aburan28/crypto/actions/runs/36528934150)
completed all 112 processes on one pinned thread of an AMD EPYC 7763 with
AVX2 and BMI2. The GF(2) table-reuse exactness tests passed. Every call
succeeded and matched canonical row-space fingerprint, rank, row and column
counts, pruning and criterion work across the three arms. Prior and new
also matched raw row fingerprint, output term count and reduction word
operations on all seven cases and four seeds.

The complete 1,040,838-byte JSON receipt is committed as
[runs/36528934150-t1.json.gz](runs/36528934150-t1.json.gz) (82,527 bytes;
compressed SHA-256
`318eaf8d0396b6ff60c6634c37e3654a665ae93f23b2fb368653a270e1be720d`,
uncompressed SHA-256
`87b67a46f3b445ad8a30ed31c2cfcd041b61d2870c22eec9dec21c1ca9399b44`).
It retains every process output/status, source and release-binary hashes,
CPU features, affinity, load, phase timings and signatures. The CI
artifact also retains the release binary.

The table reports medians of five paired **complete-call** ratios on the
primary `f5_n24_m24_d4` case, with exact five-pair bootstrap 95% intervals.
The default returns reduced rows. Prior uses selective echelon output,
fused row counting, direct packed rows and the AVX2 XOR option. New adds
table-buffer reuse. The output form change is explicit in these options;
all arms have the same canonical row space.

| Seed | Default / prior | Default / new (95% interval) | Prior / new (95% interval) | Default A/A range |
| --- | ---: | ---: | ---: | ---: |
| Frozen | 1.953× | **2.030× (2.014–2.119×)** | **1.035× (1.031–1.041×)** | 0.992–1.010× |
| Holdout A | 1.944× | 2.024× (1.995–2.040×) | 1.043× (1.026–1.139×) | 0.992–1.006× |
| Holdout B | 1.997× | 2.067× (2.044–2.088×) | 1.034× (1.028–1.040×) | 0.985–1.006× |
| Holdout C | 2.001× | 2.073× (2.060–2.117×) | 1.036× (1.031–1.044×) | 0.983–1.016× |

On the frozen primary, marginal median times were 254.87 ms default,
130.50 ms prior and 125.89 ms new. Prior/new reduction medians were
64.66/60.30 ms; build medians were 7.69/7.82 ms and unpack medians were
57.33/56.95 ms. The reduction-phase paired ratio was 1.072×. The
complete-call ratio is measured directly and is not a product of these
phase ratios. No smaller case's default/new complete-call median was
below its own A/A minimum.

The predeclared one-thread 2× gate passes: the frozen primary median and
lower interval bound exceed 2.00, each holdout median exceeds its A/A
maximum, and the smaller-case guard passes. This establishes a **2.030×
one-thread opt-in matrix-F5 solver-call gain on the measured host**.

## Preregistered four-thread control

The [four-thread CI run 36530113651](https://github.com/aburan28/crypto/actions/runs/36530113651)
completed the fixed [four-thread protocol](FOUR_THREAD_PROTOCOL.md): 112
processes, four pinned allowed CPUs, and `RAYON_NUM_THREADS=4`. All calls
were exact and successful. Its complete 1,040,578-byte JSON receipt is
[runs/36530113651-t4.json.gz](runs/36530113651-t4.json.gz) (82,874 bytes;
compressed SHA-256
`d0826506b435a4ef4917eb6e9a197e917fdf7105e508e5d24974d63043e045ff`,
uncompressed SHA-256
`75d31482a5ba2f725d7cbcde34970d6d66f3f9fe18dbfc6a2a297f97798fb357`).
It has the same full process, host, source and binary accounting as the
one-thread receipt. This runner also reported an AMD EPYC 7763, but it was
a separate instance; the four-thread comparisons below are paired within
their own run, not divided by the one-thread milliseconds.

| Four-thread seed | Default / prior | Default / new (95% interval) | Prior / new (95% interval) | Default A/A range |
| --- | ---: | ---: | ---: | ---: |
| Frozen | 1.668× | **1.755× (1.667–1.773×)** | 1.062× (0.983–1.082×) | 0.986–1.035× |
| Holdout A | 1.695× | 1.732× (1.723–1.778×) | 1.029× (1.011–1.045×) | 0.979–1.019× |
| Holdout B | 1.701× | 1.749× (1.721–1.764×) | 1.021× (1.006–1.054×) | 0.973–1.031× |
| Holdout C | 1.705× | 1.777× (1.427–1.813×) | 1.043× (0.846–1.050×) | 0.980–1.017× |

On the four-thread frozen primary, marginal median calls were 167.23 ms
default, 100.41 ms prior and 95.88 ms new. No smaller case's default/new
median fell below its A/A minimum. The combined mode passes the control's
useful-gain gate, but the frozen four-thread median and lower bound miss
2.00. The standalone prior/new table-reuse effect at four threads is also
uncertain because its lower interval bound is below 1.00. The four-thread
control therefore does not extend the 2× claim or justify a default
change. The opt-in remains available; the one-off four-thread workflow was
removed after archiving its receipt.

Neither solver-call comparison establishes one-target IC online time or a
DLP speedup. The one-thread ratio applies to the explicitly selected
echelon-output, fused-build, direct-pack, AVX2 and table-reuse combination;
ordinary default behavior remains unchanged.
