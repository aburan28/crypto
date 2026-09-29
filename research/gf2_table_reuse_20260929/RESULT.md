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
one-thread opt-in matrix-F5 solver-call gain on the measured host**. It
does not establish a four-thread gain, a default-mode improvement, one-target
IC online time or a DLP speedup. The preregistered four-thread control is
pending under [FOUR_THREAD_PROTOCOL.md](FOUR_THREAD_PROTOCOL.md); no default
change follows from the one-thread result alone.
