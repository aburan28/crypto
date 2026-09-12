# Stage 40: charged parallel precomputation control

Stage 39 reduced the fixed algebraic n=41 precompute to 4.310 seconds on one CPU. Of that, 3.893 seconds was pair-table construction, windowed relation collection, relation verification, and log solving. These paths already use Rayon, but the admitted run confined the whole scientific process to one CPU.

Stage 40 measures the available wall reduction from four CPUs while retaining total core-seconds. It runs an ABBA sequence—one core, four cores, four cores, one core—from one clean binary. Every cell uses the exact same fixed algebraic factor base, relation seed 41431, window 149, and public hash-to-curve target seeds 41401 through 41405. Each workflow starts from an empty directory and independently constructs the base, builds the pair table, derives all factor-base logs from relations, solves five unknown-scalar targets, and runs the same-target signed-Frobenius rho control.

The verifier requires the factor base, complete log table, relation/trial/lookup counts, target coordinates, recovered scalars, descent trials, rho iterations, rho additions, and verification flags to be identical across all four cells. Timing is summarized by the median of the two one-core cells and the median of the two four-core cells. The process meter and outer sampler retain wall time, total CPU, average parallelism, and peak process-tree RSS for every cell; the result also charges the clean four-core build and all four scientific cells.

This experiment can establish a finite parallel wall improvement. It does not replace the single-core Stage 39 measurement, erase parallel core cost, supply licensed Magma or external reproduction, change the asymptotic exponent, or establish a Koblitz index-calculus SOTA result.

## Hosted result

[Run `34700889826`](https://github.com/aburan28/crypto/actions/runs/34700889826) completed from merge commit `62eb3e1fc008cf33bb9e910d7a61993cd8b6921f` and passed a fresh replay after download. All four cells produced identical factor bases, log columns, 63 relations from 147,456 probes, 21,970,944 charged summand lookups, target points, recovered scalars, descent-trial counts, rho iterations, rho additions, and verification flags.

| Metric | One core median | Four core median | Change |
|---|---:|---:|---:|
| Precompute wall | 5.339 s | **2.011 s** | **2.66x faster** |
| Five-target IC wall | 5.978 s | **2.480 s** | **2.41x faster** |
| Process CPU | 5.257 core-s | 6.766 core-s | 1.287x cost |
| Amortized IC/rho | 6.59x slower | **3.61x slower** | about 1.83x better on the same instance |

The four-core cells averaged about 2.38 active cores over their complete process wall because target descent and rho retain serial work. Their sampled process-tree peaks were 317–321 MB, compared with 307–308 MB for the one-core cells. Charging the clean build and every ABBA cell gives 113.328 seconds sequential wall, 339.007 core-seconds, and 1,701,056,512 bytes maximum sampled process-tree RSS. A fresh build plus one median four-core cell is still about 139.55 times the matched rho wall.

Parallelism reduces the algorithm-only wall ratio from Stage 39's 5.03 to 3.61, and by about eightfold from Stage 35's 28.86. It spends 28.7 percent more process CPU than the one-core cell. The result therefore improves latency without satisfying the full-cost, Magma, external-review, asymptotic, or SOTA gates.
