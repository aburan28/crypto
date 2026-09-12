# Stage 40: charged parallel precomputation control

Stage 39 reduced the fixed algebraic n=41 precompute to 4.310 seconds on one CPU. Of that, 3.893 seconds was pair-table construction, windowed relation collection, relation verification, and log solving. These paths already use Rayon, but the admitted run confined the whole scientific process to one CPU.

Stage 40 measures the available wall reduction from four CPUs while retaining total core-seconds. It runs an ABBA sequence—one core, four cores, four cores, one core—from one clean binary. Every cell uses the exact same fixed algebraic factor base, relation seed 41431, window 149, and public hash-to-curve target seeds 41401 through 41405. Each workflow starts from an empty directory and independently constructs the base, builds the pair table, derives all factor-base logs from relations, solves five unknown-scalar targets, and runs the same-target signed-Frobenius rho control.

The verifier requires the factor base, complete log table, relation/trial/lookup counts, target coordinates, recovered scalars, descent trials, rho iterations, rho additions, and verification flags to be identical across all four cells. Timing is summarized by the median of the two one-core cells and the median of the two four-core cells. The process meter and outer sampler retain wall time, total CPU, average parallelism, and peak process-tree RSS for every cell; the result also charges the clean four-core build and all four scientific cells.

This experiment can establish a finite parallel wall improvement. It does not replace the single-core Stage 39 measurement, erase parallel core cost, supply licensed Magma or external reproduction, change the asymptotic exponent, or establish a Koblitz index-calculus SOTA result.
