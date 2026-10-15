#set page(paper: "us-letter", margin: 0.72in)
#set text(size: 9.5pt)

= Run report: factor-base hyperparameter sweep

Run started 2026-10-14 18:04 UTC (`data/run_start_utc.txt`). The command was
`cargo run --release --example fb_hyperparam_sweep -- research/fb_hyperparam_sweep_20261014/data`
at source commit `31199a70a`, with the protocol frozen in that commit (`PROTOCOL.md`). The master seed was
`0xFB5EE920261014`. Host: 4 vCPU Intel Xeon @ 2.10 GHz, x86_64, shared, L0 isolation, all work on 1
thread. Total wall time 291 s. Raw records are in `../../data/sweep.json` and `../../data/console.log`.

== What ran

*Binary arm.* 208 cells: $n = 8 dots 16$, $n' = 2 dots min(6, n-1)$, basis family
coordinate / random / subfield ($n' | n$), curve class generic $b$ (4 curves) / Koblitz $b = 1$.
Per cell: 4,096 yield targets, 256 brute-force cross-checks, 64 planted targets, 16 non-decomposable
targets for $d_"ff"$ and $D^*$, 4 decomposable calibration targets, $d_max = 6$.
208 of 208 cells completed.

*Prime arm.* 20-, 24- and 28-bit $a = -3$ curves with 3 secrets each. Each secret gets one rho solve, then
`ic_solve` with $m = 2$ (FB 64 … 16,384) and $m = 3$ (FB 64 … 4,096). That is 90 solves requested and 81 completed
with the correct logarithm. 9 returned no solve: $m = 3$ at FB 1,024 and 4,096 (20 bits), and at FB 4,096 (24 bits). `ic_solve` does not
distinguish an exhausted relation budget from a linear-algebra failure, so the reason is recorded as
"relation budget exhausted or LA failed".

== Controls (measured)

#table(columns: (2fr, 1fr), inset: 4pt,
  [Planted targets found by root finding], [13,312 / 13,312],
  [Planted targets found by brute force (first 16 per cell)], [3,328 / 3,328],
  [Root finding = brute force on cross-checked targets], [53,248 / 53,248],
  [Decomposable calibration targets that refuted], [0 / 832],
  [Non-decomposable targets refuted by $d_max$], [3,264 / 3,264],
)

== Measured values

#image("../../figures/yield_vs_heuristic.svg", width: 100%)
Yield per cell (4,096 targets, Wilson 95%) against $|V|^2 slash 2^n$, which is a derived reference.

#image("../../figures/dff_vs_dstar.svg", width: 100%)
Mean ± sd of $d_"ff"$ and $D^*$ over 16 non-decomposable targets per cell, for generic $b$ and $n in {10,12,14,16}$.
Every cell is in `sweep.json`.

#table(columns: (0.6fr, 1.5fr, 0.8fr, 0.7fr, 1.4fr), inset: 3pt,
  [bits], [method], [FB size], [solved], [wall ms (min–max of 3)],
  [20], [rho (negation, 1 thread)], [—], [3/3], [7.3–9.2],
  [20], [IC m=2], [64], [3/3], [20.0–27.7],
  [20], [IC m=2], [256], [3/3], [17.4–21.0],
  [20], [IC m=2], [1024], [3/3], [24.7–30.0],
  [20], [IC m=2], [4096], [3/3], [292.0–320.3],
  [20], [IC m=2], [16384], [3/3], [7578.3–7726.9],
  [20], [IC m=3], [64], [3/3], [1.7–2.1],
  [20], [IC m=3], [256], [3/3], [5.8–7.8],
  [20], [IC m=3], [1024], [0/3], [187.4–244.4],
  [20], [IC m=3], [4096], [0/3], [2610.1–2889.8],
  [24], [rho (negation, 1 thread)], [—], [3/3], [9.8–10.8],
  [24], [IC m=2], [64], [3/3], [266.9–380.5],
  [24], [IC m=2], [256], [3/3], [245.2–266.3],
  [24], [IC m=2], [1024], [3/3], [257.9–315.4],
  [24], [IC m=2], [4096], [3/3], [566.3–614.8],
  [24], [IC m=2], [16384], [3/3], [8022.9–8050.4],
  [24], [IC m=3], [64], [3/3], [21.8–22.8],
  [24], [IC m=3], [256], [3/3], [13.7–16.4],
  [24], [IC m=3], [1024], [3/3], [195.5–208.2],
  [24], [IC m=3], [4096], [0/3], [11071.4–12757.4],
  [28], [rho (negation, 1 thread)], [—], [3/3], [15.5–16.0],
  [28], [IC m=2], [64], [3/3], [4216.7–4994.7],
  [28], [IC m=2], [256], [3/3], [3792.1–4131.1],
  [28], [IC m=2], [1024], [3/3], [3588.3–3849.9],
  [28], [IC m=2], [4096], [3/3], [3694.6–4066.5],
  [28], [IC m=2], [16384], [3/3], [10903.4–11590.1],
  [28], [IC m=3], [64], [3/3], [294.1–422.5],
  [28], [IC m=3], [256], [3/3], [122.6–137.8],
  [28], [IC m=3], [1024], [3/3], [284.5–328.5],
  [28], [IC m=3], [4096], [3/3], [9643.6–10475.4],
)
#image("../../figures/prime_control.svg", width: 100%)

Every number above is measured, except the heuristic reference line. The derived cost model and the protocol's
decision rules are in the results report (`../../README.md`, `../../REPORT.pdf`).
