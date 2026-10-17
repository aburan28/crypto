#set page(paper: "a4", margin: 1.6cm)
#set text(size: 9.5pt)
#set par(justify: true)
#show table: set text(size: 8.5pt)

= Run report: lit2026 W12 zero-minor scaling, run-20261012T0009Z

*Descriptive run report (observations only).* The interpretation is in the study's `README.md` and `REPORT.pdf`.

== What ran

- *Binary.* `examples/lit2026_w12_zero_minor.rs`, built in release mode.
  - The scaling run and its replay used commit `c70d57ec7`; binary sha256 `02aaf7cb…2740`.
  - The defect sweep, `refit` and `compare` used commit `06b032f0c`. That commit changes only the fit guard and adds two analysis modes; the solver code is identical.
  - The rerun d-sweep gave the same per-run counts as the sweep that crashed after its last run (`dsweep-failed-c70d57ec7.log`, `diff` clean).
- *Host.* Shared M4 Pro (macOS, Darwin 25.6.0), isolation *L0*. `uptime` was recorded with every sample in `scaling.json` (load averages 79–106 at the start, 97–104 at the end). Wall time is reported only as L0 medians. The gate is counted operations.
- *Commands* (run inside this directory):
  - `lit2026_w12_zero_minor run --bits 16,18,20,22,24,26 --targets 8 --rho-seeds 4 --out scaling.json`
  - `lit2026_w12_zero_minor dsweep --bits 20 --ds 2,3,4 --targets 8 --out dsweep.json`
  - `lit2026_w12_zero_minor refit --in scaling.json --out fits.json`
  - The run repeated once more, then `compare --a scaling.json --b <replay>.json` (`replay-compare.txt`).
- *Curves.* The registered prime-order ladder `find_prime_order_curve(bits, 59297)`, cofactor 1. The binary rebuilds each curve from its registry parameters and asserts the ICV1 slug:
  - `icv1-fp16-t295-8d3c3165`, `icv1-fp18-tm337-d28d5e09`, `icv1-fp20-t727-cd198a38`;
  - `icv1-fp22-tm1385-475dcb5f`, `icv1-fp24-t1059-6df599df`, `icv1-fp26-tm1775-7e8fb6df`.
  - The generator is the point with the least abscissa and the smaller ordinate.
- *Targets.* 8 per curve, planted `k` from seed 20261011. The solvers never see `k`. Each answer is verified outside the solver by `[m]G = Q`.
- *Arms, interleaved target by target.*
  - `zero_minor.guess_determine`: 1 seed per target.
  - `rho.negation` (`ic_boundary::rho_walk_with` + `NegationClasses` + `RhoWalk::negation()`, cap 64×floor): 4 walk seeds per target.
- *Unit.* $F_p$ multiplications (squarings included) plus inversions, each counted once. Additions are free; this is the convention of `ic_boundary::FieldOps`. Rho's ops come from `PrimeCurve`'s counting group law. The zero-minor ops come from the example's counted `Fp`, plus `PrimeCurve`'s count for the $2ell$ scalar multiplications.

== Parameters chosen by the paper's rules

#table(columns: 8,
  [bits], [$r$], [$n'$], [$ell = 3n'$], [$ell'$], [$d$ (Eq. 4 ≥ 0.9)], [Eq. (4) success], [distinct maximal minors per guess],
  [16], [53117], [14], [42], [21], [3], [1.0000], [1.81e6],
  [18], [212437], [16], [48], [24], [3], [1.0000], [4.17e6],
  [20], [562427], [18], [54], [27], [3], [1.0000], [8.68e6],
  [22], [2435809], [20], [60], [30], [3], [0.9999], [1.67e7],
  [24], [16482251], [22], [66], [33], [3], [0.9059], [3.00e7],
  [26], [42653927], [24], [72], [36], [4], [1.0000], [3.52e9],
)

== Results: every run solved and verified (48/48 zero-minor, 192/192 rho)

#table(columns: 7,
  [bits], [rho.negation median ops (min–max), 32 runs], [zero-minor median ops (min–max), 8 runs], [zero-minor search-phase median], [paired ratio zm/rho, geomean (95% bootstrap)], [zm median wall ms (L0)], [rho median wall ms (L0)],
  [16], [1.92e3 (1.14e3–6.30e3)], [5.53e5 (3.14e5–8.30e5)], [2.22e5], [260 (215–319)], [9.0], [0.067],
  [18], [2.98e3 (1.43e3–5.13e3)], [1.48e6 (6.74e5–2.74e6)], [9.43e5], [455 (326–626)], [75], [0.11],
  [20], [4.62e3 (1.59e3–8.64e3)], [2.54e6 (9.15e5–6.30e6)], [1.71e6], [563 (350–883)], [262], [0.19],
  [22], [7.93e3 (3.12e3–2.05e4)], [9.03e6 (1.42e6–2.13e7)], [7.79e6], [785 (349–1633)], [1660], [0.44],
  [24], [1.67e4 (4.51e3–4.61e4)], [2.16e7 (6.92e6–5.82e7)], [1.98e7], [1109 (706–1743)], [4776], [0.66],
  [26], [2.49e4 (6.48e3–6.29e4)], [2.71e8 (7.76e6–5.99e8)], [2.68e8], [6715 (2263–15336)], [51628], [1.12],
)

Fits of $log_2$ cost against $log_2 r$ (`fits.json`; OLS over runs, 95% bootstrap resampling runs within each curve, t interval in brackets after it):

#table(columns: 4,
  [series], [slope], [95% bootstrap], [95% t],
  [rho.negation, field ops (n = 192)], [0.377], [0.350–0.404], [0.350–0.404],
  [rho.negation, group ops (n = 192)], [0.400], [0.371–0.430], [0.371–0.430],
  [zero-minor, field ops, all phases (n = 48)], [0.783], [0.669–0.879], [0.669–0.897],
  [zero-minor, field ops, search phases only (n = 47)], [0.905], [0.769–1.021], [0.746–1.064],
  [zero-minor, $d$-subsets covered (n = 47)], [0.929], [0.793–1.048], [0.763–1.095],
  [zero-minor, d = 3 runs only, 16–24 bits (n = 40)], [0.631], [0.542–0.717], [0.528–0.734],
)

The "search phases only" and "covered" fits have n = 47. One 16-bit target was solved by a zero entry of the dense part, before any determine step.

*Covered $d$-subsets divided by $r$*, median per curve: 1.37, 1.32, 0.92, 0.99, 0.37, 2.13 (16 → 26 bits); the ranges are in `fits.json`.

The run's own fits in `scaling.json` group by method string, and that string carries $d$. Its `d=4` entries rest on a single curve size and are degenerate; use `fits.json` instead.

*Defect sweep on `icv1-fp20-t727-cd198a38`* (8 targets per $d$, same targets). Median ops (min–max):

- $d=2$: 1.21e7 (9.58e5–4.25e7)
- $d=3$: 2.54e6 (9.15e5–6.30e6)
- $d=4$: 3.32e6 (9.11e5–9.88e6)

All 24 runs verified.

*Replay.* A second full execution of the scaling run recorded identical method, seed, op counts, group ops, search ops, verification and recovered value in 240 of 240 rows (`replay-compare.txt`).

#figure(image("scaling.svg", width: 100%), caption: [Run-20261012T0009Z. (a) Field operations per solved target, against $log_2 r$: rho.negation (32 runs per curve) and zero-minor (8 runs per curve; circles $d=3$, squares $d=4$), with the OLS fits from `fits.json`. (b) $d$-subsets of the signature matrix covered before the zero minor, with the line covered $= r$ for reference. (c) Paired per-target ratio of zero-minor to rho field ops, geometric mean with 95% bootstrap interval.])

#figure(image("dsweep.svg", width: 60%), caption: [Defect sweep, run-20261012T0009Z, 8 targets per $d$ on the 20-bit ladder curve: median field ops, with the bar spanning min to max.])

== Raw records
`scaling.json` (every run, with per-phase tallies, counters and `uptime`), `scaling.log`, `fits.json`, `dsweep.json`, `dsweep.log`, `dsweep-failed-c70d57ec7.log` (a crash in the fit step after all runs: a single-size fit gave NaN), `replay-compare.txt`, `uptime-start.txt`, `uptime-end.txt`. The plot source is `../../plot_scaling.py`; it computes nothing.
