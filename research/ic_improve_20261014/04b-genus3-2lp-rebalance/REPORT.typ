#set page(paper: "a4", margin: 1.8cm, numbering: "1")
#set text(size: 9.5pt)
#set heading(numbering: "1.")
#show table: set text(size: 8pt)

#align(center)[
  #text(size: 14pt, weight: "bold")[Genus-3 IC-2LP after the cubic GHS transfer: is the departure from $q^(4\/3)$ an artifact of the factor base or the linear algebra?]
  #v(2pt)
  research/ic_improve_20261014/04b-genus3-2lp-rebalance · 2026-10-15 · Conductor T-1 · run-20261015T0200Z \
  Protocol frozen at commit 9752513ea (sha256 633e69bb…); κ choice frozen at 288e74293
]

#block(fill: luma(240), inset: 6pt, radius: 3pt)[
  *Rebuilt baseline.* The study this follows up (`04-genus3-2lp`, its protocol, `RULES.md`, `examples/ic_improve_04_*`) is not in the repository (`problems/ic-improve-04-genus3-2lp-study-missing.md`). At the requester's choice the IC-2LP baseline was rebuilt; every arm here is paired with *that* baseline on the same instances, seeds and unit, not with the original's records.
]

= Answer

On this rebuild the departure from 4/3 is *not* an artifact of either choice.

- *Baseline (BL, κ = 2, tree-cycle elimination + Lanczos):* total cost fits $q^(1.383)$ over $2^14$–$2^20$ (99% [1.382, 1.384]); local exponents trend +0.003 per $log_2 q$ (99% [0.002, 0.004]). The relation phase fits $q^(1.332)$ (GTTD's 4/3); the linear algebra fits $q^(1.447)$, driven by the Lanczos row weight growing 7.8 → 14.4 from $2^14$ to $2^20$.
- *A1, factor base retuned per size* (κ ∈ {0.5, 0.7, 1, 1.4, 2} on training curves 5–7): verdict *reduces the exponent*, E = 1.373 [1.372, 1.374], 0.010 below BL and still 0.039 above 4/3. Tuning picked κ = 1.4 at $2^17$–$2^20$; paired saving 0.8% → 3.9%.
- *A2, structured Gaussian elimination before Lanczos* (κ = 2, BL's own relations): verdict *raises the exponent*, E = 1.419 [1.417, 1.421]; cheaper at every size (ratio 0.78 at $2^14$, 0.91 at $2^20$) but the saving shrinks as BL's weight catches up with A2's bounded ≈ 15–16.
- A3 (secondary, κ tuned on A2's pipeline) chose κ = 2 everywhere and is A2.

All 360 logarithms (180 test solves × 2 linear-algebra routes) satisfy $[log]P = Q$ and equal the planted value; 210 tuning solves and 20 rho runs also verified; no job failed (max 119 MB, 490 s).

= Unit, instances, arms

Unit: multiplications in $FF_q$ plus multiplications mod ℓ over the whole method (relation setup, walk, decomposition; elimination; Lanczos; verification); $S = "total" \/ sqrt(ell)$. Curves: hyperelliptic GHS cubic covers ($m = 3, r = 5$, T-2463's construction), $y^2 = f(x)$, deg 7 over $FF_q$, $q$ the largest prime below $2^k$, $k = 10..20$, $\#"Jac"(C)(FF_q) = h ell$ with $h in {4, 8}$; 7 curves per size, seed 20261016; test curves 1–4 × targets 1–3 (12 instances per size). Rho: r-adding walk with Brent cycle detection on Jac(C), $k <= 14$, 4 walks per size.

= Results

#table(columns: 9, align: right,
  [*k*], [*BL total*], [*A1 (κ\*)*], [*A2 total*], [*S#sub[BL]*], [*S#sub[A1]*], [*S#sub[A2]*], [*S#sub[rho]* (n=4)], [*BL/rho* raw, renorm.],
  [10], [3.142e6], [—], [2.559e6], [192.6], [—], [156.9], [550.0], [0.35, 0.58],
  [11], [6.765e6], [—], [5.323e6], [207.8], [—], [163.5], [739.1], [0.28, 0.63],
  [12], [1.560e7], [—], [1.199e7], [119.1], [—], [91.5], [439.7], [0.27, 0.35],
  [13], [3.878e7], [—], [2.997e7], [148.0], [—], [114.3], [582.6], [0.25, 0.42],
  [14], [9.686e7], [9.686e7 (2)], [7.566e7], [92.4], [92.4], [72.2], [571.4], [0.16, 0.25],
  [15], [2.472e8], [2.472e8 (2)], [1.974e8], [83.4], [83.4], [66.6], [—], [—],
  [16], [6.731e8], [6.731e8 (2)], [5.575e8], [80.3], [80.3], [66.5], [—], [—],
  [17], [1.700e9], [1.687e9 (1.4)], [1.436e9], [101.3], [100.5], [85.6], [—], [—],
  [18], [4.387e9], [4.301e9 (1.4)], [3.801e9], [92.4], [90.6], [80.1], [—], [—],
  [19], [1.167e10], [1.136e10 (1.4)], [1.036e10], [86.9], [84.6], [77.2], [—], [—],
  [20], [3.049e10], [2.930e10 (1.4)], [2.763e10], [56.8], [54.6], [51.5], [—], [—],
)
Mean over 12 test instances (95% intervals in tables.md). S moves with ℓ across curves. Rho walks took 0.30–4.0 steps/√ℓ (means 1.6–2.8 vs ideal 1.25); "renorm." rescales rho to $sqrt(pi\/2)$.

#table(columns: 6, align: (left, left, right, right, right, right),
  [*arm*], [*phase*], [*E (fit 14–20)*, 95%, 99%], [*E − 4/3*], [*E − E#sub[BL]* (99%)], [*trend T* (99%)],
  [BL], [total], [*1.383* [1.382, 1.384] [1.382, 1.384]], [+0.049], [—], [+0.003 [0.002, 0.004]],
  [BL], [relation], [1.332 [1.331, 1.334]], [−0.001], [—], [−0.001 [−0.003, 0.002]],
  [BL], [lin. algebra], [1.447 [1.445, 1.448]], [+0.113], [—], [+0.002 [0.000, 0.005]],
  [A1], [total], [*1.373* [1.372, 1.374] [1.372, 1.374]], [+0.039], [−0.010 [−0.010, −0.010]], [−0.001 [−0.002, 0.000]],
  [A2], [total], [*1.419* [1.418, 1.420] [1.417, 1.421]], [+0.086], [+0.036 [+0.035, +0.037]], [+0.001 [−0.001, 0.003]],
  [A2], [lin. algebra], [1.599 [1.596, 1.602]], [+0.266], [+0.152 [+0.149, +0.155]], [−0.027 [−0.032, −0.023]],
)

#table(columns: 8, align: right,
  [*k*], [14], [15], [16], [17], [18], [19], [20],
  [A1/BL], [1 (same run)], [1], [1], [0.992], [0.980], [0.974], [*0.961* [0.959, 0.963]],
  [A2/BL], [0.781], [0.798], [0.828], [0.844], [0.867], [0.888], [*0.906* [0.905, 0.907]],
  [A2/BL lin. alg.], [0.429], [0.505], [0.578], [0.646], [0.712], [0.766], [0.812],
)
Paired cost ratios, geometric mean over the 12 instances (99% interval where shown; all intervals in tables.md). A2/BL at k = 10–13: 0.815, 0.787, 0.768, 0.773.

#figure(image("figures/local_exponents.svg", width: 92%), caption: [Local exponents of mean total cost per size step, 12 test instances, 95% paired bootstrap intervals (B = 10000). A3 coincides with A2 (κ\*\* = 2). run-20261015T0200Z.])
#figure(image("figures/phase_exponents.svg", width: 100%), caption: [Local exponents by phase: the relation phase sits near 4/3 from $2^16$; BL's linear algebra runs at ≈ 1.45; A2's linear algebra starts steeper (it removes a falling fraction of the matrix). The 15→16 bump is a per-trial decomposition-cost step (300 → 378 mults).])
#figure(image("figures/row_weight.svg", width: 88%), caption: [Nonzeros per row of the matrix entering Lanczos. BL grows 4.9 → 14.4; A2 stays 14.6–16.2 from $2^14$; A1 at κ = 1.4 is heavier (interpretation: longer cycles in a sparser large-prime graph).])
#figure(image("figures/paired_ratios.svg", width: 88%), caption: [Paired cost ratios against BL with 95% intervals (12 instances per size).])
#figure(image("figures/total_cost.svg", width: 88%), caption: [Mean total multiplications against q with a $q^(4\/3)$ guide through BL at $2^14$ (guide only, not a fit).])
#figure(image("figures/tuning.svg", width: 88%), caption: [Tuning on training curves 5–7 × targets 1–2: mean total of the BL pipeline relative to κ = 2. The interior optimum κ = 1.4 appears from $2^17$; at $2^14$–$2^16$ the optimum is the grid edge κ = 2.])

= Decision (PROTOCOL §7) and controls

D0 baseline departs: *yes* (99% E#sub[BL] [1.382, 1.384] > 4/3; 99% T#sub[BL] [0.0016, 0.0042] > 0). A1: *reduces the exponent* (not removes: E's 99% interval excludes 4/3). A2: *raises the exponent*. Answer: *not* an artifact of either choice on this rebuild. Power: 95% half-width of E#sub[BL] 0.0009 (≤ 0.03).

Controls: A2's relation phase is BL's own record on all 132 instances (pairing). A/A: E#sub[BL] on target 1 minus targets 2–3 = +0.002, 95% [0.001, 0.004]. This excludes 0 but is 25× smaller than the departure, and every verdict's margin is ≥ 0.010; read as the bootstrap understating between-target spread at the 0.002 level with 4 curves per size. Lanczos accounting equals its closed form on all 360 runs. One operator error (an empty job list, not a test job) is logged in results/orchestration_errors.jsonl.

= Interpretation (labelled)

Lanczos costs ≈ $F^2 (2w + 7)$; with $F prop q^(2\/3)$ and $w prop q^delta$, measured $delta approx 0.15$ over $2^14$–$2^20$ (tree-cycle relations lengthen as the BFS tree deepens), so the linear algebra runs above $q^(4\/3)$; the constant term dilutes it to the measured 1.447. A heuristic balance with relation cost $prop q^2 \/ F$ gives optimal $F prop q^(2\/3) w^(-1\/3)$ and total $q^(4\/3) w^(1\/3)$, i.e. ≈ 4/3 + δ/3 ≈ 1.38 even with continuous retuning; A1 measures 1.373 (model untested beyond the grid). A2's Markowitz merges stop near the weight BL reaches by $2^20$, so its saving is a constant that shrinks with q; extrapolating the ratio trend (+0.02 per size, not measured) puts parity a few sizes past $2^20$. Only one SGE design was tested. Classification (AGENTS.md §3): A1 and A2 are engineering; neither is an advance.

Against the request's description of the original (qualitative, unpaired): local exponents 1.25 → 1.46 (original) vs 1.11 → 1.39 here; linear algebra $q^(1.52)$ vs $q^(1.447)$; relation $q^(1.30)$ vs $q^(1.332)$; row weight 7 → 25 vs 4.9 → 14.4. Same direction, smaller magnitude.

= Limits

Õ hides log factors (decomposition costs ≈ log q per trial). The κ grid is the request's; the optimum sits on its edge κ = 2 for A2 at every size and for BL at $2^14$–$2^16$. The transfer is not executed per target (T-2463 verified $"Jac"(C)[ell] = phi(E_0(FF_(q^3))[ell])$); curve construction and group order are not counted. Rho is measured only at $k <= 14$, n = 4, on Jac(C). Host: shared 4-vCPU container, no isolation; wall time is not a metric. No new theorem, so no Lean obligation.
