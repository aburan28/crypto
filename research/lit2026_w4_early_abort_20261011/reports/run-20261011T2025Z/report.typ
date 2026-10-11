#set page(paper: "a4", margin: 1.8cm)
#set text(size: 9.5pt)
#set table(stroke: 0.4pt)

= Run report: lit2026 W4 refutation-degree trace and early-abort arms, run-20261011T2025Z

*Observations only.* Interpretation is in the study report (`../../REPORT.pdf`).

== What ran

- *Binaries.* `f4_refutation_degree_probe` and `f6_ic_probe` (examples), release, built from branch `lit2026/w4-early-abort` (working tree on `origin/feat/elliptic-curve-crypto-research-88c7cc`, the commit of this PR carries the same `src/` and probe code; the summarize mode of `f4_refutation_degree_probe` was added after the runs and does not touch the solve). `f6_ic_probe.orig`: built from the program baseline worktree at `4b6e09b78` with `CARGO_TARGET_DIR=/Volumes/SSD990-2/crypto/targets/lit2026-baseline`.
- *Arms (timing).* `orig` = baseline binary; `prev` = new binary, `KIC_F4_EARLY_ABORT=0`; `cand1` = new binary, `=1` (degree-2 probe before a fresh top-degree build); `cand2` = new binary, `=2` (degree-2 probe before a child specialises its parent's bases).
- *Workloads.* Koblitz curves, standard-subspace factor base, $m = 3$, targets $[k]G$, $k = 1 dots$ in the prime-order subgroup `KoblitzCurve::new` selects, `SolverEngine::InheritedF4 { max_degree: 3 }`, engines `f4` (inherited F4, the default decomposition oracle) and `f6p` (F6-IC with the shared pair index):
  - K0/2^23 (`icv1-f2m23-t5197-69e76b73`), dim 8 (275 points), 4 targets, budget 20 000;
  - K0/2^31 (`icv1-f2m31-tm90707-c95f16f5`), dim 9 (553 points), 1 target, budget 2 000;
  - K1/2^17 (`icv1-f2m17-tm101-00378d4e`), dim 6 (63 points), 8 targets, budget 20 000;
  - trace only: K1/2^23 (`icv1-f2m23-tm5197-1f85e9e1`), dim 7, 4 targets.
- *Commands.*
  - trace: `KIC_F4_EARLY_ABORT=$k f4_refutation_degree_probe --a A --n N --dim D --targets T [--budget B] --engines f4,f6p`, for k in {0, 1, 2}, one file each `raw/trace_<cell>_mode<k>.jsonl`;
  - timing: `KIC_F4_EARLY_ABORT=$k f6_ic_probe --a A --n N --dim D --targets T [--budget B] --engines f4,f6p --reps 3`, file `raw/time_<cell>_<arm>.round<r>.jsonl`, `uptime` before each in `.uptime`;
  - summary: `f4_refutation_degree_probe --summarize raw > summary.json` (sums, medians, ratios of the recorded values; Rust).
- *Order.* Five rounds, arm order rotated: (orig, prev, cand1, cand2), (cand2, cand1, prev, orig), (prev, orig, cand2, cand1), (cand1, cand2, orig, prev), (orig, prev, cand1, cand2); within an arm the three cells in a fixed order.
- *Host and isolation.* Shared Apple M4 Pro (macOS, Darwin 25.6.0), many other sessions compiling; 1-minute load average 63.6–110.9 over the timing session (20:38–21:05 UTC). No pinning. Isolation L0. Single-threaded probes.
- *Requested / completed.* Traces: 12 requested, 12 completed. Timing: 60 probe invocations requested, 60 completed, none failed.

== Graphs

#figure(image("degrees.svg", width: 100%), caption: [Inherited F4 arm, `KIC_F4_EARLY_ABORT=0`, sums over each cell's attempts. Left: each F4 reduction (node) filed under the first from-scratch Macaulay degree (2, then 3) whose row space holds 1 (refuted) or, failing that, a forced variable; else stalled (the solver splits). Right: stage-unit word operations spent by nodes that end refuted (any cause) and by nodes that end in a split, children excluded. run-20261011T2025Z, `summary.json`.])

#figure(image("arms.svg", width: 100%), caption: [Left: `word_xors(original) / word_xors(arm)` (stage unit, deterministic, identical in all five rounds of each arm). Right: probe `total_ms` (per target best of 3 repetitions, summed) of each round (dots) and median (bar), log scale. run-20261011T2025Z, 5 rounds per arm.])

== Recorded values: where nodes are decided (f4 arm, mode 0, sums over attempts)

#table(columns: 9,
  [cell], [attempts / hits], [reductions], [splits], [refuted at deg 2], [refuted at deg 3 only], [forced, first deg], [refuted (fresh / inherited basis)], [word ops refuted / split nodes],
  [K1/2^17 d6], [8 / 1], [1201], [504], [203], [305], [189 at 3], [18 / 490], [2.60 M / 12.42 M],
  [K1/2^23 d7], [4 / 0], [1180], [500], [387], [115], [178 at 3], [6 / 496], [2.01 M / 51.26 M],
  [K0/2^23 d8], [4 / 0], [2796], [1192], [425], [770], [409 at 3], [4 / 1191], [32.25 M / 126.26 M],
  [K0/2^31 d9], [1 / 0], [1410], [602], [536], [67], [205 at 3], [2 / 601], [2.14 M / 279.44 M],
)
Derived (`summary.json`): share of degree-2 among F4 refutations 0.400 / 0.771 / 0.356 / 0.889; share of node word operations in refuted nodes 0.173 / 0.038 / 0.203 / 0.0076 (same cell order). No node was forced at degree 2. The from-scratch degree-3 space and the inherited basis refuted exactly the same nodes (`scratch_refutes_engine_does_not` = `engine_refutes_scratch_does_not` = 0 everywhere). Root systems: no attempt's root system is refuted at degree 2 or 3; at K1/2^17 none through degree 4; at the other cells the degree-4 matrix exceeded the size caps (`root_first_refutation_degree[1]`). Other refutations: 1 by linear elimination (K0/2^23), 1 (K1/2^23). F6-IC arm: the F4 bases refute no node; the geometry gate refutes them, and 66.6 / 74.0 / 70.1 / 70.7 % of its node word operations are spent in nodes that the gate eventually refutes (work done at the node before the refuting oracle call).

Probes: mode 1 ran 254 / 222 / 556 / 279 degree-2 probes and refuted 10 / 6 / 3 / 2 nodes; mode 2 ran 908 / 924 / 2111 / 1061 and refuted 193 / 381 / 422 / 534. Every solver counter and every witness is equal across modes 0, 1, 2 (one hit, K1/2^17 target 3, witness base indices [21, 48, 57]).

== Recorded values: timing arms (stage unit and wall)

#table(columns: 7,
  [cell, engine], [orig `word_xors`], [prev], [cand1 (EA=1)], [cand2 (EA=2)], [orig/cand1], [orig/cand2],
  [K1/2^17, f4], [15 015 618], [15 015 618], [14 959 088], [15 058 684], [1.00378], [0.99714],
  [K1/2^17, f6p], [12 635 914], [12 635 914], [12 636 832], [12 665 650], [0.99993], [0.99765],
  [K0/2^23, f4], [158 512 775], [158 512 775], [158 321 408], [158 846 162], [1.00121], [0.99790],
  [K0/2^23, f6p], [124 949 801], [124 949 801], [124 951 833], [125 135 729], [0.99998], [0.99851],
  [K0/2^31, f4], [281 583 026], [281 583 026], [281 350 163], [281 787 203], [1.00083], [0.99928],
  [K0/2^31, f6p], [278 237 634], [278 237 634], [278 238 729], [278 409 768], [1.00000], [0.99938],
)
`word_xors` counts the word XORs of reductions and specialisation, not the generation of Macaulay rows. Reductions, splits, hits and every per-target verdict and counter (`hit`, `exhausted`, `reductions`, `splits`, `infeasible`) are identical across all four arms and all five rounds.

#table(columns: 5,
  [cell, engine], [orig ms, median (min–max)], [prev], [cand1], [cand2],
  [K1/2^17, f4], [297.0 (231.5–752.5)], [360.9 (134.2–757.4)], [314.0 (185.9–641.2)], [256.3 (167.2–336.7)],
  [K1/2^17, f6p], [276.6 (179.8–467.7)], [121.9 (93.9–355.3)], [153.3 (113.0–242.1)], [149.1 (117.2–261.3)],
  [K0/2^23, f4], [4902.7 (1237.3–6210.8)], [5309.4 (1232.1–7334.7)], [2543.1 (1554.9–7039.1)], [3712.0 (1417.4–8124.8)],
  [K0/2^23, f6p], [2418.7 (1826.5–4304.7)], [2936.5 (1120.5–6271.0)], [1617.3 (1013.0–5156.2)], [2412.3 (1021.6–8833.4)],
  [K0/2^31, f4], [5833.6 (2078.6–8266.0)], [7484.8 (2132.7–14567.5)], [6215.3 (1762.2–11884.5)], [3433.0 (2046.7–11144.1)],
  [K0/2^31, f6p], [7781.7 (1903.5–9680.6)], [7753.0 (1855.8–8722.1)], [3963.6 (1773.6–11981.9)], [4847.6 (1825.6–7005.5)],
)
Null control, per-round wall ratio orig / prev (same tree, same stage unit; derived in `summary.json`): f4 median 1.555 (range 0.61–2.11) at K1/2^17, 0.975 (0.53–1.04) at K0/2^23, 1.104 (0.14–3.15) at K0/2^31. The round-to-round spread on this host exceeds every difference between arms.

== Records

`raw/trace_*.jsonl` (+ `.stderr` headers), `raw/time_*.jsonl`, `raw/time_*.uptime`, `summary.json`; graph source `plot.py` (reads `summary.json` only).

== Missing

No isolated host (L1/L2) and no paired bootstrap; wall times are not admissible as a comparison here. The cost of generating Macaulay rows is not in the deterministic counter.
