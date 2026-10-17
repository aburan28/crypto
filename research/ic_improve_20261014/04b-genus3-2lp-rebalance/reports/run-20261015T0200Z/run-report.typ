#set page(paper: "a4", margin: 1.8cm, numbering: "1")
#set text(size: 9.5pt)
#show table: set text(size: 8pt)

#align(center)[
  #text(size: 13pt, weight: "bold")[Run report: run-20261015T0200Z (04b genus-3 IC-2LP rebalance)]
  #v(2pt)
  Descriptive record of what ran. Interpretation is in the study README and REPORT.pdf.
]

*Spec and commit.* Protocol: `research/ic_improve_20261014/04b-genus3-2lp-rebalance/PROTOCOL.md`, frozen at commit 9752513ea (sha256 633e69bb5f71…). κ choice frozen at 288e74293 (`tuning/kappa_choice.json`). Tool: `examples/ic_improve_04b_genus3_2lp.rs` (sha256 595e3efa…), built with `cargo build --release --example ic_improve_04b_genus3_2lp` (binary sha256 8077f4d8…), rustc 1.97.0. Analysis: `examples/ic_improve_04b_analyze.rs`.

*Host and isolation.* 4 vCPU Intel Xeon @ 2.80 GHz, 15 GiB, Linux 6.18; shared cloud container, *no isolation earned*. Four jobs in parallel under `timeout 3600` and `prlimit --as=3 GiB` (`run_jobs.sh`). Start 2026-10-15T01:55Z; test records finished 04:39Z (`results/env-run-20261015T0200Z.txt`).

*Arms, seeds, workloads.* Seed 20261016. Sizes $k = 10..20$, q = largest prime below $2^k$. 7 GHS genus-3 curves per size (`results/curves/`). Tuning: $k = 14..20$ × κ ∈ {0.5, 0.7, 1, 1.4, 2} × curves 5–7 × targets 1–2. Test: curves 1–4 × targets 1–3 at κ = 2 for every k (BL and A2 from the same relation phase) and at κ = 1.4 for k = 17–20 (A1; κ\* = 2 at 14–16 reuses the κ = 2 run). Rho: $k = 10..14$, curves 1–4, target 1.

#table(columns: 5, align: (left, right, right, right, left),
  [*job set*], [*requested*], [*completed*], [*failed*], [*records*],
  [curves], [11], [11], [0], [`results/curves/` (77 curves)],
  [tuning solves], [210], [210], [0], [`results/tuning/`],
  [test solves], [180], [180], [0], [`results/test/`],
  [rho], [20], [20], [0], [`results/rho/`],
  [orchestration error], [—], [—], [1], [`results/orchestration_errors.jsonl`: an empty job list launched once (no tool arguments, exit 2); not a test job],
)

*Measured values (unit: multiplications in $FF_q$ plus mod ℓ; 12 test instances per size; 95% bootstrap intervals in `analysis.json/tables.md`).* Every one of the 360 test logarithms (180 solves × 2 linear-algebra routes), every tuning logarithm and every rho logarithm satisfied $[log]P = Q$ and equalled the planted value. Peak RSS 119 MB; longest job 490 s. One test solve replayed after the run reproduced its record's counts exactly.

#table(columns: 6, align: right,
  [*k*], [*BL mean total*], [*A1 mean total*], [*A2 mean total*], [*BL rows / nnz per row*], [*rho mean S (n = 4)*],
  [10], [3.142e6], [—], [2.559e6], [213 / 4.92], [550.0],
  [12], [1.560e7], [—], [1.199e7], [522 / 6.26], [439.7],
  [14], [9.686e7], [9.686e7], [7.566e7], [1301 / 7.76], [571.4],
  [16], [6.731e8], [6.731e8], [5.575e8], [3261 / 9.55], [—],
  [18], [4.387e9], [4.301e9], [3.801e9], [8202 / 11.71], [—],
  [20], [3.049e10], [2.930e10], [2.763e10], [20653 / 14.38], [—],
)
Even sizes shown; all sizes are in tables.md. Derived (analysis tool, OLS over $k = 14..20$): BL total exponent 1.383, A1 1.373, A2 1.419; BL relation 1.332, linear algebra 1.447.

#figure(image("../../figures/total_cost.svg", width: 90%), caption: [Recorded mean total multiplications per size, 12 test instances per arm, run-20261015T0200Z. The dashed line is a $q^(4\/3)$ guide through BL at $2^14$.])
#figure(image("../../figures/paired_ratios.svg", width: 90%), caption: [Recorded paired ratios arm/BL (geometric mean, 95% bootstrap intervals), run-20261015T0200Z.])

*Records.* `results/curves/curves_k*.jsonl`, `results/tuning/tu_*.jsonl`, `results/test/te_*.jsonl`, `results/rho/rho_*.jsonl`; analysis.json and tables.md from `ic_improve_04b_analyze report`. Missing data: none.
