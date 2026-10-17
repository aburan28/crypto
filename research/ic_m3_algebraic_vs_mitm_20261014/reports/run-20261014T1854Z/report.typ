#set page(paper: "a4", margin: 1.8cm)
#set text(size: 9.5pt)
#show table: set text(size: 8.5pt)
= Run report: claim grid 2026-10-14T18:54Z (descriptive)

*Run.* `research/ic_m3_algebraic_vs_mitm_20261014/run_grid.sh`, i.e. `ic_algebraic_floor run` for $m = 3$ and $m = 4$, run concurrently.
*Code.* Instrument and protocol at commit `76e88621e`; runner at `ef4434914`.
*Host.* Cloud container, Intel Xeon @ 2.10 GHz, 4 vCPU, Linux 6.18.44, rustc 1.97.0. Isolation L0: shared, the two grids concurrent, one thread per solve.
*Time.* 18:54–20:26 UTC.
*Workload.* Seed 20261014. Pools of 2,048 ($m = 3$) and 1,024 ($m = 4$). 6 decomposable and 6 non-decomposable targets per algebraic (arm, cell). Per target: 120 s and 4 GiB. Early stop after 4 leading censored targets.

This report records observations only; the verdicts are in the study README and `REPORT.pdf`.

== Solves by verdict

Verdicts:
- `yes`: decomposed, witness checked.
- `no`: refuted.
- `none`: no verdict.
- `timeout`: killed at budget + 15 s.

There were no crashes and 0 disagreements with ground truth.

#table(columns: 6, align: (left, left, right, right, right, right),
  [*m*], [*arm*], [*yes*], [*no*], [*none*], [*timeout*],
  [3], [f4-chain], [36], [36], [0], [0],
  [3], [f6-ic], [36], [36], [0], [0],
  [3], [s4-hybrid-f10], [24], [18], [6], [4],
  [3], [s4-hybrid-f12], [23], [18], [7], [4],
  [3], [sat], [18], [9], [13], [0],
  [4], [f4-chain], [24], [24], [3 (unsupported)], [0],
  [4], [f6-ic], [24], [24], [3 (unsupported)], [0],
  [4], [s22-hybrid-g8], [12], [6], [10], [0],
  [4], [sat], [15], [11], [14], [0],
)

== Early stops, skips and unsupported cells

*$m = 3$.*
- `sat` stopped at $n = 23$ and was skipped at 29 and 31.
- Both S#sub[4] hybrids stopped at $n = 29$ and were skipped at 31.

*$m = 4$.*
- `s22-hybrid-g8` stopped at $n = 17$.
- `sat` stopped at $n = 19$.
- F4 and F6-IC were unsupported at $n = 23, 29, 31$ (64-variable cap).

The floor equals exhaustive enumeration on 12/12 targets on every cell with $n <= 23$.

== Recorded values

#figure(image("cost_per_relation_m3.svg", width: 95%), caption: [$m = 3$, $log_2$ wall ns per relation (recorded; `ANALYSIS.json`).])
#figure(image("cost_per_relation_m4.svg", width: 95%), caption: [$m = 4$; shaded: only the floor ran.])

Raw records:
- `runs/m3/` and `runs/m4/`: `truth-*.json` per cell, and `arm-*.jsonl` per (arm, cell).
- `runs-m3.log` and `runs-m4.log`.
