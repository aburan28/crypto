#set page(paper: "a4", margin: 1.8cm)
#set text(size: 9.5pt)
#show table: set text(size: 8.5pt)
= Run report: amendment 1 grid, 2026-10-15T00:24Z (descriptive)

*Run.* `run_grid_amendment1.sh`: `ic_algebraic_floor run --m 4 --arms wide-chain` on seven cells, under `AMENDMENT_1.md` (frozen in commit `882598e55`).
*Host.* Cloud container, Intel Xeon @ 2.10 GHz, 4 vCPU. Isolation L0: the only grid running, one thread per solve.
*Time.* 00:24–01:05 UTC.
*Workload.* Seed 20261014 and a 1,024-target pool (the same targets as `runs/m4/`). 6 decomposable and 6 non-decomposable targets per cell. Per target: 120 s and 4 GiB.

This report records observations only; the verdict is in the study README and `REPORT.pdf`.

== Solves

- *Totals over 84 solves:* 40 `yes` and 30 `no`, every one agreeing with ground truth; 14 `timeout`, killed at budget + 15 s.
- *Censoring:* every non-decomposable target at $n = 29$ and $31$ timed out, plus one decomposable target at each. At $n = 31$ a further decomposable target returned `yes` only after 130 s, past the budget, so the analysis counts it as censored: 8 of 12 there, 7 of 12 at $n = 29$.
- No crashes, unsupported cells or early stops.
- The floor equals exhaustive enumeration on 12/12 targets on every cell with $n <= 23$.

#table(columns: 5, align: right,
  [*n*], [*d*], [*p*], [*floor*], [*wide-chain*],
  [11], [3], [0.009], [18.09], [27.20 [26.45, 28.42]],
  [13], [4], [0.562], [14.99], [26.05 [25.81, 26.27]],
  [17], [5], [0.181], [16.96], [29.85 [29.62, 30.10]],
  [19], [5], [0.093], [18.37], [32.27 [31.99, 32.59]],
  [23], [6], [0.026], [21.06], [35.96 [35.45, 36.63]],
  [29], [8], [0.281], [21.68], [≥ 38.58],
  [31], [8], [0.082], [23.48], [≥ 40.52],
)
Values are $log_2$ wall ns per relation, with 95% bootstrap intervals (`amendment-1/ANALYSIS.json`).

#figure(image("cost_per_relation_m4.svg", width: 100%), caption: [$m = 4$ with the wide solver; the original run's arms are shown faded. Open markers with arrows are lower bounds.])

Raw records: `runs/m4-wide/` and `runs-m4-wide.log`.
