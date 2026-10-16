#set page(paper: "a4", margin: (x: 1.6cm, y: 1.5cm))
#set text(size: 9pt)
#set par(justify: true, leading: 0.5em)

= Run report: `20261014T1800Z` (odd-characteristic GHS genus enumeration)

This report is descriptive: it records what ran and what the tool output.
Interpretation is in `../../REPORT.md`.

#table(columns: 2, inset: 3pt, stroke: 0.4pt + rgb("#bbbbbb"),
  [command], [`rustc --edition 2021 -O examples/odd_char_ghs_genus.rs -o odd_char_ghs_genus && ./odd_char_ghs_genus results/run-20261014T1800Z`],
  [commit], [`78a21243a` (tool sha256 `1e589153…a8922`, protocol sha256 `11e97a1f…def533`)],
  [host], [Intel Xeon @ 2.80 GHz, 4 logical CPUs, Linux 6.18 x86_64, rustc 1.97.0; single thread],
  [isolation], [none; this is an exact computation and no timing is claimed (wall 0.67 s, informational)],
  [inputs], [$n = 1..8$ (types and counts), $n <= 14$ (minimum genus), $q = p, p^2$ with $p in {5,...,31}$, and 81 Joux–Vitse instances ($p <= 31$, 3 $alpha$, 3 arms)],
  [exit / failures], [0 / 0 failed checks],
)

*Recorded checks* (all from the tool's own output files in `results/run-20261014T1800Z/`):
- 601 types for $n <= 8$. Riemann–Hurwitz equals the Kani–Rosen sum on every type.
- `counts.csv`: 144 of 144 `(q, n)` rows have `check_total = ok` (types add up to $q^(4n) - q^(2n)$) and `check_sigma_stable = ok`.
- `min_genus.csv`: `check_total_q5 = ok` for $n = 1..14$.
- `jv_family.csv`: on 81 of 81 rows,
  - `L_check = ok`, `jac_eq_numE = true` and `lemma1_match = true`;
  - the cycle types are `1 1 3 3` (54 rows, linear and split arms) and `2 6` (27 rows, conjugate arm).
- Unit tests (`rustc --test`): 11 of 11 passed.

*Run history.* Run `20261014T1758Z` (commit `4ad7dd5db`) was kept unchanged in its own directory. It lacks only the `check_total_q5` column, which the protocol amendment added. Its `types.csv`, `counts.csv`, `genus_by_n.csv` and `jv_family.csv` are byte-identical to this run's (checked with `cmp`). No run failed or was interrupted.

#figure(image("../../figures/genus_vs_n.svg", width: 100%),
  caption: [Recorded values from this run's `min_genus.csv`, `genus_by_n.csv` and `counts.csv`, drawn by `figures/plot.py`.])
