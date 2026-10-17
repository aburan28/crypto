#set page(paper: "a4", margin: 1.7cm)
#set text(size: 9.3pt)
#set table(stroke: 0.4pt)
#show table: set text(size: 8.3pt)
= Run report: $C_4$ ladder and ablation, 2026-10-14T20:08Z

Command: `research/gaudry_quartic_d4_cost_20261014/run_ladder.sh` (`RAYON_NUM_THREADS=3`), binary `examples/gaudry_quartic_d4` sha256 `c818a579…f524`, commit `2b49f38d1`, protocol `PROTOCOL.md` (Part A). Host: Xeon \@ 2.10 GHz, 4 vCPU, shared with the ecbench session below. Wall times are descriptive and admitted nowhere.

22 configurations ran. Each wrote `raw/<name>.jsonl` (one row per closed solve) and `raw/<name>.log`; none failed. Unsolved systems are counted from the logs ("did not close"): 20 in all, every one listed in `raw/summary.json` per configuration. All 648 closed solves agree with the meet-in-the-middle enumeration (measured). Means and intervals are derived by `examples/gaudry_quartic_d4_summary.rs`: 10,000 bootstrap resamples, seed 1.

#table(columns: 5, align: (left, right, right, right, right),
  [configuration], [closed / unsolved], [mean $C_4$], [95 % interval], [first solve],
  [p29 new planted], [17 / 3], [$5.49 dot 10^7$], [5.41–5.62], [$1.19 dot 10^8$],
  [p29 new uniform], [297 / 3], [$6.23 dot 10^7$], [6.15–6.34], [$1.19 dot 10^8$],
  [p29 committed planted], [17 / 3], [$2.35 dot 10^10$], [2.30–2.39], [],
  [p29 committed uniform], [50 / 0], [$2.36 dot 10^10$], [2.34–2.37], [],
  [p269 new planted], [19 / 1], [$6.90 dot 10^7$], [6.88–6.92], [$1.41 dot 10^8$],
  [p269 new uniform], [50 / 0], [$6.84 dot 10^7$], [6.83–6.85], [$1.40 dot 10^8$],
  [p269 committed planted], [9 / 1], [$2.55 dot 10^10$], [2.54–2.55], [],
  [p1009 new planted / uniform], [20 / 0, 50 / 0], [$8.21$ / $8.39 dot 10^7$], [], [$1.55 dot 10^8$],
  [p1009 committed planted], [5 / 0], [$2.56 dot 10^10$], [], [],
  [p4049 new planted / uniform], [10 / 0, 20 / 0], [$8.06$ / $8.47 dot 10^7$], [], [$1.59 dot 10^8$],
  [p4049 committed planted], [3 / 0], [$2.56 dot 10^10$], [], [],
)

The nine-configuration ablation at $p = 269$ is in the study's `README.md` and `figures/c4_ablation.svg`.

#figure(image("c4_ladder.svg", width: 72%), caption: [$C_4$ against $p$ (recorded means and bootstrap intervals).])
