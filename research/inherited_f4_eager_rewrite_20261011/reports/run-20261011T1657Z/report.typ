#set page(paper: "a4", margin: 2cm)
#set text(size: 10pt)
#set table(stroke: 0.4pt)

= Run report: inherited-F4 eager rewrite, `f6_ic_probe`, run-20261011T1657Z

*Observations only.* Interpretation is in the study report (`../../REPORT.pdf`).

== What ran

- *Binary:* `examples/f6_ic_probe` built in release from commit `6bdf93f4b` (branch `groebner-f4f5f6-audit-20261010`; its `src/cryptanalysis/inherited_f4.rs` is byte-identical to the PR branch `inherited-f4-eager-rewrite`). One binary for both arms; the arm is the environment variable `KIC_F4_EAGER` (`0` = lazy rows, unset/`1` = eager rewrite).
- *Commands* (one per file in `raw/`):
  - `KIC_F4_EAGER=$e f6_ic_probe --a 0 --n 23 --dim 8 --targets 4 --engines f4,f6p --reps 3`
  - `KIC_F4_EAGER=$e f6_ic_probe --a 0 --n 31 --dim 9 --targets 1 --engines f4,f6p --reps 3`
- *Order:* three rounds; in each round lazy then eager, n = 23 then n = 31 (interleaved).
- *Host and isolation:* Claude Code cloud container, 4 vCPU Intel Xeon \@ 2.10 GHz (AVX-512, BMI2), shared host, no core pinning, no frequency control. Isolation level: none (L0). Single-threaded probe.
- *Workload:* Koblitz curve $K_0$ over $FF_(2^n)$, standard-subspace factor base, $m = 3$, targets $[k]G$, $k = 1 dots$, node budget 20 000. Arms within the probe: `f4` (inherited F4) and `f6p` (F6-IC with the shared pair-sum index).
- *Requested / completed:* 12 probe invocations requested, 12 completed, no failures.
- *Unit:* `total_ms` is the probe's own wall time per arm: for each target the best of 3 repetitions, summed over targets (milliseconds, measured).

== Recorded values

#figure(image("totals.svg", width: 100%), caption: [Probe `total_ms` per arm and mode, run-20261011T1657Z, 3 rounds each (dots); bars are the median of the 3 recorded values.])

#table(columns: 6,
  [cell], [arm], [mode], [round 1], [round 2], [round 3],
  [n23 dim 8], [f4], [lazy], [1472.2], [1521.1], [1477.0],
  [n23 dim 8], [f4], [eager], [1058.4], [1068.6], [1022.3],
  [n23 dim 8], [f6p], [lazy], [998.2], [1034.8], [1046.8],
  [n23 dim 8], [f6p], [eager], [631.8], [776.6], [639.5],
  [n31 dim 9], [f4], [lazy], [1543.5], [1547.8], [1518.7],
  [n31 dim 9], [f4], [eager], [1057.2], [1036.3], [1041.3],
  [n31 dim 9], [f6p], [lazy], [1433.0], [1443.6], [1479.9],
  [n31 dim 9], [f6p], [eager], [946.8], [946.1], [911.6],
)
All values measured, ms. Derived: medians of the three rounds (shown above the bars) are n23 f4 1477.0 / 1058.4, n23 f6p 1034.8 / 639.5, n31 f4 1543.5 / 1041.3, n31 f6p 1443.6 / 946.1 (lazy / eager).

== Counters (identical in every round)

#table(columns: 6,
  [cell], [arm], [reductions], [splits], [word_xors lazy], [word_xors eager],
  [n23 dim 8], [f4], [2796], [1192], [158 512 775], [179 177 701],
  [n23 dim 8], [f6p], [1045], [636], [124 949 801], [144 396 510],
  [n31 dim 9], [f4], [1410], [602], [281 583 026], [310 429 172],
  [n31 dim 9], [f6p], [528], [323], [278 237 634], [305 397 167],
)
Reductions and splits are equal between modes; every target's verdict (`hit`, `exhausted`) is equal across modes and rounds (10 arm–target pairs, all `hit: false`, `exhausted: false`). `word_xors` differs by mode because the eager arm charges the rewrite of every kept row.

== Records

`raw/n23_dim8_eager{0,1}.round{1,2,3}.jsonl`, `raw/n31_dim9_eager{0,1}.round{1,2,3}.jsonl` (probe stdout, one JSON line per arm and target plus totals and the solver profile). Graph source: `plot_totals.py` (reads `raw/` only).

== Missing

No repeated sessions on an isolated host; no interval beyond the three rounds shown. The n23 `f6p` eager round 2 value (776.6 ms) is recorded as it came; no rerun.
