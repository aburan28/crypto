#set page(paper: "a4", margin: 1.7cm)
#set text(size: 9.5pt)
= Gate-first F6-IC on x86 — run report
Run `run-20261016T2232Z` · binary `f6_ic_probe` built from `f25d5e520` (origin/main) · 2026-10-16

== What ran
`scripts/run_arms.sh`: `f6_ic_probe --engines f4,f6p --reps 3` (`f4` inherited F4, `f6p` F6-IC with the shared pair-sum index), one binary, arms set by environment only — *off* (no switch), *gate* (`KIC_F6_GATE_FIRST=1`), *rref* (`KIC_F4_INHERIT_ROOT=rref`), *gaterref* (both) — on nine cells: the four of `research/lit2026_w5_lane_batched_f4_20261011` (K1/2^17 d6, K1/2^23 d7, K0/2^23 d8, K0/2^31 d9) and five new ones (K0/2^19 d7, K1/2^19 d6, K1/2^17 d7, K1/2^29 d8, K0/2^37 d9), m = 3, `max_degree` 3. 5 rounds, arm order rotated each round. Host: Intel Xeon @ 2.10 GHz, 4 vCPU KVM container, BMI2 present, so the inherited-F4 eager rewrite (`KIC_F4_EAGER`, default) is in use — the path the aarch64 study did not measure. Shared, unpinned, L0; load (1-min) 1.02–4.21. A repository search ran on one core during round 1 (22:33–22:36 UTC). 180 arm runs requested, 180 exited 0.

Stage unit = the probe's per-arm total `word_xors` (one pass, deterministic). CPU = the probe's `cpu_ms` (process user+system over its three repetitions). Ratios are derived: off / arm of the stage units, and of the medians over the 5 rounds for CPU; the bracket is (min off / max arm) – (max off / min arm) over rounds.

== Exactness
Every one of the 80 targets (cell × engine × target) has the same decision tuple — hit, indices, exhausted, reductions, splits, propagations, infeasible, refutations, witnesses, partial refutations, geometric additions — in every arm and every round. Every arm's stage unit is identical in all 5 rounds.

== F6-IC (`f6p`), gate-first against off
#table(columns: 7, align: (left, right, right, right, right, right, right),
  [cell], [units off], [units gate], [off/gate], [CPU off ms], [CPU gate ms], [CPU off/gate],
  [`K1_17`], [15514870], [10830979], [1.433], [258.7], [163.1], [1.59 (1.37–2.25)],
  [`K1_23`], [59881009], [38122485], [1.571], [679.1], [320.7], [2.12 (1.56–2.27)],
  [`K0_23`], [132377695], [81022008], [1.634], [1461.7], [752.8], [1.94 (1.72–2.41)],
  [`K0_31`], [285541596], [164736149], [1.733], [2267.1], [1156.9], [1.96 (1.85–2.20)],
  [`K0_19`], [13660081], [10226773], [1.336], [181.8], [116.5], [1.56 (1.25–1.75)],
  [`K1_19`], [28309592], [19117122], [1.481], [400.2], [206.5], [1.94 (1.59–2.53)],
  [`K1_17d7`], [15297452], [10867029], [1.408], [278.0], [132.3], [2.10 (1.69–2.63)],
  [`K1_29`], [182738266], [112675706], [1.622], [1554.7], [872.4], [1.78 (1.69–2.06)],
  [`K0_37`], [644803622], [376487242], [1.713], [3947.2], [2097.3], [1.88 (1.60–2.07)],
)
Gate-first skipped 447 / 432 / 840 / 405 refuted nodes on the four W5 cells, the same counts the aarch64 study recorded.

== Inherited F4 (`f4`), gate-first against off
#table(columns: 7, align: (left, right, right, right, right, right, right),
  [cell], [units off], [units gate], [off/gate], [CPU off ms], [CPU gate ms], [CPU off/gate],
  [`K1_17`], [18019804], [18019804], [1.000], [474.0], [468.7], [1.01 (0.78–1.20)],
  [`K1_23`], [62213299], [62213299], [1.000], [776.8], [842.5], [0.92 (0.81–1.08)],
  [`K0_23`], [167017420], [167017420], [1.000], [2487.9], [2425.8], [1.03 (0.98–1.12)],
  [`K0_31`], [290412310], [290412310], [1.000], [2367.7], [2439.6], [0.97 (0.85–1.09)],
  [`K0_19`], [17250812], [17250812], [1.000], [336.5], [338.9], [0.99 (0.76–1.21)],
  [`K1_19`], [30159678], [30159678], [1.000], [587.6], [567.4], [1.04 (0.92–1.34)],
  [`K1_17d7`], [22448248], [22448248], [1.000], [579.1], [491.3], [1.18 (0.75–1.35)],
  [`K1_29`], [185919270], [185919270], [1.000], [1713.3], [1747.2], [0.98 (0.90–1.13)],
  [`K0_37`], [652414757], [652414757], [1.000], [4194.7], [4128.3], [1.02 (0.93–1.15)],
)

== F6-IC, the RREF root (`KIC_F4_INHERIT_ROOT=rref`)
#table(columns: 5, align: (left, right, right, right, right),
  [cell], [units off/rref], [units off/gaterref], [CPU off/rref], [CPU off/gaterref],
  [`K1_17`], [0.931], [1.275], [0.90], [1.73],
  [`K1_23`], [1.005], [1.488], [1.12], [2.33],
  [`K0_23`], [1.102], [1.660], [1.08], [2.20],
  [`K0_31`], [1.226], [1.934], [1.18], [2.12],
  [`K0_19`], [0.911], [1.138], [1.03], [1.62],
  [`K1_19`], [0.976], [1.324], [1.01], [1.77],
  [`K1_17d7`], [0.889], [1.180], [0.98], [1.86],
  [`K1_29`], [1.101], [1.637], [1.13], [1.88],
  [`K0_37`], [1.227], [1.911], [1.17], [2.04],
)

#figure(image("fig_gate.svg", width: 100%), caption: [Gate-first ratios per cell, run-20261016T2232Z, 5 rounds.])
#figure(image("fig_arms.svg", width: 100%), caption: [F6-IC stage-unit ratios of the three switched arms, run-20261016T2232Z.])

== Observations
- Gate-first, F6-IC: 1.336–1.733× fewer stage units on all nine cells; median CPU 1.56–2.12× lower; the worst round pairing is above 1 on every cell (lowest 1.25, K0/2^19).
- Gate-first, inherited F4: stage units identical (the engine has no node oracle); CPU medians 0.92–1.18×, within this host's noise.
- RREF root: on x86 it lowers stage units on the larger cells (1.10–1.23× on K0/2^23, K0/2^31, K1/2^29, K0/2^37) and raises them on the smaller ones (0.89–0.98× on K1/2^17 d6 and d7, K0/2^19, K1/2^19); K1/2^23 is level. The aarch64 study found it lower on all four of its cells.

== Records
`raw/r<round>-<cell>-<arm>.jsonl` (first line: arm, cell, round, UTC, load), `raw/*.status`, `host.txt`, `summary.txt` and `table.tsv` (from `../../summarize.py`), `fig_*.svg|png` (from `../../plot.py`).
