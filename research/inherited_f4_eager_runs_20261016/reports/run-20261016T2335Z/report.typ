#set page(paper: "a4", margin: 1.6cm)
#set text(size: 9pt)
= Inherited F4 eager rewrite through shift runs — run report
Run `run-20261016T2335Z` · binary `f6_ic_probe` built from `1750ec759` (local branch on `f25d5e520`, with the gate-first default) · 2026-10-16

== What ran
`scripts/run_arms.sh`, one binary, arms by environment only: *bmi2* (default on this host: eager rewrite, `pdep ∘ pext` per segment), *runs* (`KIC_F4_EAGER=runs`: eager rewrite through shift runs, what a host without BMI2 now runs), *lazy* (`KIC_F4_EAGER=0`: lazy word-wise rewrite), *scatter* (`KIC_F4_REWRITE=scatter`: lazy rewrite with the bit scatter, what a host without BMI2 ran before this change). Nine cells as in `research/f6_gate_first_x86_20261016`, `--engines f4,f6p --reps 3`, 5 rounds, arm order rotated. Host: Intel Xeon @ 2.10 GHz, 4 vCPU KVM, BMI2, shared, unpinned, L0, load 1.01–1.92. 180 arm runs, all exited 0. No aarch64 host was available: *scatter* and *runs* are the before and after of a non-BMI2 host, measured here on x86.

Stage unit = per-arm total `word_xors`; CPU = probe `cpu_ms` (three repetitions). Ratios derived as in the gate-first report (medians over rounds; bracket = worst–best round pairing).

== Exactness
All 80 targets have the same decision tuple in every arm and round. Stage units are identical across rounds; *runs* charges exactly *bmi2*'s units and *lazy* exactly *scatter*'s, on every cell and engine.

== Measured values
#table(columns: 9, align: (left, left, right, right, right, right, right, right, right),
  [cell], [engine], [units runs = bmi2], [units scatter = lazy], [CPU scatter ms], [CPU runs ms], [CPU bmi2 ms], [scatter/runs], [bmi2/runs],
  [`K1_17`], [f4], [18019804], [15109049], [588.7], [574.0], [449.9], [1.03 (0.74–1.23)], [0.78],
  [`K1_17`], [f6p], [10830979], [8865230], [193.9], [183.0], [148.9], [1.06 (0.78–1.68)], [0.81],
  [`K1_23`], [f4], [62213299], [52572523], [1279.1], [920.5], [805.9], [1.39 (1.18–1.76)], [0.88],
  [`K1_23`], [f6p], [38122485], [32069197], [454.4], [340.8], [333.1], [1.33 (0.96–1.70)], [0.98],
  [`K0_23`], [f4], [167017420], [148425559], [3745.3], [2723.4], [2447.7], [1.38 (1.03–1.57)], [0.90],
  [`K0_23`], [f6p], [81022008], [67544699], [1157.9], [945.9], [746.0], [1.22 (1.00–1.54)], [0.79],
  [`K0_31`], [f4], [290412310], [264887288], [4781.3], [3013.3], [2476.9], [1.59 (1.51–1.72)], [0.82],
  [`K0_31`], [f6p], [164736149], [146257576], [2065.2], [1354.0], [1197.2], [1.52 (1.38–1.67)], [0.88],
  [`K0_19`], [f4], [17250812], [15282173], [459.7], [358.9], [291.3], [1.28 (0.97–1.55)], [0.81],
  [`K0_19`], [f6p], [10226773], [8837211], [147.2], [110.1], [89.6], [1.34 (0.66–1.89)], [0.81],
  [`K1_19`], [f4], [30159678], [25782851], [827.6], [650.7], [580.1], [1.27 (1.00–1.54)], [0.89],
  [`K1_19`], [f6p], [19117122], [16026687], [282.9], [264.6], [214.1], [1.07 (0.85–1.72)], [0.81],
  [`K1_17d7`], [f4], [22448248], [19473310], [649.8], [611.2], [562.6], [1.06 (0.83–1.27)], [0.92],
  [`K1_17d7`], [f6p], [10867029], [8987079], [167.9], [170.4], [148.3], [0.99 (0.81–1.11)], [0.87],
  [`K1_29`], [f4], [185919270], [165390163], [3235.8], [2039.6], [1753.8], [1.59 (1.44–1.71)], [0.86],
  [`K1_29`], [f6p], [112675706], [98192767], [1497.1], [1007.0], [846.5], [1.49 (1.35–1.77)], [0.84],
  [`K0_37`], [f4], [652414757], [604262254], [8779.5], [5159.3], [4033.1], [1.70 (1.56–1.86)], [0.78],
  [`K0_37`], [f6p], [376487242], [342788148], [3921.7], [2527.9], [2160.3], [1.55 (1.43–1.76)], [0.85],
)

#figure(image("fig_runs.svg", width: 100%), caption: [CPU ratios per cell, run-20261016T2335Z, 5 rounds.])

== Observations
- Against the bit scatter (the previous non-BMI2 path), shift runs: median CPU 1.22–1.70× lower on the six larger cells (K1/2^23, K0/2^23, K0/2^31, K0/2^19, K1/2^29, K0/2^37) on both engines, with the worst round pairing ≥ 1.00 on K0/2^31, K1/2^29, K0/2^37 (both engines) and K1/2^23 f4, K0/2^23 (both); on K1/2^19 1.27× (f4) and 1.07× (f6p); level (0.99–1.06×, brackets spanning 1) on K1/2^17 d6 and d7.
- Against `pdep ∘ pext` on this host, shift runs are 2–28% slower (bmi2/runs 0.78–0.98); the BMI2 path stays the default where BMI2 exists.
- Not measured: an aarch64 host. The shift-run kernel is plain integer code, so the ratio there may differ.

== Records
`raw/r<round>-<cell>-<arm>.jsonl`, `raw/*.status`, `host.txt`, `summary.txt`, `table.tsv` (`../../summarize.py`), `fig_runs.svg|png` (`../../plot.py`).
