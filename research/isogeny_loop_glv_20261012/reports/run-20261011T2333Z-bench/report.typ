#set page(paper: "a4", margin: 1.8cm)
#set text(font: "DejaVu Sans", size: 9pt)
#set par(justify: false)
#set table(stroke: 0.4pt + rgb("#c3c2b7"), inset: 4pt)
#show table: set text(7.6pt)
#text(13pt, weight: "bold")[Run report: run-20261011T2333Z-bench]

*What ran.* Five invocations of `gls_loop_bench`, each pinned to core 3 with `taskset -c 3`:
- D = −619 (`bench_run1`, start 23:33Z);
- D = −23, −71, −431, −3019 (`bench_D*`, start 23:34Z).

Each used `--pbits 127 --seed 20261012 --pairs 256 --rounds 30 --w 5`; the exact command is in each JSON. Source: the working tree that became checkpoint commit 54a16f6d4. rustc 1.97.0, release. Host: Intel Xeon \@ 2.10 GHz cloud container, not isolated. Load averages before and after each run are in `../../bench_*_env.txt` (0.02–0.39 at start).

*Arms.* wNAF (w = 5); GLS-2; GLV-GLS-4 with an affine loop; GLV-GLS-4 with a projective loop; wNAF again (A/A control). Stage diagnostics: decomposition, α′(P) affine, α′(P) projective. Unit: ns per scalar multiplication, the per-round mean over 256 pairs. Median and min of 30 rounds; all per-round values are in the JSON.

*Checks before timing (all passed, every run):*
- all arms agree on every pair;
- ψ(P) = [μ]P and α′(P) = [λ]P;
- each decomposition ≡ k (mod r);
- the baseline matches double-and-add on 8 pairs.

*Runs requested / completed.* 5 / 5. An earlier unpinned smoke run, with only a few rounds while a subagent was compiling, read A/A 0.916. It was discarded as noisy before these runs and is not recorded.

#table(columns: 6,
  [*D*], [*wNAF*], [*GLS-2*], [*GLV-GLS-4 aff.*], [*GLV-GLS-4 proj.*], [*A/A*],
  [−23], [141,473], [82,216], [71,168], [67,744], [138,830],
  [−71], [142,367], [79,139], [75,411], [70,640], [140,553],
  [−431], [140,578], [79,960], [78,460], [71,722], [138,122],
  [−619], [142,773], [81,665], [73,066], [71,172], [144,657],
  [−3019], [136,904], [80,870], [77,510], [76,152], [137,419],
)
Medians, ns. Derived ratios of medians (GLS-2 / GLV-GLS-4, projective loop): 1.214, 1.120, 1.115, 1.147, 1.062. Baseline / A/A ratios: 1.019, 1.013, 1.018, 0.987, 0.996.

#figure(image("../../figures/bench_arms.svg", width: 70%), caption: [Per-round values, D = −619, 30 rounds × 256 scalars (run-20261011T2333Z-bench, bench_run1; a later confirmation run appears beside it if recorded).])
#figure(image("../../figures/bench_trend.svg", width: 80%), caption: [Ratio of recorded medians against loop model cost, with round 2's derived model.])

*Observations only.* Interpretation: `../../README.md` §3.
