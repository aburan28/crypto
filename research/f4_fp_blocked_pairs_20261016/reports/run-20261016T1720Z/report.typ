#set page(paper: "a4", margin: 1.8cm)
#set text(size: 9.5pt)
= f4_fp: blocked Faugère–Lachartre pair reduction — run report
Run `run-20261016T1720Z` · code commit `1c55f48bb` (branch `claude/f4-f5-algorithm-speed-bvp7qd`) · 2026-10-16

== What ran
`examples/gb_cas_export` (release) on the 13 standard and 5 large (`GB_CAS_LARGE=1`) F#sub[p] systems, one case per process (`GB_CAS_ONLY`), `RAYON_NUM_THREADS=1`, one F4 run per process. Arms, interleaved per case in each of 3 rounds (`raw/run_ab.sh`):
- *blocked*: default, `F4_FP_FL_BLOCKED=1`;
- *control*: same binary, `F4_FP_FL_BLOCKED=0` (row-by-row `D − C·B'`, the code before this change);
- *msolve 0.6.5* (Ubuntu package), `msolve -g 2 -t 1` on the exported `.ms` file of the same system; process wall time, so it includes msolve's start-up (≈ 5–7 ms, visible on the small cases).
Ours is the in-process F4 time reported by the exporter. Host: Intel Xeon @ 2.10 GHz, 4 vCPU KVM cloud container, shared, unpinned (L0). Runs requested 162, completed 162, no failures.

== Correctness
The reduced basis fingerprint (`basis_fp`) is identical between blocked and control in every round for all 18 cases (`raw/ab_full.tsv`). `f4_fp` unit tests 28/28 pass, including the new `blocked_pair_reduction_matches_direct_product`.

== Measured values (median of 3 rounds, ms; ratios derived as median / median)
#table(columns: 6, align: (left, right, right, right, right, right),
  [case], [blocked], [control], [msolve], [blocked/control], [blocked/msolve],
  [`quad_n4_p31`], [1.1], [0.9], [6.8], [1.199], [0.163],
  [`quad_n5_p31`], [1.9], [2.1], [5.6], [0.911], [0.334],
  [`quad_n6_p65521`], [4.1], [5.0], [6.7], [0.817], [0.618],
  [`quad_n7_p65521`], [15.1], [16.8], [17.0], [0.898], [0.888],
  [`cubic_n3_p29`], [1.4], [1.3], [6.1], [1.065], [0.231],
  [`cubic_n4_p31`], [3.3], [4.0], [7.2], [0.825], [0.459],
  [`quart_n3_p31`], [2.2], [2.4], [5.9], [0.900], [0.372],
  [`overdet_n6_m9_p101`], [1.9], [1.5], [7.0], [1.238], [0.267],
  [`quad_n8_p65521`], [55.0], [53.2], [47.6], [1.033], [1.156],
  [`quad_n9_p65521`], [255.9], [289.8], [230.1], [0.883], [1.112],
  [`cubic_n5_p31`], [23.0], [25.1], [19.1], [0.917], [1.204],
  [`overdet_n10_m15_p65521`], [31.4], [27.4], [30.0], [1.147], [1.046],
  [`overdet_n12_m24_p65521`], [58.6], [57.6], [62.8], [1.018], [0.933],
  [`quad_n10_p65521`], [1410.5], [1803.0], [1798.2], [0.782], [0.784],
  [`quad_n11_p65521`], [9370.3], [10548.0], [12425.5], [0.888], [0.754],
  [`cubic_n6_p65521`], [261.3], [258.4], [175.0], [1.011], [1.493],
  [`quart_n4_p65521`], [15.4], [14.8], [14.7], [1.038], [1.047],
  [`overdet_n14_m28_p65521`], [238.6], [262.5], [306.0], [0.909], [0.780],
)
Per-round ranges (`table.tsv`): quad_n10 blocked 1308–1490 vs control 1752–1810; quad_n11 9299–9667 vs 10079–10820; quad_n9 241–294 vs 266–315 (overlap). Cases under ≈ 30 ms vary in both directions by ±20 % between arms and are within this host's noise.

#figure(image("fig_ratio.svg", width: 100%), caption: [Median ratios per case, 1 thread, n = 3 rounds per arm, run-20261016T1720Z.])
#figure(image("fig_times.svg", width: 90%), caption: [Median wall time for the cases where msolve takes ≥ 20 ms, run-20261016T1720Z.])

== Observations
- quad_n10 and quad_n11: blocked 0.78× and 0.89× control, per-round ranges disjoint. quad_n9 0.88×, ranges overlap.
- cubic_n6, quart_n4, quad_n8, overdet cases: within noise (0.91–1.15×).
- Against msolve on this host in this run: blocked 0.75–0.78× on quad_n10/n11, 1.11–1.20× on quad_n8/n9/cubic_n5, 1.49× on cubic_n6 (msolve's cubic_n6 ranged 153–204 ms across rounds).
- Not measured here: 4 threads. Exploratory runs before this report (not recorded as a run) showed no difference between arms at 4 threads (quad_n10 ≈ 0.91 s both).
- An extension of the blocking to the `B'` build was also tried exploratorily and showed no gain over this change; it is not in the commit.

== Records
`raw/ab_full.tsv` (round, case, arm, ms, basis_fp), `raw/run_ab.sh`, `table.tsv`, `fig_*.svg|png`; plotting source `../../plot.py`.
