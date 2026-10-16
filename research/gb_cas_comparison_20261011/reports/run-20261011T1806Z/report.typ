#set page(paper: "a4", margin: 1.8cm)
#set text(size: 9.5pt)
#set table(stroke: 0.4pt)

= Run report: Gröbner engines against msolve, Singular and PolyBoRi, run-20261011T1806Z

*Observations only.* Interpretation: `../../REPORT.pdf`.

== What ran

- *Runner:* `research/gb_cas_comparison_20261011/run.sh <bin> <this dir> 3 600 120`, with `POLYBORI_PY` set. Binaries: release builds of `examples/gb_cas_export.rs` and `examples/f6_ic_probe.rs` from commit `6bdf93f4b` plus the new exporter (the engine sources equal PR \#1671's).
- *External tools* (`raw/versions.txt`): msolve 0.6.5-1build2 and Singular 4.3.2-p10 (Ubuntu packages); PolyBoRi through `passagemath-brial` 10.8.12 (Python 3.13). Magma: *not run*. It is not installed and needs a licence; its scripts are in `inputs-fp/*.magma` and `inputs-bool-n17/*.magma`.
- *Host:* Claude Code cloud container, Intel Xeon \@ 2.10 GHz, 4 vCPU, shared, no pinning (isolation L0).
- *F_p workload:* the 13 systems of `examples/f4_fp_bench.rs` (same generator and seeds; basis fingerprints equal the bench's). Ours: `f4_fp::f4`, grevlex, the bench's degree bound, median of 3 repetitions per round, with `RAYON_NUM_THREADS` 1 and 4. msolve: `-g 2` (reduced GB), `-t 1` and `-t 4`, wall time of the process (and its own `overall(elapsed)`, 10 ms resolution). Singular: `std` under `dp` with `redSB`, timed by `rtimer` inside Singular. Three rounds, interleaved.
- *Boolean workload:* the 8 root decomposition systems of `f6_ic_probe --a 1 --n 17 --dim 6 --targets 8` (Koblitz $K_1$ over $FF_(2^17)$, $m = 3$, a 63-point factor base; 35 Boolean variables and 34 cubics each, from `KIC_F4_NODE_DUMP`). Ours: the probe's per-target decision time, best of 3 repetitions, for inherited F4 (`f4`) and F6-IC with the pair index (`f6p`). External: one full reduced Gröbner basis with the field equations (msolve `-t 1`, Singular `std`), or in the Boolean ring (PolyBoRi `groebner_basis`). Time limits: 600 s (msolve, PolyBoRi) and 120 s (Singular).
- *Requested and completed:* all F_p runs completed. Boolean: msolve and PolyBoRi completed 8 of 8 each; Singular completed 0 of 8 (all hit the 120 s limit, exit 124).

== F_p: reduced Gröbner bases (median ms of 3 rounds)

#figure(image("fp.svg", width: 100%), caption: [F_p cases with at least 10 ms of work, run-20261011T1806Z.])

#table(columns: 7,
  [case], [ours 1 thr], [ours 4 thr], [msolve 1 thr], [msolve 4 thr], [Singular], [same basis],
  [quad_n4_p31], [0.6], [0.6], [7], [8], [0], [yes],
  [quad_n5_p31], [2.8], [2.9], [7], [8], [1], [yes],
  [quad_n6_p65521], [17.0], [13.7], [8], [10], [9], [yes],
  [quad_n7_p65521], [92.3], [71.1], [14], [17], [68], [yes],
  [cubic_n3_p29], [1.4], [1.2], [6], [9], [0], [yes],
  [cubic_n4_p31], [18.8], [14.4], [9], [10], [3], [yes],
  [quart_n3_p31], [8.0], [6.4], [7], [10], [1], [yes],
  [overdet_n6_m9_p101], [1.1], [1.2], [7], [10], [1], [yes],
  [quad_n8_p65521], [694.7], [534.1], [57], [44], [567], [yes],
  [quad_n9_p65521], [5693.2], [3479.0], [246], [145], [4106], [yes],
  [cubic_n5_p31], [294.2], [234.8], [22], [22], [64], [yes],
  [overdet_n10_m15_p65521], [32.7], [26.4], [34], [35], [223], [yes],
  [overdet_n12_m24_p65521], [63.3], [48.3], [61], [48], [331], [yes],
)
All measured. "Same basis": in every round, Singular found our basis to be a reduced Gröbner basis of the same ideal as its own (`same_ideal 1`, `same_size 1`: equal size, mutual normal forms zero, equal leading terms). msolve's wall time includes process start-up and parsing (about 6–7 ms at the smallest cases).

== Boolean: deciding the n = 17 decomposition systems (ms)

#figure(image("bool_n17.svg", width: 100%), caption: [Per-target decision time, run-20261011T1806Z. Singular's points are the 120 s limit, not a result.])

#table(columns: 7,
  [target], [ours f4], [ours f6p], [msolve], [PolyBoRi], [Singular], [verdict (all that finished)],
  [1], [24.1], [16.3], [42 741], [40 421], [> 120 000], [unit ideal],
  [2], [21.6], [14.2], [44 061], [42 448], [> 120 000], [unit ideal],
  [3], [11.2], [7.2], [45 517], [52 296], [> 120 000], [decomposes],
  [4], [22.3], [13.9], [45 459], [42 692], [> 120 000], [unit ideal],
  [5], [23.2], [13.6], [45 845], [41 753], [> 120 000], [unit ideal],
  [6], [21.7], [13.3], [46 527], [43 659], [> 120 000], [unit ideal],
  [7], [21.1], [13.3], [44 485], [41 012], [> 120 000], [unit ideal],
  [8], [24.3], [15.2], [44 436], [41 005], [> 120 000], [unit ideal],
)
Verdicts agree: ours (`hit`), msolve (`[1]` or not) and PolyBoRi (`unit_ideal`) give the same answer on all 8 targets. Ours decides by splitting variables and degree-bounded inherited Macaulay bases; the external tools compute the full reduced basis of the ideal. The two are different computations that answer the same question.

== Records

`raw/fp_ours_t{1,4}.round{1,2,3}.jsonl`, `raw/fp_external.round{1,2,3}.tsv`, `raw/*.msolve_t*.out.gz` (msolve's bases), `raw/bool_n17_ours.jsonl`, `raw/bool_n17_external.tsv`, `raw/bool_n17_roots.jsonl`, `raw/versions.txt`; inputs in `inputs-fp/` and `inputs-bool-n17/` (`.ms` msolve, `.sing` Singular, `.magma` Magma). Graph source: `plot.py`, which reads `raw/` only.

== Missing

Magma timings: not run (no Magma here). No isolated host and no interval beyond the three F_p rounds; one external run per Boolean root. The n = 23 decomposition systems (47 variables) were not run externally.
