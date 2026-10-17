#set page(paper: "a4", margin: 1.8cm)
#set text(size: 9.5pt)
#set table(stroke: 0.4pt)

= Run report: `f4_fp` before and after, against msolve and Singular, run-20261011T2359Z

*Observations only.* Interpretation: `../../REPORT.pdf`.

== What ran

- *Runner:* `research/f4_fp_msolve_parity_20261011/run.sh <export_main> <export_new> <this dir> 3 300`.
- *Arms:* release builds of `examples/gb_cas_export.rs` with `src/cryptanalysis/f4_fp.rs` as on `main` (`4b6e09b78`, "main") and as on branch `f4-fp-msolve-gap` (`cba784a41`, "new"); binary SHA-256 in `raw/versions.txt`. Each at `RAYON_NUM_THREADS` 1 and 4. msolve 0.6.5 `-g 2` (reduced basis) at `-t 1` and `-t 4`, timed as process wall time (start-up included, about 6–10 ms). Singular 4.3.2 `std` under `dp`/`redSB`, timed inside Singular, 300 s limit.
- *Workload:* the 13 systems of `examples/f4_fp_bench.rs` (tier `fp`) and quad_n10, quad_n11, cubic_n6 (tier `fp-xl`): random dense systems, fixed seeds, grevlex.
- *Repetitions:* 3 interleaved rounds. Inside a round the `f4_fp` arms report the median of 3 repetitions (`fp`) or 1 (`fp-xl`).
- *Not run:* the main arm on quad_n11. On this host its quad_n10 takes 131 s single-threaded, so quad_n11 would take many minutes per repetition. Singular on quad_n11 hit the 300 s limit in all 3 rounds.
- *Host:* Claude Code cloud container, `Intel Xeon @ 2.80GHz`, 4 vCPU, shared, no pinning (isolation L0). An earlier attempt on 2026-10-11 (22:37Z) was lost to a container restart and is not used.

== Medians (ms) and agreement

#figure(image("fp_vs_msolve.svg", width: 100%), caption: [Larger systems, run-20261011T2359Z; bars at the median of 3 rounds, dots per round.])

#table(columns: 8,
  [case], [main 1], [new 1], [msolve 1], [main 4], [new 4], [msolve 4], [Singular],
  [quad_n8], [1195.2], [80.2], [78], [780.3], [82.8], [82], [1005],
  [quad_n9], [10463.8], [351.3], [429], [5332.4], [282.2], [235], [8165],
  [cubic_n5], [445.3], [27.9], [29], [353.8], [27.5], [42], [89],
  [overdet_n10_m15], [62.2], [39.9], [50], [45.5], [38.3], [50], [450],
  [overdet_n12_m24], [91.3], [74.1], [105], [79.0], [57.0], [75], [588],
  [quad_n10], [130759.5], [2417.7], [2417], [51791.8], [1285.6], [1125], [68360],
  [quad_n11], [not run], [17795.2], [17466], [not run], [6409.1], [6544], [> 300 s],
  [cubic_n6], [13864.8], [295.1], [283], [6510.8], [217.2], [242], [2778],
)
All values measured; the remaining small systems are in `medians.json`. Basis fingerprints are identical between main and new on every system. In every round Singular found the new build's basis to be its own reduced Gröbner basis (45 of 45 checks). On quad_n11, where Singular timed out, Singular compared it with msolve's reduced basis: 767 elements each, each reducing to zero against the other (`raw/quad_n11_basis_check.txt`).

#figure(image("ratios.svg", width: 100%), caption: [Ratios of the medians above. msolve's ratios on the smallest systems mostly measure its process start-up.])

== Records

`raw/ours_{main,new}_t{1,4}_{fp,fp-xl}.round{1,2,3}.jsonl`, `raw/external.round{1,2,3}.tsv`, `raw/versions.txt`, `raw/quad_n11_basis_check.txt`. The exported inputs and msolve's bases (83 MB) are not kept: `manifest.sha256` lists their hashes, and `gb_cas_export` and msolve regenerate them. `medians.json` and the graphs come from `../../plot.py`, which reads `raw/` only.
