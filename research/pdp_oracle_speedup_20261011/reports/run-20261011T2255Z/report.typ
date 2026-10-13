#set page(margin: 1.8cm, numbering: "1")
#set text(size: 9.5pt)
#set heading(numbering: "1.")
#show raw: set text(size: 8pt)

#align(center)[
  #text(size: 15pt, weight: "bold")[PS-2 run report: scalar group-law kernel, gates and program] \
  #text(size: 10pt)[PDP oracle speedup program, `research/pdp_oracle_speedup_20261011/ps2_scalar/` --- timing stage of 2026-10-11, 22:55--23:18 UTC]
]

= What ran

`ps2_scalar/run_ps2.sh` (correctness and counts, then the timing stage), `examples/pdp_speedup`, one release binary for every arm, Conductor task T-4. Host: Intel Xeon 2.80 GHz KVM guest, 4 vCPU, AVX-512F/BW, PCLMULQDQ, no GPU; timed runs on CPU 3 through `tools/isolated_bench.py run --cpus 3` (a shared container, not the repository's receipt service). Arms: `baseline` (`decompose`), `baseline2` (the same, the A/A), `so-full`, `so-liftable` (matched control), `group-scalar`. Five rounds, every round runs every arm once; workloads are the first 16 refuted PS-0 targets (6 at $l = 12$). 125 timed runs requested, 125 completed; 6 marked contended (kept), 9 refused starts retried. Raw records: `ps2_scalar/raw/{times,iso}.jsonl`, `verify_arms.jsonl`, `counts.jsonl`, `programs.jsonl`.

= Observations

#figure(
  image("ns_per_target.svg", width: 92%),
  caption: [Wall nanoseconds per refuted target by arm and rung (source `ps2_scalar/raw/times.jsonl`; hollow marks are runs the tool marked contended; black bars are per-arm medians; log scale).],
)

#table(
  columns: (auto, auto, auto, auto, auto, auto),
  align: (right, right, right, right, right, right),
  table.header[*l*][*baseline ms/target*][*so-full / baseline*][*so-liftable / baseline*][*group-scalar ms/target*][*group-scalar / baseline*],
  [8], [10.84], [1.034], [0.261], [0.337], [0.0308],
  [9], [46.99], [1.026], [0.250], [1.137], [0.0242],
  [10], [204.2], [1.000], [0.212], [3.596], [0.0177],
  [11], [876.0], [1.010], [0.247], [15.12], [0.0173],
  [12], [3795], [1.018], [0.242], [60.06], [0.0158],
)

#figure(
  image("ns_per_pair.svg", width: 92%),
  caption: [Wall nanoseconds per own pair (an arm's pairs are its own: the baseline visits every abscissa of $V$, the liftable arms only the rational base points).],
)

A/A (baseline against baseline), per target, B/A medians: 1.044, 1.023, 1.004, 1.009, 0.997 at $l = 8..12$; run-to-run spread of the baseline 8.15%, 3.44%, 3.53%, 2.77%, 3.43%.

Verdict digests equal PS-0's on every rung for all four arms (20 records, 0 mismatches); the group kernel's full root set equals the corpus root set on all 400 targets. Counted carry-less multiplies per pair: group kernel 46.9, 44.6, 43.6, 42.8, 42.4 against the baseline's 306.6, 329.8, 353.4, 377.2, 401.1. Gates of the recorded program (batch of 8 tasks): 63,630, 77,092, 85,718, 103,030, 124,886 live gates at $l = 8..12$ (hashes in `ps2_scalar/README.md`).

= Missing or not done

The incremental-update arm (Gray-code coefficient update plus the baseline's gcd test) was not built; a bound from counted operations (quartic formation is 17.2--22.5% of the baseline's carry-less multiplies) is arithmetic, not a measured time. The `semaev_gray_pairs` library module named in the protocol was not added (the kernel lives in the example crate). No corpus target is degenerate, so the baseline fallback of the kernel was exercised only by PS-1's small-field tests. Wall figures are exploratory in the repository's isolation vocabulary (shared container). Nothing end to end; nothing for ECC2K-130.
