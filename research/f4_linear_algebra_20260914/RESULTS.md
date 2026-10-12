# F4 linear algebra results — engineering, 2026-09-14

The suffix-XOR kernel removes real elimination work while preserving the exact
RREF. Across 18 matrices, the geometric mean of paired kernel time ratios is
0.717 candidate/reference (1.39x reference/candidate).
This is a kernel diagnostic. The full native solver and toy DLP timings do not
establish an end-to-end improvement.

The primary metric remains unmeasured:

| Variant | Class | Total calibrated operations | S | Matched rho ratio | Floor ratio | Correctness |
| --- | --- | --- | --- | --- | --- | --- |
| Unmodified F4 | engineering baseline | null | null | null | null | reference |
| Suffix XOR F4 | engineering | null | null | null | null | all comparisons pass |

No change to the generic-group floor or exponent is claimed. Even eliminating
this phase entirely is bounded by `1/(1-f)`, where `f` is its complete-cost share.
Matched rho and all-phase operation conversion were not measured; nulls must not
be interpreted as zero. The numerical kernel acceptance target in the frozen
contract was met, not a total-attack-cost target.

## What changed

The baseline repeatedly XORed whole rows, including zero words before each
pivot. The candidate skips those words and uses disjoint pivot/target borrows
with contiguous slices. Pivot search and elimination order remain identical.
This affects standard F4, block-constrained F4 and Macaulay profiling. It adds
no allocation or dependency. The subgroup-modulus relation-matrix solver is
unchanged; GF(2) XOR is not interchangeable with that arithmetic.

The XOR counter measures actual elimination word operations in both versions.
It excludes pivot scans, swaps, construction, extraction and other solver work.
Reducing this counter is not a calibrated reduction in total ECDLP operations.

## Correctness and provenance

- Accepted evidence: `results/run-002/`, 229 recorded processes, no missing runs.
- 18 matrix inputs, 11 interleaved pairs each: 198 exact rank and full-row comparisons.
- 304 native S3 inputs, including the original 152 and 152 fresh systems, each run
  three times in both binaries: 1,824 complete assignment checks against exhaustive truth.
- 144 full-DLP executions: 72 verified n=9 executions and 72 incomplete n=7 attempts.
  Every paired outcome, scalar, input parameter and mathematical work count matches.
- All 15 Gröbner/RREF tests pass, including 600 shape/density/seed checks and a
  padding/skipped-column check. CI now includes these tests. The library and `ic`
  also pass `cargo check --no-default-features`.
- Reference production source: `466e8b651002ad11e5de996619f9aaaedb5ac551`.
  Both release builds use Rust 1.90.0, the same saved dependency lock and features,
  caches off, one Rayon thread and the same pinned CPU. The copied reference
  kernel matches the original apart from its name and whitespace.

Run 001 is explicitly rejected and retained. A shared Cargo artifact cache
initially produced identical solver binaries, and a linear synthetic PRNG
capped dense matrix rank at 64. Run 002 uses forced clean package builds and
SplitMix64 input generation. The comparator now rejects identical binaries and
unexpectedly low dense ranks. No production algorithm or acceptance threshold
was tuned between those runs. A later attempt to use the same output directory
was refused by the runner's no-overwrite guard; the complete saved run and all
source/binary hashes were checked before acceptance.

## Kernel measurements

All rows below pass exact rank and full-RREF equality. Seed 17 is development;
937 is a fresh holdout. Setup, copying and verification are separately recorded
in the raw output. Medians and paired ratios differ because the latter are
geometric means of pairwise ratios. Intervals are the frozen paired bootstrap;
they remain shared-host timing diagnostics.

| Matrix / seed | Shape; rank | Variant | 64-bit XORs | XOR C/R | Median kernel ms | Paired time C/R [95% CI] |
| --- | --- | --- | --- | --- | --- | --- |
| deficient-1024x1024 / 17 | 1024x1024; 512 | reference | 3,149,104 | 1.000 | 3.7413 | 1.000 |
| deficient-1024x1024 / 17 | 1024x1024; 512 | candidate | 2,546,720 | 0.809 | 2.7616 | 0.732 [0.637, 0.840] |
| deficient-1024x1024 / 937 | 1024x1024; 512 | reference | 3,144,896 | 1.000 | 3.7134 | 1.000 |
| deficient-1024x1024 / 937 | 1024x1024; 512 | candidate | 2,543,562 | 0.809 | 2.8809 | 0.790 [0.696, 0.920] |
| dense-1024x1024 / 17 | 1024x1024; 1024 | reference | 8,390,656 | 1.000 | 10.4124 | 1.000 |
| dense-1024x1024 / 17 | 1024x1024; 1024 | candidate | 4,453,746 | 0.531 | 7.0649 | 0.684 [0.657, 0.715] |
| dense-1024x1024 / 937 | 1024x1024; 1023 | reference | 8,365,616 | 1.000 | 10.2157 | 1.000 |
| dense-1024x1024 / 937 | 1024x1024; 1023 | candidate | 4,448,083 | 0.532 | 6.5838 | 0.643 [0.638, 0.648] |
| dense-1024x4096 / 17 | 1024x4096; 1024 | reference | 33,467,968 | 1.000 | 25.5255 | 1.000 |
| dense-1024x4096 / 17 | 1024x4096; 1024 | candidate | 29,546,673 | 0.883 | 12.7130 | 0.510 [0.495, 0.529] |
| dense-1024x4096 / 937 | 1024x4096; 1024 | reference | 33,438,272 | 1.000 | 25.3218 | 1.000 |
| dense-1024x4096 / 937 | 1024x4096; 1024 | candidate | 29,517,311 | 0.883 | 12.4689 | 0.496 [0.480, 0.514] |
| dense-256x256 / 17 | 256x256; 256 | reference | 130,248 | 1.000 | 0.3871 | 1.000 |
| dense-256x256 / 17 | 256x256; 256 | candidate | 81,484 | 0.626 | 0.3353 | 0.948 [0.858, 1.075] |
| dense-256x256 / 937 | 256x256; 256 | reference | 130,216 | 1.000 | 0.3808 | 1.000 |
| dense-256x256 / 937 | 256x256; 256 | candidate | 81,586 | 0.627 | 0.3407 | 0.852 [0.740, 0.932] |
| dense-4096x512 / 17 | 4096x512; 512 | reference | 8,381,456 | 1.000 | 15.2685 | 1.000 |
| dense-4096x512 / 17 | 4096x512; 512 | candidate | 4,715,999 | 0.563 | 12.0671 | 0.777 [0.750, 0.797] |
| dense-4096x512 / 937 | 4096x512; 512 | reference | 8,388,480 | 1.000 | 15.4025 | 1.000 |
| dense-4096x512 / 937 | 4096x512; 512 | candidate | 4,718,382 | 0.562 | 12.0250 | 0.775 [0.750, 0.793] |
| dense-512x1024 / 17 | 512x1024; 512 | reference | 2,094,656 | 1.000 | 2.7172 | 1.000 |
| dense-512x1024 / 17 | 512x1024; 512 | candidate | 1,636,540 | 0.781 | 1.8978 | 0.690 [0.662, 0.720] |
| dense-512x1024 / 937 | 512x1024; 512 | reference | 2,090,192 | 1.000 | 2.6607 | 1.000 |
| dense-512x1024 / 937 | 512x1024; 512 | candidate | 1,633,199 | 0.781 | 1.7670 | 0.654 [0.640, 0.666] |
| macaulay-n5-l2-m3-d4 / 17 | 395x482; 333 | reference | 52,424 | 1.000 | 0.2230 | 1.000 |
| macaulay-n5-l2-m3-d4 / 17 | 395x482; 333 | candidate | 20,963 | 0.400 | 0.2214 | 1.059 [0.962, 1.242] |
| macaulay-n5-l2-m3-d4 / 937 | 394x482; 333 | reference | 47,192 | 1.000 | 0.2163 | 1.000 |
| macaulay-n5-l2-m3-d4 / 937 | 394x482; 333 | candidate | 19,424 | 0.412 | 0.2109 | 0.980 [0.957, 1.009] |
| macaulay-n9-l6-m2-d4 / 17 | 920x784; 712 | reference | 1,849,627 | 1.000 | 3.0990 | 1.000 |
| macaulay-n9-l6-m2-d4 / 17 | 920x784; 712 | candidate | 683,168 | 0.369 | 2.4016 | 0.710 [0.628, 0.778] |
| macaulay-n9-l6-m2-d4 / 937 | 711x764; 648 | reference | 1,293,348 | 1.000 | 2.2668 | 1.000 |
| macaulay-n9-l6-m2-d4 / 937 | 711x764; 648 | candidate | 468,017 | 0.362 | 1.7277 | 0.768 [0.753, 0.785] |
| sparse-1024x2048 / 17 | 1024x2048; 1024 | reference | 16,208,864 | 1.000 | 14.8227 | 1.000 |
| sparse-1024x2048 / 17 | 1024x2048; 1024 | candidate | 12,281,058 | 0.758 | 8.4147 | 0.567 [0.560, 0.575] |
| sparse-1024x2048 / 937 | 1024x2048; 1024 | reference | 16,316,960 | 1.000 | 15.1631 | 1.000 |
| sparse-1024x2048 / 937 | 1024x2048; 1024 | candidate | 12,387,440 | 0.759 | 8.4572 | 0.556 [0.552, 0.560] |

The larger real F4 holdout (711x764, rank 648) uses
63.8% fewer word XORs and takes
0.768x the baseline kernel time
(1.30x faster). The smaller real Macaulay case
has no resolved timing improvement despite fewer XORs. Its development timing
regresses in the point estimate, and that regression remains in the table.
Dense/wide synthetic cases are not evidence that real F4 matrices have the same
structure or that the attack gains the same speedup.

## Complete native solver stage

Cold times include field setup, encoding and complete F4 solving. Independent
exhaustive truth evaluation and a second diagnostic root reduction are measured
separately, not double-charged to the solver. The solver already checks its own
answers. Values are medians of three runs of the entire target batch.

| n / ell / m | Split | Variant | Cold stage ms | Diagnostic root XORs | Paired time C/R [95% CI] |
| --- | --- | --- | --- | --- | --- |
| 3 / 1 / 2 | frozen | reference | 0.0811 | 98 | 1.000 |
| 3 / 1 / 2 | frozen | candidate | 0.0788 | 98 | 1.007 [0.972, 1.053] |
| 3 / 1 / 2 | fresh | reference | 0.0677 | 104 | 1.000 |
| 3 / 1 / 2 | fresh | candidate | 0.0714 | 104 | 1.041 [1.004, 1.065] |
| 3 / 2 / 2 | frozen | reference | 0.2293 | 443 | 1.000 |
| 3 / 2 / 2 | frozen | candidate | 0.2292 | 443 | 0.985 [0.947, 1.010] |
| 3 / 2 / 2 | fresh | reference | 0.2351 | 488 | 1.000 |
| 3 / 2 / 2 | fresh | candidate | 0.2370 | 488 | 1.003 [0.926, 1.077] |
| 5 / 2 / 2 | frozen | reference | 0.7666 | 3,567 | 1.000 |
| 5 / 2 / 2 | frozen | candidate | 0.7539 | 3,567 | 0.990 [0.982, 0.996] |
| 5 / 2 / 2 | fresh | reference | 0.7658 | 3,612 | 1.000 |
| 5 / 2 / 2 | fresh | candidate | 0.7305 | 3,612 | 0.940 [0.887, 1.007] |
| 7 / 3 / 2 | frozen | reference | 2.1101 | 16,093 | 1.000 |
| 7 / 3 / 2 | frozen | candidate | 2.1534 | 16,093 | 1.022 [1.012, 1.028] |
| 7 / 3 / 2 | fresh | reference | 2.0748 | 16,496 | 1.000 |
| 7 / 3 / 2 | fresh | candidate | 2.0495 | 16,496 | 0.993 [0.983, 1.006] |
| 9 / 3 / 2 | frozen | reference | 1.2108 | 21,445 | 1.000 |
| 9 / 3 / 2 | frozen | candidate | 1.2196 | 21,445 | 1.002 [0.952, 1.047] |
| 9 / 3 / 2 | fresh | reference | 1.2535 | 22,168 | 1.000 |
| 9 / 3 / 2 | fresh | candidate | 1.2344 | 22,168 | 0.987 [0.962, 1.007] |
| 3 / 1 / 3 | frozen | reference | 0.3566 | 623 | 1.000 |
| 3 / 1 / 3 | frozen | candidate | 0.3603 | 623 | 1.013 [1.009, 1.016] |
| 3 / 1 / 3 | fresh | reference | 0.8036 | 559 | 1.000 |
| 3 / 1 / 3 | fresh | candidate | 0.7199 | 559 | 0.965 [0.888, 1.016] |
| 5 / 2 / 3 | frozen | reference | 65.2960 | 46,919 | 1.000 |
| 5 / 2 / 3 | frozen | candidate | 67.4180 | 32,599 | 1.029 [1.018, 1.040] |
| 5 / 2 / 3 | fresh | reference | 66.1600 | 47,140 | 1.000 |
| 5 / 2 / 3 | fresh | candidate | 67.2988 | 32,744 | 1.023 [1.002, 1.050] |

Most frozen systems have one-word rows and cannot benefit from skipping earlier
words. For n=5, ell=2, m=3 the diagnostic root XOR count falls about 30.5%, but
complete stage time is about 2–3% higher in this run. Kernel savings therefore do
not establish a complete-solver gain on this corpus.

## Full DLP attempts

Every process starts cold, including setup and final scalar recovery. Raw JSON
retains the phase timings and process CPU/RSS. Phase timers are subsets, not an
exhaustive calibrated budget. The paired intervals below include one for every
full-DLP group, so they do not support a runtime improvement claim either.

| n / solver | Split | Variant | Total elapsed ms | Verified / executions | Paired time C/R [95% CI] |
| --- | --- | --- | --- | --- | --- |
| 7 / groebner | frozen | reference | 435.725 | 0 / 12 | 1.000 |
| 7 / groebner | frozen | candidate | 437.076 | 0 / 12 | 1.003 [0.983, 1.025] |
| 7 / groebner | fresh | reference | 214.975 | 0 / 6 | 1.000 |
| 7 / groebner | fresh | candidate | 214.638 | 0 / 6 | 0.998 [0.994, 1.003] |
| 7 / sat | frozen | reference | 336.883 | 0 / 12 | 1.000 |
| 7 / sat | frozen | candidate | 334.680 | 0 / 12 | 0.994 [0.970, 1.013] |
| 7 / sat | fresh | reference | 168.000 | 0 / 6 | 1.000 |
| 7 / sat | fresh | candidate | 166.591 | 0 / 6 | 0.991 [0.958, 1.036] |
| 9 / groebner | frozen | reference | 66.308 | 12 / 12 | 1.000 |
| 9 / groebner | frozen | candidate | 66.494 | 12 / 12 | 1.004 [0.991, 1.014] |
| 9 / groebner | fresh | reference | 36.528 | 6 / 6 | 1.000 |
| 9 / groebner | fresh | candidate | 37.730 | 6 / 6 | 1.024 [0.969, 1.091] |
| 9 / sat | frozen | reference | 29.117 | 12 / 12 | 1.000 |
| 9 / sat | frozen | candidate | 28.603 | 12 / 12 | 0.987 [0.946, 1.026] |
| 9 / sat | fresh | reference | 12.467 | 6 / 6 | 1.000 |
| 9 / sat | fresh | candidate | 12.525 | 6 / 6 | 1.004 [0.980, 1.030] |

The n=7 factor base has no usable columns; those attempts remain incomplete in
both versions. The SAT path is a control and does not call the changed kernel.
All n=9 planted scalars verify `[k]P = Q`, including fresh seeds 503 and 607 with
scalars 7 and 11. These are toy groups and do not establish scaling.

The measured improvement is confined to elimination. Profiling larger real
Macaulay workloads is the next useful step before choosing a more complex
blocked elimination method or optimizing the final sparse relation solve.

## Evidence

- [Frozen contract](contract.json), [runner](run.py), [comparison](compare.py).
- [Accepted comparison](results/run-002/summary.json),
  [input/source/binary hashes](results/run-002/hashes.json),
  [hardware](results/run-002/hardware.json),
  [all process outcomes](results/run-002/processes.json).
- [Rejected first run](results/run-001/REJECTED.md),
  [tests](validation/tests.txt), [dependency lock](validation/dependencies.lock.txt).
- [Canonical scoreboard](../../docs/index-calculus-scoreboard.html#f4-linear-algebra-20260914).
