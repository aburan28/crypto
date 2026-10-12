#set page(paper: "a4", margin: (x: 1.8cm, y: 2cm), numbering: "1")
#set text(size: 9.5pt)
#set table(stroke: 0.4pt, inset: 4pt)
#show heading: set block(above: 1.2em, below: 0.6em)

#align(center)[
  #text(size: 15pt, weight: "bold")[Run report: run-20261011T1955Z-verify] \
  Eigen-Séta v2: correctness rates, constant-time probes at Toy16 and L1, timing test, filter rates (descriptive record) \
  #text(size: 8.5pt)[Study: research/eigen_seta_rust_20261011 · sources: commit d504af28e (filter step: 72051ac67)]
]

= What ran

#table(columns: (auto, 1fr),
  [Times (UTC)], [Start 19:56:36Z, end of phase 1 20:29:49Z, end 20:37:48Z (2026-10-11). Filter step 20:38:30Z–20:46:52Z.],
  [Sources], [Commit d504af28e with a clean worktree. The filter step used 72051ac67, which differs only by capping the public basis search and adding the filter example.],
  [Build / tools], [rustc 1.97.0, release with debuginfo. valgrind-3.22.0 (callgrind), glibc malloc.],
  [Host], [Intel Xeon \@ 2.80 GHz, 4 logical cores, Linux 6.18.44 x86_64, shared cloud VM],
  [Phase 1], [4 callgrind probes at L1 in parallel with native failure-rate runs: Toy16 ×3000, Toy32 ×600, L1 ×20],
  [Phase 2], [12-run Toy16 probe grid: keys 1–4 × (valid 1, valid 2, foreign)],
  [Phase 3], [`dudect_decaps 1000`, alone on the machine],
  [Probe], [`valgrind --tool=callgrind --toggle-collect='*probe_core*' ./ct_probe <key> <ct> <mode> <set>`. Only `decaps_core` is measured. Allocator = glibc malloc and arena code, sbrk, brk, morecore, mem\*.],
  [Records], [`records/` (failrate/, toyprobe/, l1probe/ with callgrind/\*.xz, dudect_output.txt, f_filter_rate.txt, verify_run.sh, ctsummary.sh, build.log, commit.txt)],
)
Requested and completed: 4 + 12 probes, 3,620 trials, 1 timing test, 5 filter measurements. Failures: none.

= Decapsulation correctness

#table(columns: 4,
  [*Set*], [*Trials*], [*Correct*], [*Failures*],
  [Toy16], [3,000], [3,000], [0],
  [Toy32], [600], [600], [0],
  [L1], [20], [20], [0],
)

= Instructions in `decaps_core`

#table(columns: 7,
  [*set*], [*key*], [*ct*], [*acc.*], [*total*], [*allocator*], [*rest*],
  [Toy16], [1], [valid 1], [yes], [1,147,950,979], [31,741,763], [1,116,209,216],
  [Toy16], [1], [valid 2], [yes], [1,147,950,979], [31,741,763], [1,116,209,216],
  [Toy16], [1], [foreign], [no], [1,148,059,231], [31,850,015], [1,116,209,216],
  [Toy16], [2], [valid 1], [yes], [1,147,388,400], [31,179,184], [1,116,209,216],
  [Toy16], [2], [valid 2], [yes], [1,147,388,400], [31,179,184], [1,116,209,216],
  [Toy16], [2], [foreign], [no], [1,147,693,051], [31,483,835], [1,116,209,216],
  [Toy16], [3], [valid 1], [yes], [1,147,459,087], [31,249,871], [1,116,209,216],
  [Toy16], [3], [valid 2], [yes], [1,147,459,087], [31,249,871], [1,116,209,216],
  [Toy16], [3], [foreign], [no], [1,147,849,768], [31,640,552], [1,116,209,216],
  [Toy16], [4], [valid 1], [yes], [1,147,477,368], [31,268,152], [1,116,209,216],
  [Toy16], [4], [valid 2], [yes], [1,147,477,368], [31,268,152], [1,116,209,216],
  [Toy16], [4], [foreign], [no], [1,147,787,126], [31,577,910], [1,116,209,216],
  [L1], [1], [valid 1], [yes], [124,250,795,358], [3,254,365,748], [120,996,429,610],
  [L1], [1], [valid 2], [yes], [124,250,795,358], [3,254,365,748], [120,996,429,610],
  [L1], [1], [foreign], [no], [124,257,288,150], [3,260,858,540], [120,996,429,610],
  [L1], [2], [valid 1], [yes], [124,262,712,389], [3,266,282,779], [120,996,429,610],
)
Observation: at each set the rest takes a single value. The allocator count varies with key and mode, and is equal for the two valid ciphertexts of a key. A first L1 summary that did not count sbrk and brk as allocator code showed a 174-instruction difference between the keys, all of it in those libc heap-growth routines.

#figure(image("ct_probe.svg", width: 92%), caption: [Instruction counts per run, Toy16 and L1. Measured; records/{toyprobe,l1probe}/ctprobe_summary.txt.])

= Timing test (ns per call, 1,000 samples per class)

```
samples/class 1000; median ns: fixed 122373850, random 120844155
t(raw) = +1.17
t(crop <= p50) = +0.14  (n = 473, 528)
t(crop <= p75) = +1.46  (n = 730, 771)
t(crop <= p90) = +2.36  (n = 899, 902)
```

= Filter acceptance rates (fraction of odd f candidates)

#table(columns: 5,
  [*Set*], [*Drawn*], [*Accepted*], [*Rate*], [*log2*],
  [Toy16], [12,408], [300], [0.02418], [−5.37],
  [Toy32], [24,956], [300], [0.01202], [−6.38],
  [L1], [55,933], [200], [0.00358], [−8.13],
  [L3], [39,987], [100], [0.00250], [−8.64],
  [L5], [51,887], [100], [0.00193], [−9.02],
)
Rejections: gcd(f,21) ≠ 1 (about 43%) and the two-squares step (52–56% of all drawn). The schedule condition rejects about half of the candidates that reach it.

#figure(image("f_filter.svg", width: 70%), caption: [log2 of the acceptance rate; error bars are approx. 1/(ln 2 · √k). Measured; records/f_filter_rate.txt.])

= Scope

Toy16 and L1 probes only, not L3 or L5. One compiler and one target. Instruction counts, not memory addresses. One key in the timing test. Zero failures in 3,000 trials bounds the Toy16 failure rate only by about 10#super[−3] at 95% confidence.
