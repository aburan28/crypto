#set page(paper: "a4", margin: (x: 1.8cm, y: 1.8cm), numbering: "1")
#set text(size: 9.5pt)
#set heading(numbering: "1.")
#show table: set text(size: 8pt)

#align(center)[
  #text(size: 14pt, weight: "bold")[Run report: endo-chain-bench counting and timing runs, 2026-10-16]
  #v(2pt)
  `research/endo_native_bench_20261016/` · first run started 2026-10-16 17:59 UTC · commit `33e2dd10bd8b75a13122683576bf5a81f1583577`
]

Descriptive record of the six runs (one counting and one timing run per
curve). Observations only; the interpretation is in the study's results
report (`../../REPORT.pdf`) and `README.md`. Every value is copied from the
committed records; the graphs are drawn from them by `../../figures/plot.py`.

= What ran

- *Command.* `endo-chain-bench` (harness `85793758e`), protocol committed in
  `PROTOCOL.md` at `33e2dd10b` before any run on these curves. Counting:
  `target/opcount/release/endo-chain-bench --count-only --chains … --expect-sha256 … --model … --reference … --pairs 256 --seed 20261016`.
  Timing: the same binary without the opcount feature, `--rounds 20`, under
  `isolated_bench run --wait --cpus 3 --settle 5 --max-other-cpu 0.5`. The
  full command list is in `README.md`, "Commands".
- *Build.* `rustc 1.97.0 (2d8144b78 2026-07-07)`, release, no `target-cpu`;
  field over `Uint<L>::mont_mul` with 4 / 6 / 13 limbs.
- *Host.* Intel(R) Xeon(R) Processor @ 2.10GHz, cloud VM, 4 logical CPUs,
  kernel 6.18.44-fc-v139. Isolation earned: L0 (shared VM, timing pinned to
  CPU 3 and admitted by `isolated_bench`; counting runs not isolated, counts
  do not depend on contention).
- *Workloads.* 256 $(k, P)$ pairs per curve, seed 20261016; timing 20 rounds,
  arm order forward on even and reversed on odd rounds; statistic median and
  minimum of round wall time / pairs.
- *Arms.* Baselines (width-$w$ NAF, Jacobian and affine tables, $w = 3..7$),
  GLV-2 grid on the declared chain (generic / optimised evaluator × table ×
  width), GLV elements at `optimised/jacobian/w5`, $phi(P)$ and
  decomposition stage diagnostics, one A/A arm. 38 / 96 / 180 arms.
- *Controls.* A/A arm (the reference re-timed); positive control before the
  protocol: the generic harness's counts equal the CryptoPro-B sweep's
  committed counts on all 97 arms (`positive_control/`).

= Runs requested and completed

#table(columns: (auto, auto, auto, auto, auto, auto, auto, auto), stroke: 0.4pt + gray,
  [*curve*], [*arms*], [*all correct (count / time)*], [*timing wall s (isolated)*], [*`contended`*], [*other CPU-s*], [*invol. switches*], [*major faults*],
  [Bandersnatch], [38], [true / true], [23.3], [false], [0.48], [55], [0],
  [Tom-384], [96], [true / true], [98.9], [false], [1.79], [136], [0],
  [CP6-782 G1], [180], [true / true], [570.2], [false], [11.11], [621], [0],
)

6 of 6 runs completed, each admitted at the first attempt. One earlier
attempt to launch the counting runs under `/usr/bin/time` failed before the
bench started because that tool is not installed; it produced no record.

= Measured values

#figure(image("../../figures/fig_counted_vs_model.svg", width: 100%),
  caption: [Counting runs. Unit: field operations $M_"eq" = M + S$ per $k dot P$, mean over 256 pairs. (a) Derived: min counted $M_"eq"$ over baseline arms / min over GLV arms, and the same from the model column. (b) Measured deviation from the frozen model for every total arm.])

#figure(image("../../figures/fig_timed.svg", width: 100%),
  caption: [Timing runs, n = 20 rounds per arm. (a) Recorded `ratio_reference_over_arm` of the declared GLV arm (counted, median, min); hatched band: counted ± $|1 - "A/A median ratio"|$. (b) Per-round ns/op ÷ reference median.])

#figure(image("../../figures/fig_grid.svg", width: 100%),
  caption: [Grid of the declared chain: counted mean $M_"eq"$ and median µs per $k dot P$ against width.])

#figure(image("../../figures/fig_phi.svg", width: 100%),
  caption: [$phi(P)$ stage diagnostics: counted $M_"eq"$ per call (identical on every input) against median µs per call.])

Recorded values for the declared GLV arm (README "Headline"): counted
$M_"eq"$ 2 158.6 / 3 240.4 / 2 991.5 against reference 3 289.8 / 4 929.6 /
4 462.4; median ns 94 059 / 310 273 / 1 366 596 against 141 443 / 466 490 /
2 175 424; A/A median ratios 0.997 / 1.001 / 1.038. Per-round spread
(max/min − 1) on the reference 38 %, 32 %, 49 %; on the declared arm 85 %,
55 %, 66 %.

= Raw records

`../../bandersnatch.counts.json`, `../../tom384.counts.json`,
`../../cp6_782.counts.json`; `../../<curve>.timing.json`;
`../../<curve>.isolated_bench_record.jsonl`; `.stdout.md` tables beside
each; inputs in `../../frozen_inputs/` (`SHA256SUMS`).
