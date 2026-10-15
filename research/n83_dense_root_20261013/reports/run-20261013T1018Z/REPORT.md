# N83 dense-root 64-point panel (2026-10-13 10:18–10:50 UTC)

K1798 dense met the frozen lowest-mean-probe rule among completed arms:
169,069.8125 probes per target across the same 64 public points, versus
300,021.90625 for K1500 dense and 371,728.25 for both K1271 layouts. The
paired K1500-minus-K1798 difference is 130,952.09375 probes, with a
10,000-resample 95% percentile interval [52,057.703125, 209,990.828125].
All 256 rows in the four completed arms passed independent generic-curve
scalar and four-point relation replay. The repeated K1500 regular control
reached exit 124 at the 900-second cap before a target row.

![Panel probes, capacity and exploratory online time](panel.svg)

![Paired probe differences and bootstrap intervals](panel_paired.svg)

| Arm | Actual B | Columns | Table allocation | Mean probes | Mean online ms (exploratory) | Full process s | Status |
|---|---:|---:|---:|---:|---:|---:|---|
| K1271 regular | 210,986 | 1,271 | 8 GiB | 371,728.25 | 171.706 | 177.480 | 64/64 independent replay pass |
| K1271 dense | 210,986 | 1,271 | 4 GiB | 371,728.25 | 151.812 | 183.176 | 64/64 independent replay pass |
| K1500 dense | 249,000 | 1,500 | 8 GiB | 300,021.90625 | 129.283 | 242.784 | 64/64 independent replay pass |
| K1500 regular | missing at cap | missing at cap | no completed receipt | missing | missing | 901.21 wall from `time -l` | timeout 124; max RSS 12,899,270,656 B |
| K1798 dense | 298,468 | 1,798 | 8 GiB | 169,069.8125 | 76.774 | 391.965 | 64/64 independent replay pass |

The 64 public point coordinates were frozen before any arm using SHA-256
`e9ed7ddb0aeafdbfc7354d67279a8712572eea9d2aee3896f99633bdc2e937cc`.
Their exact BLAKE3-indexed scalar law and separate answer sidecar are in
`../../PROTOCOL.md`. The producer never receives the sidecar. All arms use
the same degree-83, `a=b=1` Koblitz curve, rank seed 29, 12 rank threads,
one target thread, source binary SHA-256
`caa723d4e2894b8573c04f6ec2ec474ce4fb41eb8a2c452b793b3d2d286aa351`,
and 900-second process cap. The four-arm run order was K1500 regular, K1500
dense, K1271 dense, K1271 regular; the K1798 dense capacity extension was
committed before this panel started and ran after the four arms. Source
freezes are commits `4489ffdba3`, `48f3892294`, and `0f8f99ab63`.

The native analyzer (`../../analyzer/`) replayed all 64 recovered scalars
for every completed arm using the generic curve library. It checked every
archived base representative in the prime subgroup, all 166 signed
Frobenius-orbit points and labels per representative, point distinctness,
and each target's four-point sum against its public input. The archived
bases are the K1271 regular, K1500 dense and K1798 dense runs; matching
base hashes bind the alternate K1271 dense arm to the replayed base. K1271
regular and dense also returned identical first-witness indices, x codes,
probe counts and scalars on 64/64 paired targets. K1500 table parity is
censored by the regular-arm timeout. The largest target phase-sum residual
among completed arms was 0.000501 ms.

The paired-probe differences (left minus right; positive favors right) are:

| Pair | Mean probes | 95% paired bootstrap interval |
|---|---:|---:|
| K1271 dense − K1500 dense | 71,706.34375 | [−24,081.171875, 168,322.5] |
| K1500 dense − K1798 dense | 130,952.09375 | [52,057.703125, 209,990.828125] |
| K1271 dense − K1798 dense | 202,658.4375 | [121,448.546875, 288,107.375] |
| K1271 regular − K1271 dense | 0 | [0, 0] |

The native analyzer used 10,000 paired bootstrap resamples, SplitMix64 seed
`20261013256`, index `random_u64 % 64`, and sorted positions 250 and 9750.
These intervals describe this deterministic target law and run; they do not
turn shared-host wall times into controlled speedups. CPU isolation was not
established on this Apple M4 Pro host, and another research process was
active. The one-target IC/rho online comparison remains in the primary
[run report](../run-20261013T0942Z/REPORT.md); this panel is a secondary
factor-base/probe study.

Raw evidence: `../../runs/panel/` and `../../analysis_panel.json`.
Archived base SHA-256 values are K1271
`affbdaa1b0437daee986f7b13dd16347052f67dcbebd3772248b33fa7e326128`,
K1500 `125c5d2198994e80c7b685db942c0da81d078ec0a3fb4199ae908e536aee62cd`,
and K1798 `a9707636d4bef65e80b3c772ce29b29f2259c4b8e5e2eedfade8b14c44ab22ab`.
The arithmetic capacity theorem and its formalization status are in
`../../PROOF_MAP.md`.
