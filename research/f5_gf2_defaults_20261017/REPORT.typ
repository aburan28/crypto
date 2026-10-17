#set page(paper: "a4", margin: 1.5cm)
#set text(size: 8.5pt)
= Matrix F5 / GF(2) defaults: direct pack on; direct unpack and table reuse stay opt-in
Study `research/f5_gf2_defaults_20261017` · runs `run-20261017T0331Z` (3 rounds) and `run-20261017T0341Z` (15 rounds; F4 10 rounds) · binaries `f4_f2_bench`, `f6_ic_probe` built from `78a90d27f` · 2026-10-17

== Question
Three exact switches had each passed an incremental gate on x86 CI and were left opt-in only because their protocols tied a default change to a 2× gate: `KIC_F5_DIRECT_PACK` (`research/f5_direct_pack_20260929`, 1.21×), `KIC_F5_UNPACK_DIRECT` (`research/f5_further_2x_20260929`, 1.04–1.05×), `KIC_GF2_REUSE_TABLE` (`research/gf2_table_reuse_20260929`, 1.035×). Can they be on by default without regressing on another host and on the other `gf2_elim` users?

== What ran
One binary per tool (all three switches on at `78a90d27f`), arms by environment only: *on*; *nopack*, *nounpack*, *noreuse* (that switch `=0`); *off* (all three `=0`, the previous defaults). `f4_f2_bench 3 24 f5` (seven matrix-F5 cases, n ≤ 24, `Reduced` = production form and `SelectiveEchelon`, 1 and 4 threads), `f4_f2_bench f4` (matrix F4) and, in run 1, `f6_ic_probe --engines f4` (inherited F4) on three cells. Arm order rotated each round. Host: Intel Xeon @ 2.10 GHz, 4 vCPU KVM, AVX-512F, shared, unpinned, L0, load 0.36–3.23. Run 1: 84 runs; run 2: 260 runs; all exited 0. Times are the bench's per-case median of 3 repetitions (wall) and the probe's process CPU; ratios are medians over rounds (arm / on, > 1 means the defaults-on arm is faster), the range is the worst–best round pairing.

== Exactness
Every F5 case (row and row-space fingerprints, rank, output terms, reduction word operations), every matrix-F4 case (basis fingerprint, word XORs) and every inherited-F4 target decision is identical in every arm and round of both runs.

== Decision rule
Run 1 left direct unpack and table reuse ambiguous (on `f5_n20_m20_d4` single switches off looked faster than `on` while `off` was slower). Before run 2 the rule in `scripts/DECISION_RULE.md` was committed: a switch stays on only if in every (form, threads) cell its arm/on median is ≥ 0.97 on `f5_n20_m20_d4` and `f5_n24_m24_d4` and ≥ 1.00 on `f5_n24_m24_d4`, and, for table reuse, ≥ 0.97 on every matrix-F4 case with n ≥ 16.

== Outcome (run 2)
- *Direct unpack*: fails: arm/on 0.986 on `f5_n24_m24_d4`, SelectiveEchelon, 4 threads (other cells 0.974–1.060). Stays opt-in.
- *Table reuse*: fails: arm/on 0.955 on matrix-F4 `n16_m24` (F5 cells 0.995–1.070 pass). Stays opt-in.
- *Direct pack*: on by default (run 1, 3 rounds; not repeated in run 2: arm/on 1.006–2.173 on every F5 cell; on `f5_n24_m24_d4` the worst round pairing is ≥ 1.02 in all four cells).
- Previous defaults vs all three on (off ÷ on), run 2: 1.11–1.27× on the two largest F5 cases in every cell.

== Run 2, larger cases
#table(columns: 9, align: (left, left, left, right, left, right, right, right, right),
  [arm], [tool], [form], [thr], [case], [on ms], [arm ms], [arm/on], [round range],
  [noreuse], [f5], [form0], [1], [`f5_n20_m20_d4`], [72.50], [73.94], [1.020], [0.75–1.44],
  [nounpack], [f5], [form0], [1], [`f5_n20_m20_d4`], [72.50], [72.08], [0.994], [0.77–1.57],
  [off], [f5], [form0], [1], [`f5_n20_m20_d4`], [72.50], [82.71], [1.141], [0.82–1.47],
  [noreuse], [f5], [form0], [1], [`f5_n24_m24_d4`], [273.83], [285.74], [1.044], [0.82–1.21],
  [nounpack], [f5], [form0], [1], [`f5_n24_m24_d4`], [273.83], [277.97], [1.015], [0.84–1.28],
  [off], [f5], [form0], [1], [`f5_n24_m24_d4`], [273.83], [315.05], [1.151], [0.96–1.37],
  [noreuse], [f5], [form0], [4], [`f5_n20_m20_d4`], [47.76], [51.11], [1.070], [0.73–1.26],
  [nounpack], [f5], [form0], [4], [`f5_n20_m20_d4`], [47.76], [50.11], [1.049], [0.73–1.27],
  [off], [f5], [form0], [4], [`f5_n20_m20_d4`], [47.76], [55.12], [1.154], [0.83–1.43],
  [noreuse], [f5], [form0], [4], [`f5_n24_m24_d4`], [165.71], [171.92], [1.038], [0.83–1.16],
  [nounpack], [f5], [form0], [4], [`f5_n24_m24_d4`], [165.71], [167.08], [1.008], [0.84–1.17],
  [off], [f5], [form0], [4], [`f5_n24_m24_d4`], [165.71], [183.82], [1.109], [0.92–1.25],
  [noreuse], [f5], [form2], [1], [`f5_n20_m20_d4`], [48.79], [48.57], [0.995], [0.77–1.41],
  [nounpack], [f5], [form2], [1], [`f5_n20_m20_d4`], [48.79], [47.54], [0.974], [0.79–1.36],
  [off], [f5], [form2], [1], [`f5_n20_m20_d4`], [48.79], [55.29], [1.133], [0.79–1.41],
  [noreuse], [f5], [form2], [1], [`f5_n24_m24_d4`], [145.42], [149.97], [1.031], [0.84–1.31],
  [nounpack], [f5], [form2], [1], [`f5_n24_m24_d4`], [145.42], [154.11], [1.060], [0.85–1.23],
  [off], [f5], [form2], [1], [`f5_n24_m24_d4`], [145.42], [183.47], [1.262], [1.05–1.59],
  [noreuse], [f5], [form2], [4], [`f5_n20_m20_d4`], [33.32], [34.62], [1.039], [0.86–1.25],
  [nounpack], [f5], [form2], [4], [`f5_n20_m20_d4`], [33.32], [33.03], [0.991], [0.80–1.26],
  [off], [f5], [form2], [4], [`f5_n20_m20_d4`], [33.32], [42.41], [1.273], [1.04–1.60],
  [noreuse], [f5], [form2], [4], [`f5_n24_m24_d4`], [101.08], [105.84], [1.047], [0.93–1.27],
  [nounpack], [f5], [form2], [4], [`f5_n24_m24_d4`], [101.08], [99.66], [0.986], [0.87–1.18],
  [off], [f5], [form2], [4], [`f5_n24_m24_d4`], [101.08], [122.68], [1.214], [1.03–1.51],
  [noreuse], [f4], [n/a], [1], [`n12_m12`], [7.35], [7.56], [1.028], [0.76–1.53],
  [noreuse], [f4], [n/a], [1], [`n14_m14`], [24.61], [25.93], [1.054], [0.85–1.35],
  [noreuse], [f4], [n/a], [1], [`n16_m16`], [298.01], [293.50], [0.985], [0.83–1.32],
  [noreuse], [f4], [n/a], [1], [`n16_m24`], [77.73], [74.22], [0.955], [0.75–1.08],
  [noreuse], [f4], [n/a], [1], [`n18_m18`], [1776.17], [1772.28], [0.998], [0.92–1.08],
)

== Run 1, larger cases and the other engines
#table(columns: 9, align: (left, left, left, right, left, right, right, right, right),
  [arm], [tool], [form], [thr], [case], [on ms], [arm ms], [arm/on], [round range],
  [nopack], [f5], [form0], [1], [`f5_n20_m20_d4`], [71.80], [82.29], [1.146], [1.00–1.52],
  [noreuse], [f5], [form0], [1], [`f5_n20_m20_d4`], [71.80], [70.80], [0.986], [0.86–1.05],
  [nounpack], [f5], [form0], [1], [`f5_n20_m20_d4`], [71.80], [72.30], [1.007], [0.86–1.09],
  [off], [f5], [form0], [1], [`f5_n20_m20_d4`], [71.80], [85.87], [1.196], [0.95–1.36],
  [nopack], [f5], [form0], [1], [`f5_n24_m24_d4`], [266.60], [324.23], [1.216], [1.08–1.36],
  [noreuse], [f5], [form0], [1], [`f5_n24_m24_d4`], [266.60], [285.15], [1.070], [0.91–1.19],
  [nounpack], [f5], [form0], [1], [`f5_n24_m24_d4`], [266.60], [272.34], [1.022], [0.87–1.09],
  [off], [f5], [form0], [1], [`f5_n24_m24_d4`], [266.60], [314.27], [1.179], [1.02–1.24],
  [nopack], [f5], [form0], [4], [`f5_n20_m20_d4`], [52.85], [56.72], [1.073], [0.98–1.22],
  [noreuse], [f5], [form0], [4], [`f5_n20_m20_d4`], [52.85], [46.68], [0.883], [0.80–0.99],
  [nounpack], [f5], [form0], [4], [`f5_n20_m20_d4`], [52.85], [45.74], [0.865], [0.83–1.02],
  [off], [f5], [form0], [4], [`f5_n20_m20_d4`], [52.85], [57.70], [1.092], [1.05–1.35],
  [nopack], [f5], [form0], [4], [`f5_n24_m24_d4`], [159.30], [175.58], [1.102], [1.02–1.11],
  [noreuse], [f5], [form0], [4], [`f5_n24_m24_d4`], [159.30], [157.93], [0.991], [0.90–1.05],
  [nounpack], [f5], [form0], [4], [`f5_n24_m24_d4`], [159.30], [161.41], [1.013], [0.93–1.10],
  [off], [f5], [form0], [4], [`f5_n24_m24_d4`], [159.30], [173.86], [1.091], [1.02–1.28],
  [nopack], [f5], [form2], [1], [`f5_n20_m20_d4`], [52.97], [53.32], [1.006], [0.99–1.09],
  [noreuse], [f5], [form2], [1], [`f5_n20_m20_d4`], [52.97], [48.05], [0.907], [0.88–0.96],
  [nounpack], [f5], [form2], [1], [`f5_n20_m20_d4`], [52.97], [45.98], [0.868], [0.80–0.92],
  [off], [f5], [form2], [1], [`f5_n20_m20_d4`], [52.97], [56.60], [1.068], [1.03–1.29],
  [nopack], [f5], [form2], [1], [`f5_n24_m24_d4`], [145.97], [164.72], [1.129], [1.02–1.20],
  [noreuse], [f5], [form2], [1], [`f5_n24_m24_d4`], [145.97], [153.58], [1.052], [0.96–1.18],
  [nounpack], [f5], [form2], [1], [`f5_n24_m24_d4`], [145.97], [149.02], [1.021], [0.90–1.03],
  [off], [f5], [form2], [1], [`f5_n24_m24_d4`], [145.97], [196.95], [1.349], [1.11–1.48],
  [nopack], [f5], [form2], [4], [`f5_n20_m20_d4`], [32.01], [42.74], [1.335], [1.29–1.57],
  [noreuse], [f5], [form2], [4], [`f5_n20_m20_d4`], [32.01], [33.16], [1.036], [0.96–1.14],
  [nounpack], [f5], [form2], [4], [`f5_n20_m20_d4`], [32.01], [33.77], [1.055], [1.04–1.12],
  [off], [f5], [form2], [4], [`f5_n20_m20_d4`], [32.01], [42.91], [1.341], [1.06–1.59],
  [nopack], [f5], [form2], [4], [`f5_n24_m24_d4`], [98.92], [112.48], [1.137], [1.03–1.26],
  [noreuse], [f5], [form2], [4], [`f5_n24_m24_d4`], [98.92], [105.06], [1.062], [0.95–1.14],
  [nounpack], [f5], [form2], [4], [`f5_n24_m24_d4`], [98.92], [100.44], [1.015], [0.98–1.09],
  [off], [f5], [form2], [4], [`f5_n24_m24_d4`], [98.92], [122.11], [1.234], [1.18–1.32],
  [noreuse], [f4], [n/a], [1], [`n12_m12`], [7.84], [7.00], [0.893], [0.76–1.00],
  [noreuse], [f4], [n/a], [1], [`n14_m14`], [25.41], [24.78], [0.975], [0.92–1.00],
  [noreuse], [f4], [n/a], [1], [`n16_m16`], [293.54], [308.31], [1.050], [0.95–1.18],
  [noreuse], [f4], [n/a], [1], [`n16_m24`], [81.93], [76.85], [0.938], [0.73–1.11],
  [noreuse], [f4], [n/a], [1], [`n18_m18`], [1793.38], [1919.97], [1.071], [1.03–1.09],
  [noreuse], [f4], [n/a], [1], [`n20_m20`], [19337.74], [19604.29], [1.014], [0.99–1.04],
  [noreuse], [f4], [n/a], [1], [`n20_m30`], [1478.42], [1371.63], [0.928], [0.87–1.05],
  [noreuse], [inh], [n/a], [1], [`K0_23`], [2403.36], [2409.95], [1.003], [0.87–1.10],
  [noreuse], [inh], [n/a], [1], [`K0_31`], [2483.26], [2508.61], [1.010], [0.92–1.08],
  [noreuse], [inh], [n/a], [1], [`K1_23`], [810.36], [825.74], [1.019], [1.00–1.04],
)

#figure(image("fig_f5.svg", width: 100%), caption: [Median arm/on ratios on the two largest F5 cases, both runs.])

== Limits
One shared 4-vCPU host, L0, wall times; run-to-run ranges are wide (single runs can move ±25 %). Direct pack's measured gain here is largest on the smaller F5 cases. No IC-versus-rho statement.

== Records
`reports/<run>/raw/*.jsonl` (first line: tool, arm, form, threads, round, UTC, load), `raw/*.status`, `host.txt`, `summary.txt`, `table.tsv` (`summarize.py`), `fig_f5.svg|png` (`plot.py`); scripts in `scripts/`.
