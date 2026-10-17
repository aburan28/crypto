#set page(paper: "a4", margin: (x: 2cm, y: 2cm), numbering: "1")
#set text(size: 9.5pt)
#set heading(numbering: "1.")
#show table: set text(size: 8pt)

#align(center)[
  #text(size: 15pt, weight: "bold")[The F4 tower engine over a generic prime field]
  #v(2pt)
  Lane F, lead 4 · 2026-10-15 · branch `claude/f4-generic-field` · `research/f4_generic_field_20261015/`
]

*Summary.* `f4_tower_generic<F: PrimeField>` runs the PKM tower F4 engine of
`f4_fp_tower` over any prime-field backend. At `Fp32` it reproduces
`f4_fp_tower` field for field (182 rows, 0 differences) and is not slower
(paired median ratio 1.09 and 1.06); above $2^31$ it is 2.18–2.24× faster; at
45 bits the new AVX-512 IFMA backend `Fp52` is 2.57× `Fp64`. All five
pre-registered speed claims pass. With it, the Kummer tower runs at
$v_2(p-1) = 40, 48, 56$ and at the NIST P-224 prime ($v_2 = 96$), and the
Dickson torus tower at P-192 and P-256. In all 74 valued cells the solving
degree $D$ — and the width and step count — equal the recorded 31–32-bit
values (H_same supported, none differing); on the fixed binary even the
multiply-add count is identical across the five Kummer primes. Every finished
row is verified (exhaustive decomposition count; points re-checked in
`BigUint`). Host: 4 vCPU Xeon 2.10 GHz, shared with another lane (L0).

= Engine and backends

- `f4_tower_generic`: the `u32` engine's pairs, criteria, symbolic
  preprocessing, Faugère–Lachartre elimination and full-rank exit, with
  coefficients `F::Elem` and dense kernels on Stage 0's lazy accumulators
  (`row_load`/`row_sub`/`row_reduce`, sparse `row_scatter_axpy`); per-row
  product counting against `lazy_capacity`; block sizes independent of the
  backend, so the report is too.
- `Fp52` ($p < 2^52$): `Fp64`'s Montgomery elements; accumulator
  $"lo" + "hi" dot 2^52$ filled by `vpmadd52luq/huq`, 4 entries per 2
  instructions; capacity 4095 products.
- `pkm_tower_pilot --engine generic --backend auto|…` (`Fp32` < $2^32$,
  `Fp52` < $2^52$, `Fp64` < $2^63$, `FpP192`, `FpWide<3>`, `Fp224`,
  `FpWide<4>`), `--p P-192|P-224|P-256|proth:K`; above $2^63$ a builder
  generic over the backend (`pkm_tower_wide`) with `--verify`.

= Correctness

#table(columns: (auto, 1fr), stroke: 0.4pt + gray,
  [*gate*], [*result*],
  [unit tests], [7/7 + 7/7: every backend (`Fp32` vector/scalar, `Fp64`, `FpWide<1>`, `FpWide<4>`, `Fp224`, `Fp52` IFMA/scalar, capacity 1–3 wrappers) = `f4_fp_tower` report and basis on random rings at 5 primes; brute-force refutation and point counts at a 40-bit prime and P-224; backend agreement at 40, 48, 56 bits and P-224; `Fp52` vs `BigUint` through full capacity],
  [G2 equality], [`u32` vs generic `Fp32`: 42 + 28 + 56 + 56 rows (Kummer p₂, split Dickson p₁, torus T31, T24; tower and nulls; $m = 2, 3, 4$; $N <= 18$), *0 differences*; Kummer p₂ at `Fp64`, `Fp52`, `FpWide<4>`: 22 rows each, 0 differences],
  [G3], [arms of every speed cell agree on every field but the clock],
  [G4], [`pkm_tower_check verify`: 97/97 u64-builder rows agree; wide `--verify`: 82/82 rows ok (refutation ⇔ no decomposition; every decomposition and planted point satisfies every input equation and basis element in `BigUint`)],
)

= Speed (pre-registered)

7 interleaved repetitions, one thread, one core; F4 time per system.

#table(columns: (auto, 1fr, auto, auto, auto, auto, auto), stroke: 0.4pt + gray,
  [*id*], [*claim*], [*base ms*], [*cand. ms*], [*ratio*], [*min–max*], [*verdict*],
  [S1a], [split Dickson p₁, $m=2, N=16$: `Fp32` vs `f4_fp_tower` (≥ 0.95)], [3005], [2791], [1.091], [0.900–1.123], [PASS],
  [S1b], [torus T31, $m=2, N=16$: same (≥ 0.95)], [2988], [2784], [1.058], [0.969–1.096], [PASS],
  [S2a], [Kummer p₂ > $2^31$, $m=2, N=16$ (≥ 1.5)], [8921], [4009], [2.175], [2.034–2.306], [PASS],
  [S2b], [Kummer p₂, $m=3, N=12$ (≥ 1.5)], [15578], [7001], [2.240], [2.114–2.270], [PASS],
  [S3], [Kummer 45-bit, $m=2, N=16$: `Fp52` vs `Fp64` (≥ 1.8)], [10333], [4074], [2.568], [2.342–2.653], [PASS],
)

#figure(image("fig_claims.svg", width: 85%), caption: [Paired ratios (dot: median, bar: min–max, tick: threshold).])

#figure(image("fig_cost.svg", width: 85%), caption: [ns per $F_p$ multiply-add in the engine: 0.48 (`Fp32`, 31-bit), 0.71 (`Fp32`, 32-bit; `f4_fp_tower` 1.58), 0.72 (`Fp52`, 45-bit; `Fp64` 1.83), 2.0 (`Fp64`, 61-bit), 8.5 (P-192), 14.5 (P-224) — 30× the 31-bit cost.])

A 4-limb IFMA accumulator for P-224/P-256 was tried (exploratory) and dropped:
correct, but 16.9 vs 12.4 ns per entry against `FpWide<4>` on this host.

= Measured cells

Two random targets per cell; `D / width / steps`; reference: round 8 (Kummer p₂ $= 3 dot 2^30 + 1$; torus T31) and round 2 (Kummer p₁, $m = 4$).

#table(columns: 7, stroke: 0.4pt + gray,
  [*$m$*], [*$N$*], [*recorded*], [*$v_2=40$*], [*$v_2=48$*], [*$v_2=56$*], [*P-224*],
  [2], [10], [4/337/10], [4/337/10], [4/337/10], [4/337/10], [4/337/10],
  [2], [12], [5/1231/13], [5/1231/13], [5/1231/13], [5/1231/13], [5/1231/13],
  [2], [14], [5/2304/15], [5/2304/15], [5/2304/15], [5/2304/15], [5/2304/15],
  [2], [16], [5/4789/17], [5/4789/17], [5/4789/17], [5/4789/17], [5/4789/17],
  [2], [18], [5/8644/20], [5/8644/20], [5/8644/20], [5/8644/20], [5/8644/20],
  [3], [9], [6/1186/16], [6/1186/16], [6/1186/16], [6/1186/16], [6/1186/16],
  [3], [12], [6/5019/18], [6/5019/18], [6/5019/18], [6/5019/18], [6/5019/18],
  [4], [8], [6/1049/21], [6/1049/21], [6/1049/21], [6/1049/21], [6/1049/21],
)
$N = 2 dots 8$ ($m = 2$) and $N = 6$ ($m = 3$) also equal (README §4).

#table(columns: 5, stroke: 0.4pt + gray,
  [*$m$*], [*$N$*], [*torus T31 (recorded)*], [*P-192*], [*P-256*],
  [2], [12], [5/1223/13], [5/1223/13], [5/1223/13],
  [2], [16], [5/4786/17], [5/4786/17], [5/4786/17],
  [2], [18], [5/8639/20], [5/8639/20], [—],
  [3], [9], [6/1184/16], [6/1184/16], [6/1184/16],
  [3], [12], [6/5018/18], [6/5018/18], [—],
)

*Extension cells* (one target each; PROTOCOL §8 B), all equal to the recorded $p_1$ values:
Kummer $m=2, N=20$ at $v_2 = 40, 48, 56$: 6/33446/24 (4.82·10¹¹ multiply-adds; 166, 216, 614 s);
$m=3, N=15$ at $v_2=40$: 7/31397/24 (3.19·10¹², 1330 s);
$m=4, N=12$ at $v_2=40$ and P-224: 7/9452/29 (5.16·10¹⁰; 16 s and 383 s).
Kummer $m=2, N=20$ at P-224: 6/33446/24 (4.82·10¹¹, 3690 s, `--verify` ok). Torus $m=2, N=20$ at P-192: 6/33426/24 (4.79·10¹¹, 1905 s, `--verify` ok). The $N = 20$ Kummer multiply-add count is identical at all four Kummer primes (481,961,879,776).

#figure(image("fig_d_vs_n.svg", width: 95%), caption: [Solving degree against $N$; line: recorded small-prime $D$.])
#figure(image("fig_time_vs_n.svg", width: 80%), caption: [F4 time per system against $N$, $m = 2$ (two threads, shared host).])

= Reading and disclosures

H_same is supported in every series (PROTOCOL §6): no valued cell differs in
$D$; width and steps agree too, and on the fixed binary so does the
multiply-add count. The Kummer system's coefficients do not involve $zeta$, so
$v_2(p-1)$ can enter only through accidental cancellation mod $p$ — not seen.

Deviations: the first depth run sized S-row blocks by accumulator width, so
multiply-adds depended on the backend; fixed and re-run (v1 vs v2: 0
differences in $D$, width, steps over 197 rows). The `FpWide<4>` equality arm
on the old binary showed 4 multiply-add differences; re-run on the fixed
binary, 0. A container restart interrupted the extension run, which was
resumed. Conductor was unavailable.
