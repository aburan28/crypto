# CryptoPro-B chain sweep: every chain realisation and scalar-multiplication parameter of the GLV endomorphisms

A constructive measurement on GOST R 34.10-2001 **CryptoPro-B** (RFC 4357
§11.4.4; `p = 2^255 + 3225`, prime order `n`, CM by the maximal order of
discriminant −619, class number 5).  It follows
`research/cryptopro_b_glv_chain_20261005` (PR #1408), which measured GLV-2
with one chain realisation of one endomorphism at one setting, and sweeps
the parameters around it: **which chain** (22 orderings of 11 elements, every
realisation within 2× of the cheapest), **which chain evaluator** (PR #1408's
projective one and a Jacobian one), and **which scalar-multiplication
settings** (wNAF width 3–7, tables affine or Jacobian), for the baseline as
well as for GLV.

The hypothesis, arms, units, predictions, success and stop conditions and
commands were committed in [`PROTOCOL.md`](PROTOCOL.md) (commit `df397ae9`)
before either run, and that file is unchanged.  Everything below was
produced by one counting run and one timing run of
`cryptopro-b-chain-sweep` at that commit.

This is **not an attack**: it computes `k·P` faster on the same curve.  Class
**engineering** (§3).

## Verdict

| preregistered condition | result | |
|:--|:--|:--|
| 1. every correctness check passes in both runs | every baseline and GLV arm equals the textbook `BigUint` `k·P` on all 256 pairs, every `φ(P)` equals `λ·P` on all 256 points and reproduces its frozen test vectors, every decomposition is valid | **met** |
| 2. `φ(P)` counts equal the model exactly; totals within 1 % | all 44 `φ(P)` arms exact (`M`, `S`, no inversion, identical on every input); worst total 0.26 % from the model's expectation | **met** |
| 3. counted: cheapest GLV ≤ cheapest baseline / 1.35 | 1 993.7 against 2 807.6 M_eq: **1.408** | **met** |
| 4. wall time: the declared GLV configuration beats the fastest baseline by > 1.252 (median) and > 1.337 (min), and the continuity arm by more than the A/A noise floor | the fastest baseline by median is the reference itself (`jacobian/w5`); **1.417** (median), **1.497** (min); 1.131× the speed of PR #1408's configuration on medians (11.6 % less time) against a 0.7 % A/A floor | **met** |
| 5. the lowest `φ(P)` median is `4+1w/7.5.5` with the optimised evaluator | 3 095 ns, the lowest of the 44 stage arms (next: `4+1w/5.7.5` optimised, 3 350) | **met** |

## Results

The bench's own table, every arm a row, is [`sweep.stdout.md`](sweep.stdout.md)
(full per-round record: [`sweep.json`](sweep.json); counting run:
[`counts.stdout.md`](counts.stdout.md), [`counts.json`](counts.json)).  The
rows that answer the question, verbatim from it (counted M_eq = field
multiplications + squarings per scalar multiplication, mean over the 256
pairs; ns per scalar multiplication over 20 rounds; ratios against the
reference row):

| arm | class | counted M_eq per k·P | ratio reference/arm, counted | ns per k·P, median of 20 | min of 20 | ratio reference/arm, median | ratio reference/arm, min | correctness |
|:--|:--|--:|--:|--:|--:|--:|--:|:--|
| baseline: width-5 NAF, Jacobian table (the reference) | reference | 2 807.6 | 1.000 | 117 103 | 111 538 | 1.000 | 1.000 | yes |
| baseline: width-6 NAF, Jacobian table | baseline variant | 2 835.2 | 0.990 | 117 640 | 112 867 | 0.995 | 0.988 | yes |
| baseline: width-5 NAF, affine table (#1408's baseline) | baseline variant | 2 916.5 | 0.963 | 123 945 | 118 152 | 0.945 | 0.944 | yes |
| GLV-2, 4 + ω as 7·5·5, optimised evaluator, Jacobian tables, width 5 (declared in the protocol) | engineering | 1 993.7 | 1.408 | 82 670 | 74 489 | 1.417 | 1.497 | yes |
| GLV-2, same chain and evaluator, width 4 | engineering | 2 004.6 | 1.401 | 83 199 | 78 876 | 1.408 | 1.414 | yes |
| GLV-2, same chain, generic evaluator, Jacobian tables, width 5 | engineering | 2 041.7 | 1.375 | 85 381 | 77 414 | 1.372 | 1.441 | yes |
| GLV-2, same chain, optimised evaluator, affine tables, width 4 | engineering | 2 069.6 | 1.357 | 89 466 | 84 666 | 1.309 | 1.317 | yes |
| continuity: #1408's configuration (5·5·7, generic evaluator, affine tables, width 5) | engineering | 2 205.3 | 1.273 | 93 486 | 88 415 | 1.253 | 1.262 | yes |
| A/A control: the reference timed again | A/A control | 2 807.6 | 1.000 | 117 952 | 111 120 | 0.993 | 1.004 | yes |

**The best configuration found is the one the protocol declared**: `4 + ω`
evaluated as `7·5·5` with the Jacobian evaluator, Jacobian tables, width 5 —
1.408× the best baseline in counted operations and 1.417× (median) in wall
time.  It is also the fastest GLV arm of the 31 measured.

### The grid on the cheapest chain, and the baseline grid

Counted M_eq / median µs per scalar multiplication (a view of the same
table):

| GLV-2 via 4 + ω (7·5·5) | w = 3 | w = 4 | w = 5 | w = 6 | w = 7 |
|:--|--:|--:|--:|--:|--:|
| optimised evaluator, jacobian tables | 2 146 / 89.6 | 2 005 / 83.2 | 1 994 / 82.7 | 2 151 / 87.5 | 2 585 / 104.2 |
| optimised evaluator, affine tables | 2 119 / 92.4 | 2 070 / 89.5 | 2 157 / 91.4 | 2 456 / 101.1 | 3 138 / 124.5 |
| generic evaluator, jacobian tables | 2 194 / 90.0 | 2 053 / 85.7 | 2 042 / 85.4 | 2 199 / 91.1 | 2 633 / 106.6 |
| generic evaluator, affine tables | 2 167 / 94.4 | 2 118 / 91.6 | 2 205 / 94.2 | 2 504 / 103.9 | 3 186 / 127.0 |
| baseline, jacobian table | 3 056 / 126.4 | 2 886 / 119.2 | 2 808 / 117.1 | 2 835 / 117.6 | 3 012 / 125.4 |
| baseline, affine table | 3 018 / 129.7 | 2 924 / 125.2 | 2 916 / 123.9 | 3 030 / 126.7 | 3 342 / 138.3 |

* The baseline's optimum is flat between widths 5 and 6 (2 808 / 2 835
  M_eq), GLV's between 4 and 5 (2 005 / 1 994): its scalars are half as
  long, so widths 6–7 pay for tables a 128-bit scalar cannot amortise.
* **Jacobian tables beat affine tables for both arms** at every width ≥ 4,
  in counts and in wall time: a Fermat inversion costs 264 M_eq in this
  field, more than the 5 M_eq a full addition costs over a mixed one
  (11M + 5S against 7M + 4S) summed over the ≈ 43 additions.  At width 3
  the counts favour affine tables slightly and the wall time does not.  GLV gains more
  from dropping it (its batched conversion covers sixteen points, the
  baseline's eight), which is why the tuned comparison is *more* favourable
  to GLV than PR #1408's, not less.
* The evaluator is worth exactly 48 M_eq at every width and table setting
  (the `generic` and `optimised` rows differ by the chain alone).  Its stage
  diagnostic is 1.4 µs (4.54 → 3.10 µs); the totals differ by 0.4–3.6 µs
  across the ten settings, which is that difference plus round-to-round
  noise.

### Every element at the declared configuration

| element | order | chain, counted M_eq | total, counted M_eq | ratio reference/arm, counted | ns, median | ns, min | ratio reference/arm, median |
|:--|:--|--:|--:|--:|--:|--:|--:|
| 4 + 1ω | 7·5·5 | 80 | 1 993.7 | 1.408 | 82 670 | 74 489 | 1.417 |
| 9 + 1ω | 7·5·7 | 95 | 2 009.2 | 1.397 | 83 784 | 80 180 | 1.398 |
| 15 + 2ω | 7·5·5·5 | 111 | 2 022.2 | 1.388 | 84 507 | 78 070 | 1.386 |
| 2 + 1ω | 23·7 | 112 | 2 023.8 | 1.387 | 83 761 | 79 139 | 1.398 |
| 0 + 1ω | 31·5 | 121 | 2 035.6 | 1.379 | 84 820 | 80 714 | 1.381 |
| 20 + 1ω | 23·5·5 | 128 | 2 041.8 | 1.375 | 84 089 | 80 804 | 1.393 |
| 5 + 1ω | 37·5 | 139 | 2 053.0 | 1.368 | 84 973 | 79 744 | 1.378 |
| 54 + 1ω | 5·5·5·5·5 | 136 | 2 053.3 | 1.367 | 85 397 | 79 768 | 1.371 |
| 25 + 1ω | 23·5·7 | 143 | 2 057.1 | 1.365 | 85 110 | 82 416 | 1.376 |
| 39 + 1ω | 7·5·7·7 | 141 | 2 058.5 | 1.364 | 84 482 | 80 684 | 1.386 |
| 37 + 3ω | 23·5·5·5 | 159 | 2 070.6 | 1.356 | 85 333 | 80 970 | 1.372 |

The chain choice moves the total by at most 77 M_eq (3.9 %); the measured
medians span 82.7–85.4 µs, in nearly the counted order (the gaps between
neighbours are within twice the A/A floor, so the order inside the pack is
not separately established).  `4 + ω` is the cheapest element by count and
by median.

### `φ(P)` alone, every chain and both evaluators (stage diagnostics)

| chain | element | order | generic: counted M + S | generic: ns median / min | optimised: counted M + S | optimised: ns median / min |
|:--|:--|:--|--:|--:|--:|--:|
| `4+1w/7.5.5` | 4 + 1ω | 7·5·5 | 124 + 4 | 4 536 / 4 346 | 78 + 2 | 3 095 / 2 852 |
| `4+1w/5.5.7` | 4 + 1ω | 5·5·7 | 124 + 4 | 4 563 / 4 353 | 87 + 2 | 3 381 / 3 253 |
| `4+1w/5.7.5` | 4 + 1ω | 5·7·5 | 124 + 4 | 4 533 / 4 371 | 87 + 2 | 3 350 / 3 213 |
| `9+1w/7.5.7` | 9 + 1ω | 7·5·7 | 139 + 4 | 5 093 / 4 833 | 93 + 2 | 3 623 / 3 353 |
| `9+1w/7.7.5` | 9 + 1ω | 7·7·5 | 139 + 4 | 5 075 / 4 776 | 93 + 2 | 3 637 / 3 350 |
| `9+1w/5.7.7` | 9 + 1ω | 5·7·7 | 139 + 4 | 5 093 / 4 815 | 102 + 2 | 3 932 / 3 782 |
| `15+2w/7.5.5.5` | 15 + 2ω | 7·5·5·5 | 159 + 5 | 5 892 / 5 426 | 108 + 3 | 4 221 / 4 026 |
| `2+1w/23.7` | 2 + 1ω | 23·7 | 224 + 3 | 8 310 / 7 998 | 111 + 1 | 4 539 / 4 345 |
| `15+2w/5.5.5.7` | 15 + 2ω | 5·5·5·7 | 159 + 5 | 5 828 / 5 472 | 117 + 3 | 4 521 / 4 343 |
| `15+2w/5.5.7.5` | 15 + 2ω | 5·5·7·5 | 159 + 5 | 5 877 / 5 569 | 117 + 3 | 4 518 / 4 368 |
| `15+2w/5.7.5.5` | 15 + 2ω | 5·7·5·5 | 159 + 5 | 5 810 / 5 447 | 117 + 3 | 4 520 / 4 197 |
| `0+1w/31.5` | 0 + 1ω | 31·5 | 269 + 3 | 9 978 / 9 609 | 120 + 1 | 4 948 / 4 784 |
| `20+1w/23.5.5` | 20 + 1ω | 23·5·5 | 244 + 4 | 8 978 / 8 256 | 126 + 2 | 5 099 / 4 615 |
| `54+1w/5.5.5.5.5` | 54 + 1ω | 5·5·5·5·5 | 179 + 6 | 6 622 / 5 755 | 132 + 4 | 5 093 / 4 715 |
| `5+1w/37.5` | 5 + 1ω | 37·5 | 314 + 3 | 11 633 / 11 059 | 138 + 1 | 5 774 / 5 517 |
| `39+1w/7.5.7.7` | 39 + 1ω | 7·5·7·7 | 189 + 5 | 6 920 / 6 522 | 138 + 3 | 5 372 / 4 960 |
| `39+1w/7.7.5.7` | 39 + 1ω | 7·7·5·7 | 189 + 5 | 6 942 / 6 507 | 138 + 3 | 5 392 / 4 958 |
| `39+1w/7.7.7.5` | 39 + 1ω | 7·7·7·5 | 189 + 5 | 7 019 / 6 502 | 138 + 3 | 5 345 / 4 669 |
| `25+1w/23.5.7` | 25 + 1ω | 23·5·7 | 259 + 4 | 9 521 / 8 910 | 141 + 2 | 5 631 / 5 259 |
| `25+1w/23.7.5` | 25 + 1ω | 23·7·5 | 259 + 4 | 9 543 / 9 051 | 141 + 2 | 5 658 / 5 315 |
| `39+1w/5.7.7.7` | 39 + 1ω | 5·7·7·7 | 189 + 5 | 6 948 / 6 470 | 147 + 3 | 5 656 / 5 288 |
| `37+3w/23.5.5.5` | 37 + 3ω | 23·5·5·5 | 279 + 5 | 10 429 / 9 647 | 156 + 3 | 6 302 / 5 872 |

* Counted `M` and `S` equal the chain sweep's model for all 44 arms, and the
  medians follow the counts (Spearman rank correlation 0.976 over the 44
  arms; 35–42 ns per counted M_eq).
* Under the `generic` evaluator the orderings of one element cost the same
  (identical counts, medians within 1.5 %); under the `optimised` one the
  largest prime belongs first: `0+1w/31.5` costs 121 M_eq / 4.9 µs where PR
  #1408 measured `5·31` at 272 M_eq / 9.5 µs.
* The decomposition (`num-bigint` Babai rounding, no field operations, so
  absent from the counts) takes 1.07–1.25 µs for every element, ≈ 1.4 % of
  the GLV total.

### Where PR #1408's 1.25× went

PR #1408's configuration is the continuity arm here.  Against its own
baseline (`affine/w5`) it measures 1.326 (median) / 1.336 (min) in this run
(counted: 1.323); PR #1408 measured 1.252 / 1.337.  The minima agree; PR
#1408's median ratio sat lower in a noisier run (its A/A median ratio was
0.936, this run's is 0.993).  From that configuration to the best one, in
counted M_eq:

| step | GLV arm | baseline arm | ratio |
|:--|--:|--:|--:|
| PR #1408's configuration (`generic`, affine, w 5 against affine w 5) | 2 205.3 | 2 916.5 | 1.323 |
| tables kept Jacobian, both arms (w 5) | 2 041.7 | 2 807.6 | 1.375 |
| + the optimised chain evaluator on `7·5·5` | 1 993.7 | 2 807.6 | **1.408** |

## Model against measurement

The prediction file (`frozen_inputs/chainsweep_model.json`) gave every
total arm's expected counted M_eq from the same wNAF recoding and lattice
bases; the counting run is within 0.26 % of it everywhere and exact for
every chain.  Wall time per counted M_eq is 39.7–43.6 ns across the total
arms and 35.4–41.5 ns across the `φ(P)` arms; the counts leave out
additions and subtractions, recoding, table lookups and the decomposition,
which is what that spread is.  The counted headline (1.408) and the median
headline (1.417) agree to within 1 %.

## Host, build and isolation

* Host: `Intel(R) Xeon(R) Processor @ 2.10GHz`, a cloud VM with 4 logical
  CPUs (`pclmulqdq sse4_2 popcnt aes avx bmi1 avx2 bmi2 avx512f adx`),
  kernel `6.18.44-fc-v70`, Linux x86-64.  The bench's host line says
  `logical cpus: 1` because it ran pinned to CPU 3.
* Build: `rustc 1.97.0 (2d8144b78 2026-07-07)`, `release` profile, no
  `target-cpu` flags; the counting binary is the same source built with
  `--features cryptopro-b-opcount` into `target/opcount`.
* Commit at run time: `df397ae9794f8f464cbefa997564b6cb88149de6` for both
  runs.  Both records say `git_dirty: true` only because the runs' own output
  files were untracked at that moment; no tracked file differed from the
  commit, and the outputs were committed next.
* Isolation (§10): the timing run went through the native
  `isolated_bench run` on CPU 3 (`--max-other-cpu 0.5`, as declared), which
  admitted it at the first attempt and recorded
  ([`isolated_bench_record.jsonl`](isolated_bench_record.jsonl)):
  `contended: false`, other processes 0.42 CPU-seconds during the 25.8 s
  run (the agent runtime `claude[95]` 0.37 of it), 62 involuntary context
  switches, no major faults, PSI `some avg10` 0.16 before and 0.08 after, 53
  threads moved off the CPU.  The counting run is deterministic and was not
  isolated (counts do not depend on contention).
* Noise: A/A median ratio 0.993 and min ratio 1.004.  Per-round spread
  (max/min − 1 over the 20 rounds): reference 12.7 %, A/A 24.6 %, the
  declared GLV arm 39.1 %, the continuity arm 21.7 %, the cheapest `φ(P)`
  arm 43.2 %.  Medians are the robust statistic; the minima are reported as
  the least-contended ones and the declared arm's minimum ratio (1.497) is
  less reliable than its median ratio.
* One hardware class: Linux x86-64 on this VM.  Not measured on Arm64, other
  x86-64 hosts, other field arithmetic or a constant-time implementation.

## Honesty notes

* **Variable time everywhere**; a measurement harness, never for secret
  scalars (see the module docs).
* **The table conclusion depends on the inversion.**  This field inverts by
  Fermat (256 S + 8 M = 264 M_eq).  With an inversion several times cheaper
  (binary extended Euclid, safegcd) the batched affine conversion would cost
  less than the full-addition penalty and affine tables would win again for
  both arms.  The chain conclusions (which element, which order, which
  evaluator) do not depend on it.
* **Field arithmetic is the crate's generic `U256::mont_mul`** (no dedicated
  squaring, no assembly).  A faster field moves both arms; the decomposition
  (≈ 1.1 µs) does not scale with it.
* **The chain set is bounded by its generator**: chains outside it are not
  measured.  Their counted cost is bounded below by the chain sweep's
  certificate (≥ 173 M_eq against 80 for the cheapest), which is a property
  of the evaluator formulas this crate's counting build confirms on all 22
  measured chains; the bound itself was computed by the autoresearcher.
* **Counts are of field operations**, not instructions; additions and
  subtractions are counted in the JSON (`add`) but not in M_eq.

## Commands

Exactly what was run, in this order, at commit `df397ae9`:

```sh
cargo test --release --lib cryptopro_b                                    # 41 passed
cargo test --release --features cryptopro-b-opcount --lib cryptopro_b     # 43 passed
cargo build --release --bin cryptopro-b-chain-sweep --bin isolated_bench
CARGO_TARGET_DIR=target/opcount cargo build --release \
  --features cryptopro-b-opcount --bin cryptopro-b-chain-sweep

D=research/cryptopro_b_chain_sweep_20261006
SHA=683f92d0e56f5606f78b6c685f747a7b67106a7b5e5883e81869c8d6199d7336

target/opcount/release/cryptopro-b-chain-sweep --count-only \
  --chains $D/frozen_inputs/chains.constants.json --expect-sha256 $SHA \
  --pairs 256 --seed 20261006 --json $D/counts.json > $D/counts.stdout.md

target/release/isolated_bench run --wait --cpus 3 --settle 5 --max-other-cpu 0.5 \
  --label cryptopro-b-chain-sweep --out $D/isolated_bench_record.jsonl -- \
  target/release/cryptopro-b-chain-sweep \
  --chains $D/frozen_inputs/chains.constants.json --expect-sha256 $SHA \
  --pairs 256 --rounds 20 --seed 20261006 --counts $D/counts.json \
  --json $D/sweep.json > $D/sweep.stdout.md
```

## Files

| file | what |
|:--|:--|
| `PROTOCOL.md` | the protocol, committed before the runs |
| `README.md` | this record |
| `frozen_inputs/chains.constants.json`, `frozen_inputs/chainsweep_model.json`, `frozen_inputs/SHA256SUMS` | the frozen chain set and the prediction (Python provenance, generated in aburan28/crypto-autoresearcher) |
| `counts.json`, `counts.stdout.md` | the counting run |
| `sweep.json`, `sweep.stdout.md` | the timing run (every per-round time, host, build, correctness) |
| `isolated_bench_record.jsonl` | the isolation record of the timing run |
| `src/ecc/cryptopro_b_chain_set.rs` | the chain-set loader and its tests |
| `src/ecc/cryptopro_b_point.rs` | the evaluators, table modes and their tests |
| `src/ecc/cryptopro_b_field.rs` | the `cryptopro-b-opcount` counter |
| `src/bin/cryptopro_b_chain_sweep.rs` | the benchmark |
