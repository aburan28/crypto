# Protocol: the F4 tower engine over a generic prime field (2026-10-15)

Frozen before any measured cell. Calibration runs made before this file
(scratch, not committed, disclosed in §7) set the thresholds.

## 1. What is built

- `src/cryptanalysis/f4_tower_generic.rs`: `f4_tower_generic<F: PrimeField>`,
  the algorithm of `f4_fp_tower::f4_tower` with coefficients `F::Elem` and the
  dense kernels on `F::Wide` lazy accumulators (`row_load` / `row_sub` /
  `row_reduce` / `row_reduce_scaled`, sparse `row_scatter_axpy`), with
  per-row lazy-capacity counting. The `u32` engine is unchanged.
- `src/cryptanalysis/fp_backend_ifma.rs`: `Fp52`, `p < 2^52`, `Fp64`'s
  elements with an AVX-512 IFMA row kernel (capacity 4095 products).
- `examples/pkm_tower_pilot.rs`: `--engine generic`, `--backend
  auto|fp32|fp52|fp64|wide1|wide3|p192|fp224|wide4`, `--p` accepting
  `P-192`, `P-224`, `P-256`, `proth:K` (smallest prime `c·2^K + 1`, `c` odd)
  or a decimal prime; `auto` picks `Fp32` below `2^32`, `Fp52` below `2^52`
  with IFMA, `Fp64` below `2^63`, `FpP192` at P-192, `FpWide<3>` to 192 bits,
  `Fp224` to 255 bits, `FpWide<4>` above. Above `2^63` (or with `--builder
  wide`) systems are built by `examples/pkm_tower_wide/wide.rs`, generic over
  the backend (Kummer and Dickson torus), with `--verify`.

## 2. Primes

| label | p | v₂(p − 1) | v₂(p + 1) | backend (auto) |
|:--|:--|--:|--:|:--|
| T24 | 14,680,063 | 1 | 21 | Fp32 |
| p₁ | 2,013,265,921 = 15·2^27 + 1 | 27 | 1 | Fp32 |
| T31 | 2,046,820,351 = 61·2^25 − 1 | 1 | 25 | Fp32 |
| p₂ | 3,221,225,473 = 3·2^30 + 1 | 30 | 1 | Fp32 |
| v40 | 29,686,813,949,953 = 27·2^40 + 1 (the backlog card's prime) | 40 | 1 | Fp52 |
| v48 | `proth:48`, smallest prime c·2^48 + 1 | 48 | 1 | Fp52 |
| v56 | `proth:56`, smallest prime c·2^56 + 1 | 56 | 1 | Fp64 |
| P-224 | 2^224 − 2^96 + 1 | 96 | 1 | Fp224 |
| P-192 | 2^192 − 2^64 − 1 | 1 | 64 | FpP192 |
| P-256 | 2^256 − 2^224 + 2^192 + 2^96 − 1 | 1 | 96 | FpWide<4> |

## 3. Correctness gates (must pass before any speed claim is read)

- G1 unit tests (`cargo test --release --lib -- f4_tower_generic
  fp_backend_ifma`): every backend that holds the prime (`Fp32`, scalar
  `Fp32`, `Fp64`, `FpWide<1>`, `FpWide<4>`, `Fp224`, `Fp52` IFMA and scalar,
  and capacity-1/2/3 wrappers) reproduces `f4_fp_tower` report and basis on
  random rings at 97, 786433, 2013265921, 3221225473, 4294967291; Fp52
  against `BigUint` at the full capacity; exact refutation and point count by
  brute force at v40 and P-224; backends agree above `2^32` (56-bit, 40/48-bit,
  P-224).
- G2 equality runs (`run_eq.sh`): every round-8 series both engines run with
  `N ≤ 18` (`m = 2`), `N ≤ 12` (`m = 3, 4`), two random targets, tower and
  nulls: `pkm_tower_check compare-rows` of the `u32` rows against the generic
  `Fp32` rows (engine label rewritten, `backend`, `builder` ignored besides the
  wall clock) must report **0 differences** (D, width, steps, multiply-adds,
  every counter). The Kummer p₂ series at `Fp64`, `Fp52`, `FpWide<4>` (same
  prime, `N ≤ 16`) must equal the `Fp32` rows the same way.
- G3 the speed arms' rows of a cell must agree on every field but the clock
  and the backend label (checked by the report tool).
- G4 measured cells: `pkm_tower_check verify` on every finished `u64`-builder
  row (exhaustive over `V^m`, independent code) and the wide builder's
  `--verify` on every finished wide row (exhaustive decomposition count by
  backend square roots, verdict agreement, every decomposition point and the
  planted point re-checked in `BigUint` against every input equation and every
  basis element) must pass; `planted_ok` must be true on every planted row.

## 4. Pre-registered speed claims

Harness `run_speed.sh BIN runs/speed 7`: one system per run (random target 0),
`RAYON_NUM_THREADS=1`, `taskset -c 3`, 4 GB address-space cap, 7 repetitions,
arms of a cell adjacent, list order reversed on even repetitions. Time = the
row's `ms` (F4 only). Statistic: median over repetitions of the paired ratio
`t(baseline) / t(candidate)`; PASS iff it is at least the threshold. Spread
reported as min–max of the paired ratios. The host is shared with another
lane (not isolated, L0).

| id | cell | baseline → candidate | threshold |
|:--|:--|:--|--:|
| S1a | C1: split Dickson, p₁, m = 2, N = 16 | f4_fp_tower → generic Fp32 | 0.95 (no regression) |
| S1b | C2: torus Dickson, T31, m = 2, N = 16 | f4_fp_tower → generic Fp32 | 0.95 (no regression) |
| S2a | C3: Kummer, p₂, m = 2, N = 16 | f4_fp_tower → generic Fp32 | 1.5 |
| S2b | C4: Kummer, p₂, m = 3, N = 12 | f4_fp_tower → generic Fp32 | 1.5 |
| S3 | C5: Kummer, v40, m = 2, N = 16 | generic Fp64 → generic Fp52 (IFMA) | 1.8 |

Descriptive arms (no threshold): C6 Kummer v56 `Fp64` m = 2 N = 16; C7
Kummer P-224 `Fp224` m = 2 N = 14; C8 torus P-192 `FpP192` m = 2 N = 14.
Reported as ns per multiply-add and as a multiple of the `Fp32` cost.

## 5. Measured cells (`run_depth.sh BIN runs/depth`)

Generic engine, full-rank exit on, degree bound `n + d + 6`, seed
`0x504B4D54`, two random targets per cell, `RAYON_NUM_THREADS=2`, 4 GB cap,
budget 5400 s per system.

1. Kummer at v40, v48, v56 (u64 builder): `m = 2`, `N = 2…18`; `m = 3`,
   `N = 6, 9, 12`; `m = 4`, `N = 8`.
2. Kummer at P-224 (wide builder, `--verify`): the same cells.
3. Dickson torus at P-192 (`--verify`): `m = 2`, `N = 2…18`, tower and null;
   `m = 3`, `N = 6…12`. At P-256: `m = 2`, `N = 2…16`; `m = 3`, `N = 6, 9`.
4. Planted runs (one planted target, `t ≤ 4`, `m = 2, 3`) at v40, v48, v56,
   P-224 and P-192 (torus).
5. Optional, run last and only if time remains: `m = 2`, `N = 20` at v40, one
   target.

A cell beyond reach (the budget or the 4 GB cap) is reported as stopped, with
the step it reached; the series stops there.

## 6. Hypothesis and decision rule

H_same: at fixed `(m, N)` the solving degree `D` of the Kummer tower system
does not depend on `v₂(p − 1)` or on the size of `p`, and the torus `D` at
P-192 / P-256 equals the split / torus `D` at 31-bit primes.

Reference values: the recorded rows of `research/pkm_dickson_torus_20261015/`
(Kummer at p₂, p₁, p₀; split Dickson at p₁; torus at T31 / T24) and
`research/pkm_tower_round2_20260925` values quoted there (m = 4).

- **Supported** in a series if every valued cell's `D` (both targets) equals
  the reference `D`.
- **Refuted** if any cell's `D` differs in any target (one difference is
  enough; the cell is reported).
- Width (`max_cols_to_solution`) and step count are compared and reported
  the same way but decide nothing (round 8 saw prime-dependent widths at
  24 bits from cancellations).
- Cells not run or stopped are "not valued", never counted for or against.

## 7. Calibration (scratch, before this file)

One random system each, single thread, host load 2.3–4.5: Kummer m = 2,
N = 16: p₁ u32 2.91 s / generic Fp32 2.47 s; p₂ u32 9.35 s / Fp32 4.33 s;
v40 Fp64 11.2 s / Fp52 4.45 s; v56 Fp64 11.9 s; P-224 N = 14 Fp224 7.96 s
(FpWide<4> 7.92 s); torus P-192 N = 14 FpP192 5.07 s (FpWide<3> 5.21 s). The
P-224 spare-bit and P-192 Solinas backends do not change the row kernel's
cost in the engine, so no claim compares them.
