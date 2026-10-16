# Protocol: variable-time prime-field backends (`fp_backend`)

Frozen before the measured run. Pilots 1–6 (exploratory: 3 repetitions,
8–10 ms per measurement, earlier versions of the code) chose the
implementation and the thresholds below. No verdict uses them.

## What is measured

The backends of `src/cryptanalysis/fp_backend.rs` (trait `PrimeField`):
`Fp32` (canonical `u32`, Barrett), `Fp64` (Montgomery `u64`, `p < 2^63`),
`FpWide<N>` (Montgomery CIOS over `N` limbs), `FpWide<N, true>` (carry-free
CIOS, needs a spare top bit; alias `Fp224`), `FpP192` (Solinas, NIST P-192
only), against the code they replace or are compared with:

- **f4_fp_tower u32 path**: verbatim copies, in `examples/fp_backend_bench.rs`,
  of `Fp::reduce`/`Fp::mul`, `lazy_axpy` with its AVX2/AVX-512 dispatch,
  `Kernel::axpy`/`Kernel::sub`, and the arithmetic and bitset work of
  `Acc::sub_tail` in `src/cryptanalysis/f4_fp_tower.rs` (the engine does not
  export them). The F4 inner-loop mix is reproduced as
  - `axpy_dense_L`: load a row of length `L` into a `u64` accumulator, 16
    `Kernel::sub` updates from 4 pivot rows, then reduce every entry (the
    `bred` construction in `eliminate`), per entry update;
  - `scatter_sparse`: 256 sparse rows of ≤ 64 entries over 8,192 columns
    subtracted into a dense accumulator with the bitset updates
    (`sub_tail`), per touched entry.
- **prime_fast::Fp64** (`p < 2^62` Montgomery): mul, inverse.
- **ct_bignum::MontgomeryContext<3/4>** (constant time): `mont_mul`,
  `mont_sqr` at P-192 (3 limbs), P-224 and P-256 (4 limbs).
- **num_bigint::BigUint**: `(a·b) % p` and `modpow(p − 2)`.

Primes: `p31 = 2013265921` (15·2^27 + 1, the PKM test prime p₁), `p32 =
4294967291`, `p61 = 2^61 − 1`, `p63 = 2^63 − 25`, NIST P-192, P-224, P-256,
and `p200` = the first prime above 2^199 + 0x1234567.

Operations (ns per operation): `mul_tp` (4,096 independent products),
`mul_lat` / `sqr_lat` (dependent chain), `add_tp`, `inv`, `batch_inv` (per
element of 4,096), `sqrt` (64 squares), `legendre`, `axpy_lazy_L` /
`axpy_eager_L` (the dense F4 row kernel through the lazy accumulator, or
reduced after every product; `L` = 10^4 and 10^5), `dot_lazy_1e4` /
`dot_eager_1e4`. Comparison arms for the algorithm choices: `inv_fermat`,
`inv_plain_binary`, `legendre_euler`, `sqrt_ts_plain`, and
`scatter_sparse_scalar` (per-entry `wide_mul_add`).

## Correctness gate (before any timing)

- `cargo test --release --lib fp_backend` passes: every operation of every
  backend against `BigUint` on edge values (0, 1, 2, p − 1, p − 2, (p ± 1)/2,
  2^(32k) ± 1) and random values, at p = 3, 5, 7, 13, 17, five 31/32-bit
  primes, five 61–63-bit primes (two with v₂(p − 1) = 40), P-192, a 190-bit
  prime, P-224, P-256, the 200-bit prime, 2^127 − 1, 2^128 − 159, one-limb
  primes up to 2^64 − 59 and a 512-bit prime (inverse and Jacobi); lazy rows,
  sparse scatter, `LazyRow` and dot products equal the eagerly reduced
  result including all-(p − 1) worst cases; the binary-GCD inverse never
  falls back.
- The benchmark asserts, before timing, that each backend's products and
  square roots match `BigUint`, that lazy rows equal eager rows, that the
  `Fp32` dense and sparse kernels produce exactly the f4_fp_tower copies'
  results, that `ct_bignum` agrees with `BigUint`, and that the comparison
  arms (Fermat, plain binary, Euler, plain Tonelli–Shanks) agree with the
  default methods.
- Debug-assertion build: the same tests also run with
  `debug-assertions = on` (the capacity asserts) once before the report.

## Measurement

- Command: `taskset -c 3 nice -n 0 prlimit --as=4294967296 --
  fp_backend_bench --reps 9 --target-ms 40 --out runs/final`, release
  build, one thread, at the commit that adds this protocol (recorded in
  `runs/final/meta.txt` with `uptime`, CPU model and `rustc -V`).
- Each arm's iteration count is calibrated once (≈ 40 ms per measurement),
  then frozen. Repetition `r` runs every arm once, in an order rotated by
  `7r` and reversed on odd `r` (interleaved arms). Host: 4 vCPU, shared
  with another lane, not isolated (wall clock, L0).
- Analysis: `fp_backend_report runs/final` (Rust). For a claim with
  baseline `B` and candidate `C`, the paired ratio of repetition `r` is
  `B_r / C_r` (ns/op; above 1 means `C` is faster). Verdict: **PASS** when
  the median of the 9 paired ratios reaches the threshold, else **FAIL**.
  The minimum and maximum paired ratios are reported as the spread, and
  per-arm medians with min/max. One run; no reruns to change a verdict.

## Pre-registered claims

| id | claim (ratio = baseline / candidate ≥ threshold) | threshold |
|:--|:--|--:|
| C1a, C1b | `Fp32` dense row kernel vs f4_fp_tower, p31, L = 10^4, 10^5 (no regression) | 0.95 |
| C2a, C2b | `Fp32` sparse scatter (`row_scatter_axpy`) vs `sub_tail`, p31 and p32 (no regression) | 0.95 |
| C3a, C3b | `Fp32` dense row kernel vs f4_fp_tower at p32 (eager there), L = 10^4, 10^5 | 2.0 |
| C4 | `Fp32` mul throughput vs f4_fp_tower `Fp::mul`, p31 | 0.95 |
| C5a, C5b | `Fp64` mul throughput and latency vs `prime_fast::Fp64`, p61 | 0.95 |
| C5c | `Fp64` inverse (one-word almost-Montgomery inverse) vs `prime_fast::Fp64`, p61 | 1.4 |
| C6a, C6b | `FpP192` mul throughput / latency vs ct_bignum<3>, P-192 | 1.5 / 2.0 |
| C6c, C6d | `FpWide<3>` mul throughput / latency vs ct_bignum<3>, P-192 | 1.0 / 1.25 |
| C6e, C6f | `FpWide<4, true>` mul throughput / latency vs ct_bignum<4>, P-224 | 1.3 / 1.8 |
| C6g, C6h | `FpWide<4>` mul throughput / latency vs ct_bignum<4>, P-256 | 1.1 / 1.3 |
| C7a, C7b | wide mul vs `BigUint` `(a·b) % p` at P-192, P-224 | 5.0 |
| C7c | `FpWide<4>` mul vs `BigUint` at P-256 | 4.0 |
| C7d, C7e | wide inverse vs `BigUint::modpow(p − 2)` at P-192, P-224 | 5.0 |
| C8a | `Fp32` lazy vs eager row update, p31, 10^4 | 5.0 |
| C8b | `Fp64` lazy vs eager row update, p61, 10^4 | 2.0 |
| C8c–g | `FpWide` lazy vs eager row update: P-192 (10^4, 10^5), P-224 spare (10^4, 10^5), P-256 (10^4) | 1.5 |
| C9a, C9b | lazy vs eager dot product, P-192, P-224 spare | 1.5 |
| C10a, C10b | binary-GCD inverse vs plain binary / Fermat, `FpWide<4>` at P-224 | 1.3 / 2.5 |
| C11a, C11b | binary Jacobi vs Euler's criterion, P-192 (`FpWide<3>`), P-224 (`FpWide<4>`) | 3.0 |
| C12 | windowed vs plain Tonelli–Shanks, `FpWide<4>` at P-224 (s = 96) | 3.0 |
| C13 | Solinas `FpP192` vs Montgomery `FpWide<3>` mul throughput, P-192 | 1.3 |
| C14 | `FpWide<4, true>` vs `FpWide<4>` mul throughput, P-224 | 1.15 |

The exact arm pairs are the `CLAIMS` table in `examples/fp_backend_report.rs`.
A FAIL is reported as such; a failed specialisation claim (C13, C14) means
the specialisation is not recommended in the module documentation.

## Out of scope

Constant-time behaviour; multi-threaded kernels; whole-F4 or whole-curve
end-to-end timings (lanes F and G measure those once they port); SIMD
(AVX-512 IFMA) wide-field kernels.
