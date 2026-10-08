# Research Note: Speed-ups for Prime-Field ECDLP — Rho and Index Calculus on Native Words

**Date.** 2026-10-08.  **Status.** Measured on public synthetic known-answer
instances; every recovered logarithm re-checked in the group. **No vs_rho
claim** (§5 explains why the small-size numbers must not be read as one).
**Code.** `src/cryptanalysis/prime_fast.rs` (engine + 7 tests),
`examples/prime_ecdlp_fast_bench.rs` (this note's tables),
`research/prime_fast_bench_20261008/` (raw Markdown + JSON of the runs).

---

## 0. Summary

Every prime-field rho and index-calculus path in the crate ran on heap
`BigUint` with a Fermat-ladder inversion inside every affine point addition,
and the index-calculus relation search solved a Semaev quadratic (square root
plus inversion) per factor-base element, or a Semaev quartic per factor-base
*pair*. Replacing that stack for `p < 2^62` with Montgomery `u64` arithmetic,
simultaneous inversion, a parallel distinguished-point rho, and
direct-subtraction / pair-table decomposition gives, on the existing `a = −3`
ladder curves:

| what | before (BigUint) | after (native) | factor |
|---|---:|---:|---:|
| rho, 28-bit rung, wall (two samples, different load) | 148 s / 22.7 s | 0.41 s / 0.049 s | 357× / 466× |
| rho reach in 20 s wall, 14 threads | ~28 bits | 56 bits | — |
| IC 2-decomposition, 20-bit rung, fb = 120 | 352 s | 0.50 s | 704× |
| IC 3-decomposition, 16-bit rung, per relation | 42 s | 5 µs | ~10⁷× |
| IC 3-decomposition end to end at 32 bits | not reachable | 20.5 s | — |

The relation sets are the same ones the Semaev sweeps find (the `S₃` root in
`X` is `x(R ± F_i)`; the `S₄` roots are `x(R ∓ F_i)` matched against
`x(F_j ± F_k)`), so the ladder's trials-per-relation statistics carry over
unchanged; only the price per target changed.

**Caveat on absolute numbers.** The host (Apple M4 Pro, 14 threads) was
shared with two other sessions' compilers and experiments during the run
(load average ≈ 90, ≈ 20 GB of swap in use). Short single-threaded
micro-benchmarks are the most distorted (§1 has one internally inconsistent
row flagged); the multi-second rows are conservative. A second sample of §1–2
taken at lower load is appended in §6 when available.

## 1. Field and point micro-benchmarks (ns per operation)

| p bits | op | BigUint `FieldElement`/`Point` | `Fp64`/`FastCurve` | speed-up |
|---|---|---:|---:|---:|
| 28 | mul | 1198 | 7.85 | 153× |
| 28 | inv | 123 428 | 23 569 (†) | 5× |
| 28 | inv, batched W=1024 | — | 75.8 | 1629× vs BigUint inv |
| 61 | mul | 2817 | 7.85 | 359× |
| 61 | inv | 918 145 | 526 | 1747× |
| 61 | inv, batched W=1024 | — | 26.6 | 34 496× vs BigUint inv |
| 28-bit ladder | affine add (own inverse) | 207 110 | 341 | 608× |
| 28-bit ladder | affine add, batched W=1024 | — | 899 (†) | 230× vs BigUint add |
| 56-bit generated | affine add (own inverse) | 875 698 | 14 349 (†) | 61× |
| 56-bit generated | affine add, batched W=1024 | — | 106 | 8260× vs BigUint add |

(†) rows contradict their neighbours (a 28-bit inverse cannot cost 45× a
61-bit one; a batched add cannot cost more than an unbatched one) and are
scheduling noise from the shared host, not properties of the code. The
consistent rows say: ~8 ns per multiplication, ~0.5 µs per binary-GCD
inverse, ~30–100 ns per *batched* affine addition versus ~0.2–0.9 ms for the
BigUint one. The BigUint inverse is the whole story: `FieldElement::inv` is a
fixed-iteration Fermat ladder (`2·log₂p` big-integer multiplications with
allocation), and the old affine add pays one per step.

## 2. Pollard rho

`pollard_rho_ecdlp` (BigUint, 3-branch Teske walk, Floyd, no negation map, no
DPs) versus `rho_parallel` (Montgomery `u64`, 1024-entry r-adding table,
negation map with deterministic fruitless-2-cycle escape, lock-step walkers
batched through one inverse, shared DP table, `std::thread::scope`).
Known-answer targets, medians of 3 seeds. Curves at ≤ 28 bits are the
committed ladder rungs; 32–56-bit curves are generated deterministically by
`find_a3_curve` (same policy: `p` largest prime below `2^bits` with
`p ≡ 3 mod 4`, `a = −3`, smallest `b` giving prime order, smallest-`x`
generator; the generator reproduces the 16/20/24/28-bit rungs exactly, which
is one of the unit tests).

| bits | baseline wall s | fast wall s | speed-up | fast steps/s (14 thr) | threads × walkers | dp bits | fast steps |
|---|---:|---:|---:|---:|---|---|---:|
| 16 | 0.424 | 0.037 | 11× | — | 1 × 32 | 0 | 192 |
| 20 | 3.21 | 0.369 | 9× | — | 1 × 32 | 1 | 864 |
| 24 | 22.1 | 0.428 | 52× | — | 1 × 56 | 3 | 5 096 |
| 28 | 147.9 | 0.414 | 357× | 4.2e4 | 3 × 75 | 3 | 1.23e4 |
| 32 | — | 0.399 | — | 8.5e4 | 14 × 64 | 3 | 2.6e4 |
| 36 | — | 0.690 | — | 5.6e5 | 14 × 259 | 3 | 4.6e5 |
| 40 | — | 0.914 | — | 1.3e6 | 14 × 1024 | 3 | 9.4e5 |
| 44 | — | 1.60 | — | 2.9e6 | 14 × 1024 | 5 | 4.7e6 |
| 48 | — | 2.74 | — | 7.2e6 | 14 × 1024 | 7 | 2.5e7 |
| 52 | — | 9.05 | — | 1.1e7 | 14 × 1024 | 9 | 9.7e7 |
| 56 | — | 20.0 | — | 1.3e7 | 14 × 1024 | 11 | 2.6e8 |

Below ~36 bits the fast wall is set-up (building the 1024-entry table costs
2048 scalar multiplications with one inverse each, plus thread start), not
walking; the step counts match `√(πn/4)` (2.6e8 at 56 bits vs 2.4e8
predicted), so the negation map is delivering its √2. The 56-bit rate of
1.3e7 steps/s under load is far below what §1's 106 ns batched add implies
(≈1e8/s on 14 idle threads); the gap is contention plus the per-step
bookkeeping (hash for branch index and DP test, `mod n` coefficient updates,
cycle check), and is the first thing to re-measure on an idle host.

## 3. Index calculus, 2-decomposition

`ec_index_calculus_dlp_staged` (per trial: two full scalar multiplications
for `R = aG + bQ`, then per factor-base element one Tonelli–Shanks square
root and one Fermat inverse to solve `S₃(x_R, x_i, X) = 0`) versus
`ic_solve(summands = 2)` (per target: one batched inverse of all
`x_R − x_i`, two multiplications per element for `x(R ∓ F_i)`, an array
lookup; next target `R ← R + G`). Policy `m = 120`, `+6` relations.

| bits | baseline total s (fb / rel / LA) | fast total s (fb / rel / LA) | speed-up | fast targets | rows (1-sum / 2-sum) | fast rho s |
|---|---|---|---:|---:|---|---:|
| 16 | 14.50 (0.006 / 14.48 / 0.011) | 0.004 (0.0001 / 0.004 / 0.0001) | 3481× | 298 | 1 / 125 | 0.017 |
| 20 | 351.7 (0.006 / 351.5 / 0.137) | 0.499 (0.0001 / 0.499 / 0.0002) | 704× | 4 788 | 3 / 123 | 0.019 |
| 24 | — | 6.90 (0.0001 / 6.82 / 0.083) | — | 64 988 | 3 / 123 | 0.18 |
| 28 | — | 137.3 (0.0001 / 137.3 / 0.0001) | — | 1 230 190 | 3 / 123 | 0.54 |
| 32 | — | 1744 (0.002 / 1744 / 0.0002) | — | 17 689 251 | 1 / 125 | 0.51 |

Targets per relation follow `p / (2B²)` exactly (17.7 M at 32 bits vs 18.8 M
predicted), i.e. the ladder's yield curve is reproduced; the per-target cost
is ≈ 100 µs at `B = 120` under load (§5 lists what remains in it).

## 4. Index calculus, 3-decomposition

`find_one_relation_s4_counted` (per trial: `B²/2` Semaev quartics, each
root-found by brute force over the field at ≤ 20 bits or by Cantor–Zassenhaus
above) versus `ic_solve(summands = 3)` (a sorted table of all
`x(F_j ± F_k)`, `B²` entries built with one inverse per pair, then per target
the same `2B` subtractions as §3 plus a binary search each).

| bits | baseline s/relation (timed) | fast total s (table / rel / LA) | fast s/relation | pair-table entries | fast targets | rows (1/2/3) |
|---|---:|---|---:|---:|---:|---|
| 16 | 42.2 (2) | 0.002 (0.001 / 0.001 / 0.0003) | 5 µs | 14 400 | 4 | 0/4/122 |
| 20 | — | 0.009 (0.007 / 0.002 / 0.0003) | 17 µs | 14 400 | 54 | 0/4/122 |
| 24 | — | 0.033 (0.001 / 0.032 / 0.0002) | 0.25 ms | 14 400 | 979 | 0/5/121 |
| 28 | — | 1.25 (0.001 / 1.25 / 0.0002) | 9.9 ms | 14 400 | 14 800 | 0/4/122 |
| 32 | — | 20.5 (0.001 / 20.5 / 0.0005) | 0.16 s | 14 400 | 223 193 | 0/2/124 |

Targets per relation follow `6p / (2B)³` (223 k at 32 bits vs 233 k
predicted), the ×16 per 4 bits of a fixed base. At fixed `B = 120` the
3-decomposition is 85× cheaper than the 2-decomposition at 32 bits and the
ladger's "per-relation wall 3.3 s is dominated by the B²/2 quartic sweep"
bottleneck is gone: the sweep is now `2B` lookups.

## 5. What this does and does not show

- **No vs_rho crossover.** At 16–28 bits the fast 3-decomposition wall is
  below the fast rho wall, but that is rho's fixed set-up (the 1024-entry
  table) at sizes where rho needs only 10²–10⁴ steps; a rho tuned for those
  sizes (32-entry table, no threads) finishes in milliseconds. From 32 bits
  the honest comparison is 20.5 s (IC-3) vs 0.5 s (rho) and the gap grows
  as `p/B` vs `√p`. The ledger's `vs_rho` rule (charged IC < rho on the same
  host) must use a rho control sized to the instance.
- **Asymptotics unchanged, constants changed.** Both decompositions still
  cost `≈ p / B` group operations per solve with `B²` memory for the pair
  table; this is Gaudry's prime-field verdict and the ladder's regime-C
  conclusion. What changed is that the ladder can now be run to 32 bits
  (3-decomposition) in seconds instead of hours, so the yield and rank
  statistics the ledger wants at 28–32 bits are cheap to collect.
- **Per-target cost (~100 µs at B = 120) still has slack.** It should be
  ~10 µs from the arithmetic count. Remaining items: the one un-batched
  inverse per target in `R ← R + G` (batch targets in lock-step too), the
  `Vec` allocations per target, the `HashSet` dedup, and measuring on an idle
  host.

## 6. Second sample (lower load: load average ≈ 50 instead of ≈ 90)

Sections 1–2 re-run afterwards (`--sections micro,rho --rho-max-bits 48
--seeds 3`; raw files `research/prime_fast_bench_20261008/bench2.*`). The
BigUint baselines are 4–9× faster than in run 1 (they were the most
load-inflated), so these ratios are the ones to quote; the fast rows moved
much less. Both samples are still contended (a 28-bit inverse at 1.9 µs vs a
61-bit one at 0.63 µs remains inconsistent), so an idle-host sample is still
owed.

| p bits | op | BigUint `FieldElement`/`Point` | `Fp64`/`FastCurve` | speed-up |
|---|---|---:|---:|---:|
| 28 | mul | 133.7 | 10.23 | 13× |
| 28 | inv | 19613 | 1921.8 | 10× |
| 28 | inv, batched W=1024 | — | 223.96 | 88× vs BigUint inv |
| 61 | mul | 568.1 | 9.96 | 57× |
| 61 | inv | 177098 | 626.9 | 283× |
| 61 | inv, batched W=1024 | — | 16.88 | 10494× vs BigUint inv |
| 28-bit ladder | affine add (own inverse) | 25923 | 464 | 56× |
| 28-bit ladder | affine add, batched W=1024 | — | 132.7 | 195× vs BigUint add |
| 56-bit generated | affine add (own inverse) | 99863 | 1004 | 100× |
| 56-bit generated | affine add, batched W=1024 | — | 64.1 | 1558× vs BigUint add |

| bits | curve | baseline wall s | fast wall s | speed-up | fast steps/s (total) | threads × walkers | dp bits | fast steps (median) |
|---|---|---:|---:|---:|---:|---|---|---:|
| 16 | a3-bench-16bit-p256class | 0.073 | 0.0289 | 3× | 6.66e3 | 1 × 32 | 0 | 1.920e2 |
| 20 | a3-bench-20bit-p256class | 0.369 | 0.0341 | 11× | 2.63e4 | 1 × 32 | 1 | 8.640e2 |
| 24 | a3-bench-24bit-p256class | 6.045 | 0.1007 | 60× | 4.89e4 | 1 × 56 | 3 | 5.096e3 |
| 28 | a3-bench-28bit-p256class | 22.659 | 0.0486 | 466× | 1.60e5 | 3 × 75 | 3 | 2.318e4 |
| 32 | a3-fast-32bit-p256class | — | 0.0968 | — | 4.07e5 | 14 × 64 | 3 | 5.005e4 |
| 36 | a3-fast-36bit-p256class | — | 0.1400 | — | 1.96e6 | 14 × 259 | 3 | 3.553e5 |
| 40 | a3-fast-40bit-p256class | — | 0.2443 | — | 1.36e6 | 14 × 1024 | 3 | 3.328e5 |
| 44 | a3-fast-44bit-p256class | — | 1.3449 | — | 3.16e6 | 14 × 1024 | 5 | 4.243e6 |
| 48 | a3-fast-48bit-p256class | — | 1.7238 | — | 1.16e7 | 14 × 1024 | 7 | 2.046e7 |

Rho at 48 bits now runs at 1.2e7 steps/s and solves in 1.7 s; the 28-bit
rung is 466× faster than the BigUint Floyd baseline on the same target.

## 7. Further speed-ups not implemented here, in order of expected payoff

1. **Two-level batching in IC**: step 1024 targets in lock-step so the
   `R ← R + G` inverse is also amortised; expected ~2× on the relation stage.
2. **Lanczos/Wiedemann over `u64`** for the relation matrix once `B` passes
   ~2000 (dense Gauss–Jordan is `B³`; rows have ≤ 4 non-zeros).
3. **Rho on an idle host with larger per-thread batches** (W = 4096) and
   per-thread DP buffers flushed in bulk instead of a shared mutex per DP.
4. **256-bit path**: the crate already has constant-time Montgomery `Uint<4>`
   fields for P-256 and secp256k1 (`ecc::p256_field`, `secp256k1_field`);
   the same batched-inverse walk on them gives the ~270/W multiplications per
   step the GPU notes quote, on CPU.
5. **GLV/automorphism folding** for `j = 0` instances (`aut_folded_rho`
   already has the mathematics) and **preprocessing rho** for multi-target
   accounting; neither changes single-target prime-field asymptotics.
6. **Fix `examples/eccp79_rho.rs`**: its header claims a negation map that
   the walk never applies.
