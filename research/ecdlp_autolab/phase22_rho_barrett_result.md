# Phase 22.6 Barrett → rho integration: empirical result

**Date:** 2026-06-09
**Status:** COMPLETED — end-to-end speedup measured, arithmetic verified
**Program:** `phase22_rho_barrett.c`

## What was built

The v7 roadmap flagged two engineering wins that had been measured in
isolation but never integrated into the affine Pollard rho walk:

1. **Barrett reduction** for the fixed 81-bit benchmark prime. Phase 22.1
   proved 9.47× on `mulmod` in isolation, but `phase21_rho_v3.c` still ran
   on GMP `mpz_mul + mpz_mod`.
2. **Batched modular inversion** (Montgomery's trick). The affine walk pays
   one `mpz_invert` per step — the open caveat #2 in
   `phase22_barrett_result.md`. The inversion dominates the per-add cost
   (~80 muls via Fermat).

`phase22_rho_barrett.c` fuses both. Each thread runs `W` independent walks.
Per step it gathers `W` denominators and inverts them with a single field
inversion plus ~3W multiplies (Montgomery's trick). All field arithmetic is
`__uint128_t` Barrett — no GMP in the hot loop. This is the CPU embodiment of
the GPU strategy (many independent walks + special-form reduction + one
shared inversion per batch).

## Correctness

Two self-checks run before every benchmark and must pass or the program
aborts:

1. **Batched inverse vs GMP.** 64 random field elements inverted via
   Montgomery's trick; each compared against `mpz_invert`. Match.
2. **Walk trajectory vs GMP.** A single Barrett lane is run for 2000 steps
   and compared point-for-point against an independent GMP reference walk
   (same partition function, same negation map). Final points match exactly.

So the speedup is not from changing the walk — the trajectory is identical to
the v3 reference; only the arithmetic underneath is faster.

## Result (67.a1 benchmark prime, 81-bit, 4-core machine)

All runs on the same host (`nproc = 4`), curve
`p = 1208925819614629469615699, a = 12`, base point from `QUICK_START.md`.

| Implementation | Threads | Throughput | Per-core |
|----------------|--------:|-----------:|---------:|
| v3 (affine + GMP) baseline | 4 | 6.14e6 ops/s | 1.54e6 |
| **Barrett + batched inverse** | 4 | **1.87e7 ops/s** | **4.67e6** |
| v3 baseline | 1 | 1.75e6 ops/s | 1.75e6 |
| **Barrett + batched inverse** | 1 | **5.25e6 ops/s** | **5.25e6** |

**End-to-end speedup: ~3.0×** (4-thread: 1.87e7 / 6.14e6 = 3.04×;
single-core: 5.25e6 / 1.75e6 = 3.00×).

### Batch-size sensitivity

Single-core throughput vs `W` (walks per thread):

| W | ops/sec |
|---:|--------:|
| 64 | 5.25e6 |
| 256 | 4.81e6 |
| 1024 | 4.36e6 |
| 4096 | 4.12e6 |

Small batches win: `W ≈ 64–256` keeps the working set (W points + scratch)
in L1/L2. Large batches amortize the single inversion better in theory but
thrash cache in practice, so the inversion saving is already saturated by
W=64. Recommended default: `W = 256` (4-thread optimum, 1.87e7).

## Honest placement against the budget

The v6 gap to 80-bit feasibility was 5.5×. This delivers 3.0×, so the gap
narrows to ~1.8×. That is **below** the Phase 22.1 projection of 5–7×
end-to-end and above the pessimistic 1.28× floor.

Why 3.0× and not the projected 5–7×: the projection assumed inversion was
~25% of the add cost and Barrett would cut the remaining 75% by 9×. In
practice (a) Montgomery's trick removes most of the inversion cost but not
the per-lane bookkeeping, and (b) the `__uint128_t` Barrett `mulmod` in the
full point-arithmetic context (with the carry-handling in `mul_128`) runs
slower than the 9.4 ns microbenchmark because of register pressure and the
batched-inversion memory traffic. 3.0× is the measured, verified number.

## Remaining gap to 80-bit feasibility (~1.8×)

Still on the table from the v7 plan, in rough order of expected return:

1. **Montgomery-form arithmetic** to drop the conditional subtractions in
   Barrett and the `mul_128` carry path.
2. **SIMD (AVX2)** across lanes — the batched walks are already SIMD-shaped,
   so this is a natural fit (Phase 22.2).
3. **Multi-target rho** (√t speedup for t simultaneous targets), orthogonal
   to this work and multiplies directly on top of the 3.0×.
4. **GPU port** — the batched-walk structure here maps directly onto one
   thread per walk; this CPU result is the reference implementation.

## Reproduce

```bash
cd research/ecdlp_autolab
gcc -O3 -o phase22_rho_barrett phase22_rho_barrett.c -lgmp -lpthread
./phase22_rho_barrett 1208925819614629469615699 12 5 \
  942402533238732514353214 807227645196494903606326 100000 4 256
# self-check OK ... ~1.87e7 ops/sec on 4 threads
```
