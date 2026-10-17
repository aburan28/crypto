# Strong-rho canonicalization, round 2: protocol

Declared on 2026-10-16, before any change in this round was written or
measured. Parent: `72843c83a` (`main`, which carries #1999).

## Where the cost is

A line-table build of the parent was profiled under Callgrind on one n53
fixture (`runs/profile/`).
- `canonicalize` is still about 430 instructions per call, about 45% of a
  938-instruction walk step.
- The byte-table lookups are now about 3 instructions each.
- `MulMod::mul` (the orbit multiplier applied to `a` and `b`) spends about 27
  instructions per call on one line: `u64`→`f64` and `f64`→`u64` casts. Rust
  lowers those to multi-instruction sequences on baseline x86-64.
- `rotate` tests `k == 0` on every call.

## Candidate changes (no change to any value computed)

1. **`MulMod::mul`'s fast path** (modulus below 2^50, so both operands and the
   quotient are below 2^50): cast through `i64`. The quotient estimate is the
   same number, computed with single-instruction conversions.
2. **`NormalBasis::rotate`** drops the `k == 0` branch. For `v < 2^n` and
   `n < BITS`, the general expression already returns `v` at `k = 0`. A
   `debug_assert` pins the precondition.
3. **`canonicalize` reuses rotations.** Outside the tie branch, `x_k` is the
   least rotation already found and `y_k` is the rotation the sign test
   computed; they are no longer recomputed.
4. **The multiplier table also holds `−λ^k`**, so a negated step reads it
   instead of computing `0 − λ^k`.

## Identity gate

- Every output line of `koblitz_rho_fixture n 0 signed_frobenius 4 strong`,
  for n = 41, 53 and 61, is identical once `*_ms` and `peak_rss_bytes` are set
  aside.
- All strong-rho unit tests pass, with debug assertions on.

Any difference stops the round.

## Measurements

- **Primary: Callgrind instructions per walk step.** n53, 4 fixtures, one run
  per arm. Deterministic.
- **Native time:** ns per walk step at n = 53 and 61.
  - `isolated_bench run --cpus 3 --max-other-cpu 0.3`, `RAYON_NUM_THREADS=1`;
  - an A/A pass (baseline × 6), then 5 interleaved A/B rounds;
  - contended runs excluded and counted.
- **Rule.** A gain inside the A/A spread is not established. This VM's A/A
  spread was 11–18% in the last round, so the instruction count carries the
  result.
- **IPC.** `tools/perf_ipc.sh` is run once per arm to record whether the host
  can count. This VM has no PMU, so an `unsupported` record is the expected
  outcome; no IPC is claimed without counters.

## Addendum (written after changes 1–4 were measured, before change 5 was)

Changes 1–4 cut instructions per step from 938 to 850 at n53 (`runs/callgrind/`),
with identical outputs at n = 41, 53 and 61.

A `perf record -e cpu-clock` sample of that build (`runs/profile/`) still puts
the `MulMod` quotient line at about 9% of `canonicalize`'s time. Its float
multiply and conversions form one long dependency chain.

5. **Shoup's fixed-operand product for the orbit multiplier.**
   - Every multiplier is one of the 2n values `±λ^k`, so each gets a
     precomputed `w' = ⌊w·2^64/m⌋`.
   - `a·w mod m` is then one 64×64→128 high product, two low products and at
     most one subtraction: exact for `a, w < m < 2^63`.
   - Larger moduli keep the old path.
   - A new unit test checks it against the plain product on moduli up to and
     past `2^63`.

The final candidate carries changes 1–5. It is measured against the same
baseline binary, under the same identity gate and the same measurements.

## Addendum 2 (change 6; its Callgrind count was taken before this entry)

6. **`canonicalize` is `#[inline(always)]`** into `solve`'s three call sites.
   That removes six callee-saved register saves and restores per call, and
   passes the 40-byte walk state in registers instead of memory.
   - Deviation from the protocol: one Callgrind run of this build (2,369,298,652
     instructions at n53 over 4 fixtures, outputs identical) was taken before
     this entry was written.
   - It joins the final candidate because no value changes.
   - The native A/B below is the first time measurement of the final candidate
     that includes it.
