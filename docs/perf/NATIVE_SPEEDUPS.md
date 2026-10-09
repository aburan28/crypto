# Native speedups: assembly, FFI, and what actually pays

A decision record and todo list for "should the hot arithmetic be
assembly or FFI?", written 2026-10-09 against the ecbench harness and
the ECDLP research code in `src/cryptanalysis/`.  Items marked **TODO**
are open; each names the constraint that decides how it may be done.

## 1. Decision: no assembly now, intrinsics only behind a profile

**Not yet.**  Three facts about this repository decide it:

1. **The research unit does not move with codegen.**  ecbench scores a
   method in counted group additions over `√r` (`docs/ecbench/README.md`
   §7).  Wall time and instruction counts are controls, not the result.
   Hand-tuned multiplication changes no frontier, bound or verdict; it
   only matters where throughput caps the trials one can afford (the
   ECC2K-130 campaign, isolab wall-clock arms).
2. **The cheap native wins were already in.**  secp256k1 and P-256 run
   Montgomery form on four 64-bit limbs with `u64×u64→u128` products;
   `Gf2` (`semaev_decomp.rs`) multiplies with `pclmulqdq` / `pmull`,
   folds the reduction, keeps squaring chains in a vector register and
   batches inversions with Montgomery's trick where method identity
   allows (`koblitz_fast::add_many*`, `koblitz_strong_rho`); SHA-256
   uses SHA-NI; the `cryptanalysis` C repository already ships Rust FFI
   bindings.  There is no large pure-codegen gap left in the field layer.
3. **Assembly costs replayability.**  Sealed records and host-class
   hashing assume one binary behaves identically everywhere; isolab
   runs x86_64 and aarch64.  A per-architecture kernel doubles what must
   be audited, and a bug in it invalidates sealed evidence.  The one
   real argument for assembly, that a compiler cannot reintroduce
   branches into constant-time code, applies to the library's signing
   paths, not the research harness.

**Revisit** only when a callgrind profile of a wall-clock-bound task
shows most instructions inside one field multiply.  Even then, prefer
`core::arch` intrinsics (`_mulx_u64`, `_addcarry_u64`, `pclmulqdq`)
behind runtime feature detection with the portable path kept as the
replay oracle, exactly as `Gf2` does today.  No inline `asm!`.

## 2. What the profile found (2026-10-09)

Profiled with `ecbench profile-input` piped into
`valgrind --tool=callgrind ecbench exec` at the last compiling commit
on `main` (`93ff0c72`; `main` itself did not build that day, see §4).

The ecbench hot paths are **already single-word arithmetic**, not
`BigUint`: `PrimeCurve` is `u64` affine, `BinaryGroup` is
`koblitz_fast::FastPoint` on `Gf2`.  The cost sat one level down:

| site | before | why it was slow |
|---|---|---|
| `ic_boundary::mulmod` | `(a as u128 * b as u128) % m as u128` | `__umodti3` library call on every multiplication |
| `ic_boundary::addmod` | `(a as u128 + b as u128) % m as u128` | `__umodti3` on every addition (additions are "free" in the unit, so nobody looked) |
| `ic_boundary::invmod` | extended Euclid on `i128` | `__divti3` library call per Euclid step; one inversion per affine group operation |
| `ecbench/generic.rs` | private copies of the first two | same, in BSGS and kangaroo coefficient updates |

Measured on the same `ecbench exec` child inputs (calibration spec
`research/ecbench_calibration_20261002`, target 0, round 0), callgrind
`Ir` between the solve markers, collection toggled inside the walk for
the per-function rows.  The recovered scalar and every counted phase are
identical before and after; only `*_ns` timings and placement differ.

| child | before (solve Ir) | after (solve Ir) | ratio |
|---|---|---|---|
| `rho.negation`, prime `p = 42652151` (26-bit, `r ≈ 2^25.3`) | 4,853,494 | 1,929,629 | 0.40 |
| `rho.signed_frobenius`, Koblitz `a = 0, n = 41` (`r ≈ 2^40`) | 131,784,217 | 129,040,220 | 0.98 |

Where the prime walk's instructions went, before → after:

| function | before | after |
|---|---|---|
| `PrimeCurve::add` (group law, inversion inlined) | 31.7 % | 33.2 % |
| `__divti3` + `u128_div_rem` + `__umodti3` + `__modti3` | 40.0 % | 0 |
| `rho_walk_with` (walk, table, canonicalisation) | 11.2 % | 29.8 % |
| `malloc_consolidate` + `unlink_chunk` | 11.7 % | 29.7 % |
| `PrimeCurve::double_raw` | 4.8 % | 5.6 % |

The Koblitz walk is bounded elsewhere: `Gf2::inv_clmul` 35.6 %,
`TunedWalker::walk` 27.6 %, `SignedFrobeniusClasses::canon_parts`
25.4 %, `FastCurve::add` 4.0 %, and `__umodti3` 3.5 % → 1.7 % (its
coefficient modulus `r ≈ 2^40` is above the 32-bit fast path; only the
additions sped up).  Wall time per child is single milliseconds at
these sizes and was not the measurement.

Replay: `ecbench verify --replay-all` on the committed `prime` and
`koblitz` calibration sessions passes with the patched binary (every
measured run re-executed with identical counts).

Tests: `cargo test --release --lib -- ic_boundary::tests::word_sized…
ic_boundary::tests::invmod_inverts… ecbench` at `93ff0c72` + this patch:
67 passed, 1 failed.  The failure,
`ecbench::record::tests::legacy_record_deserializer_reproduces_v1_decimal_semantics`,
is serde float parsing and needs serde_json's `arbitrary_precision`
feature, which `main`'s `Cargo.toml` enables and `93ff0c72`'s does not;
it touches nothing this change does.

## 3. What changed

- `mulmod`, `addmod`, `invmod` in `src/cryptanalysis/ic_boundary.rs`
  are word-sized with the 128-bit forms as the fallback for unreduced
  operands or a modulus above 32 bits.  Results are bit-identical (the
  inverse is unique; the fast paths are proved equal by the fallback
  condition), and `tests::word_sized_field_helpers_match_the_128_bit_forms`
  pins them to the old formulas on 200 000 random and edge inputs.
  `ecbench/generic.rs` imports them instead of keeping copies.
- Sealed ecbench records are unaffected: counted operations, keys,
  walks and recovered scalars are unchanged, only the instructions
  behind them.  `ecbench verify --replay-all` on an earlier session is
  the check.

### Release profile

`Cargo.toml` has no `[profile.release]` section, so the default (16
codegen units, no LTO) is what every one of the 56 workflows that run
`cargo build --release` ships.  Fat LTO with one codegen unit is being
measured against it on the same two child inputs (solve `Ir` and clean
build time); the decision and its numbers follow in this section.
`target-cpu=native` is not on the table, see §4.

## 4. TODO / feature list

- **TODO (method identity): lockstep walks with batched inversion for
  `rho.plain` / `rho.negation` / `kangaroo.vow`.**  Running `k` walks in
  lockstep and inverting their `k` denominators with Montgomery's trick
  turns `k` inversions into one plus `3(k−1)` multiplications, the
  single largest remaining wall-clock lever on both prime and binary
  curves.  It changes the order in which distinguished points reach the
  table, so which collision is found first changes and so do the counted
  steps.  Under `docs/ecbench/README.md` that is a **new method id**
  (`rho.plain_lockstep:k=…`), with its own bound and frontier entry,
  never a change to an existing method.  `koblitz_strong_rho` already
  does this for the signed-Frobenius walk and is the template.
- **TODO: allocation inside the prime rho walk.**  After the division
  fix, about 30 % of the prime walk's instructions are glibc
  `malloc_consolidate` / `unlink_chunk`: something in `rho_walk_with` or
  the distinguished-point table allocates per step or per walk.  Find
  it with `--toggle-collect` and hoist the buffer; this changes no
  count.
- **TODO: Montgomery or Barrett `mulmod` for 33–64-bit moduli.**  The
  fast path above stops at 32-bit moduli because that is every curve
  ecbench registers.  `roster_prime_instance` can carry wider primes; a
  per-modulus precomputed reciprocal (Möller–Granlund 2011) gives exact
  `u128 mod u64` in two multiplications.  Needs a place to keep the
  reciprocal: `PrimeCurve` is a serialised `{p, a, b}` and its hash feeds
  the curve id, so the reciprocal must live beside it, not in it.
- **TODO: `invmod` by binary GCD or Bernstein–Yang.**  Removes the
  remaining hardware division per Euclid step.  Only worth it once the
  lockstep item above has not removed most inversions.
- **TODO: aarch64 parity for `Gf2`'s folded reduction.**  `mul_fold`,
  `sqr_k_fold`, `inv_clmul` and `batch_inv_clmul` are x86_64 only;
  aarch64 takes the table reduction.  isolab's aarch64 hosts therefore
  measure a different instruction mix for the same counted work.
- **TODO: a perf-index kernel for the counted prime-curve rho.**
  `examples/perfbench` has no kernel through `ic_boundary::PrimeCurve`'s
  group law (`dlp/*` covers the Z/pZ, bsgs_fast and Koblitz paths), so
  `scripts/perf/perfindex.py` cannot see regressions in it.
- **TODO: constant-time inversion in the library's signing paths uses
  `BigUint` (`ecc::field::FieldElement`)**.  That is the one place where
  assembly's branch-freedom argument applies; the first step is the
  existing `U256` Montgomery type, not assembly.
- **Not planned: `target-cpu=native`.**  It breaks cross-host
  comparability of the same binary; keep runtime feature detection.
