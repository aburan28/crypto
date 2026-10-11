# Side-channel assessment of the Eigen-Séta implementation

Scope: this crate (`eigen_seta/`) at the commit that adds this file. This is a
research implementation. Nothing here has been reviewed by a third party.

## 1. Threat model

* **In scope.** Timing, including remote timing, and microarchitectural timing
  from secret-dependent branches or memory addresses. That covers the cache,
  the branch predictor and the instruction cache.
* **Out of scope.** Power and EM analysis, fault injection, speculative
  execution beyond the absence of secret-indexed loads, and the operating
  system (swap, core dumps).
* **Platform assumption.** 64×64→128-bit multiplication has data-independent
  latency. This holds on current x86-64 and on AArch64 application cores. It
  does not hold on some embedded cores, for example Cortex-M3.

## 2. What is secret

| Object | Lifetime | Where |
|---|---|---|
| `lam_s`, `lam_t`, `cp2` (mod 2^E), `c3`, `cp3` (mod 3^F3) | long-term | `SecretKey` |
| `z` (implicit-rejection seed) | long-term | `SecretKey` |
| secret degree f, α, the Kani kernel | key generation only | not stored |
| s (the KEM seed), masks (a, b), φ_s, E_C basis images | per encapsulation | stack and heap temporaries |
| everything computed from the secret key inside decapsulation | per decapsulation | stack and heap temporaries |

Public data: the public key, the ciphertext, and all parameters. That
includes the bit length of f, which is fixed by the range [f_lo, 2 f_lo).

## 3. Structure: public phase vs. secret phase

Decapsulation is split so that each phase can be checked on its own.

1. `decaps_prepare(sk, ct)` is **variable time**. It depends only on the public
   key and the ciphertext. It parses both, lifts S and T from x-coordinates,
   fixes their relative sign with x(S − T), and computes the deterministic
   canonical bases of E_C[3^F3] and E_A[3^F3] (counter-driven point search).
   A malformed ciphertext returns the implicit-rejection key here. That
   reveals only that the public input was malformed.
2. `decaps_core(sk, prep, ct)` is **intended constant time**. It covers the
   trapdoor inversion, the deterministic re-encryption, the comparison and
   the constant-time selection between the real key and the rejection key.
   It is `#[inline(never)]` so that tools can measure it alone.

Encapsulation runs `enc_det_core`. That is constant time in s, after the
same public canonical-basis computation on E_A.

## 4. Per-component audit of the secret phase

"CT" here means no branch and no memory index depends on secret data, and
loop counts are public. Each row was checked by reading the code. Section 5
gives the tool-based check.

| Component | File | Status | Notes |
|---|---|---|---|
| F_p mul (CIOS; Montgomery-friendly SOS when p + 1 = 0 mod 2^64) | `fp.rs` | CT | Fixed loops. Reduction uses a borrow mask. The friendly path skips the zero limbs of p + 1, which is a public constant. |
| F_p add, sub, neg | `fp.rs` | CT | Borrow-masked final subtraction. |
| F_p inversion: Bernstein–Yang safegcd | `safegcd.rs` | CT, after a fix | Public divstep count ⌊(49B + 57)/17⌋, rounded up to a multiple of 62. Swaps and negations are masks. 0⁻¹ = 0. Checked against Fermat at 1, 2, 9, 21 and 42 limbs (`tests/field.rs`). The sign masks needed `black_box` (§5). |
| F_p sqrt, Legendre | `fp.rs` | CT | Fixed public exponents. |
| F_p² mul, square, inversion via the norm | `fp2.rs` | CT | |
| F_p² sqrt (Adj–Rodríguez-Henríquez Alg. 9) | `fp2.rs` | CT | Both branches computed, result chosen by `conditional_select`. |
| Montgomery ladder, 3-point ladder | `mont.rs` | CT | Fixed-width scalars. The `nbits` iteration count is public. `cswap` uses masks. |
| y-recovery (Okeya–Sakurai), `smul_ct` | `mont.rs`, `kem.rs` | CT | See residual issue R3. |
| 3-isogeny chain (x-only, fixed strategy) | `isog3.rs` | CT | The strategy depends only on F3. |
| Weierstrass add, double, dbl_iter | `weier.rs` | CT | Doubling and addition both computed, then selected. Infinity flags are `Choice`. |
| Vélu 2-isogeny, isomorphism `WIso::find` | `weier.rs` | CT | `WIso` picks a square root in CT. The sign is resolved later by CT selection over both candidates. |
| Product-kernel step | `chain.rs` | CT | The 2-torsion generator is chosen by `conditional_select` on `inf`. |
| Gluing, Richelot, splitting (Castryck–Decru formulas) | `kani.rs` | CT | Fixed-degree polynomial arithmetic. The inverse modulo the quartic uses a 4×4 adjugate, not a Euclidean algorithm. Batch inversion uses fixed products. |
| Decapsulation chain schedule | `chain.rs` | CT | Graph, product kernel, glue, Richelot×(E−6), split, product kernel, graph. The schedule is public. Key generation rejects keys for which it does not hold. |
| Factor identification (swap, ±) | `kem.rs` | CT | All 2×2×2 candidates evaluated, the right one kept by `conditional_select`. |
| Weil pairing (Miller loop over public bits of 3^F3) | `pairing.rs` | CT | Projective, exact line scaling. One inversion at the end. |
| DLP in μ_{3^k} (Pohlig–Hellman) | `pairing.rs` | CT | k fixed rounds. Digits found by `ct_eq` and masks. Inversion via Frobenius. |
| `OddModulus` reduce, mul (mod 3^F3) | `scalar.rs` | CT | Fixed widths, masked final subtraction. |
| Scalars mod 2^E | `scalar.rs` | CT | Fixed width. `mod_2k` branches only on public k. Inversion by Newton with a public iteration count. |
| SHAKE256 (G, K, R, P) | `sha3` crate | CT | Inputs have public lengths. |
| FO comparison and key selection | `ct.rs`, `kem.rs` | CT | `subtle::ConstantTimeEq` over the whole ciphertext, then byte-wise `conditional_select`. |

`choice_from_bit` and the masks pass through `core::hint::black_box` so that
LLVM cannot turn the selections into branches. Rust gives no guarantee of
this, which is why section 5 checks the compiled binary.

## 5. Tool-based check: instruction counts of the secret phase

`examples/ct_probe.rs` builds keys and ciphertexts from fixed seeds with a
test RNG. It runs `decaps_prepare` outside the measured region, then runs
`decaps_core` inside `probe_core`. Callgrind counts the instructions executed
in `probe_core` only:

```
CARGO_PROFILE_RELEASE_DEBUG=true cargo build --release --example ct_probe
valgrind --tool=callgrind --toggle-collect='*probe_core*' \
    ./target/release/examples/ct_probe <key_seed> <ct_seed> <valid|foreign>
```

`foreign` passes a well-formed ciphertext made for a different public key.
It goes through the whole trapdoor and is rejected by the FO
re-encryption check.

**Results.** These come from run `run-20261011T1955Z-verify`, built from
commit d504af28e with a clean worktree, using rustc 1.97.0 in release mode
with debuginfo, on x86-64 with glibc malloc. Records are in
`research/eigen_seta_rust_20261011/reports/run-20261011T1955Z-verify/records/`.

* **Toy16.** 4 key seeds, each with two valid ciphertexts and one foreign
  ciphertext.
* **L1, the NIST level-1 parameters.** Two keys, two valid ciphertexts under
  key 1, and one foreign ciphertext.

"Allocator" means glibc's malloc and arena code, heap growth (`sbrk`, `brk`,
`morecore`) and the `mem*` routines.

| set | key seed | ciphertext | accepted | total instructions | allocator | non-allocator |
|---|---|---|---|---|---|---|
| Toy16 | 1 | foreign | no | 1,148,059,231 | 31,850,015 | **1,116,209,216** |
| Toy16 | 1 | valid 1 | yes | 1,147,950,979 | 31,741,763 | **1,116,209,216** |
| Toy16 | 1 | valid 2 | yes | 1,147,950,979 | 31,741,763 | **1,116,209,216** |
| Toy16 | 2 | foreign | no | 1,147,693,051 | 31,483,835 | **1,116,209,216** |
| Toy16 | 2 | valid 1 | yes | 1,147,388,400 | 31,179,184 | **1,116,209,216** |
| Toy16 | 2 | valid 2 | yes | 1,147,388,400 | 31,179,184 | **1,116,209,216** |
| Toy16 | 3 | foreign | no | 1,147,849,768 | 31,640,552 | **1,116,209,216** |
| Toy16 | 3 | valid 1 | yes | 1,147,459,087 | 31,249,871 | **1,116,209,216** |
| Toy16 | 3 | valid 2 | yes | 1,147,459,087 | 31,249,871 | **1,116,209,216** |
| Toy16 | 4 | foreign | no | 1,147,787,126 | 31,577,910 | **1,116,209,216** |
| Toy16 | 4 | valid 1 | yes | 1,147,477,368 | 31,268,152 | **1,116,209,216** |
| Toy16 | 4 | valid 2 | yes | 1,147,477,368 | 31,268,152 | **1,116,209,216** |
| L1 | 1 | foreign | no | 124,257,288,150 | 3,260,858,540 | **120,996,429,610** |
| L1 | 1 | valid 1 | yes | 124,250,795,358 | 3,254,365,748 | **120,996,429,610** |
| L1 | 1 | valid 2 | yes | 124,250,795,358 | 3,254,365,748 | **120,996,429,610** |
| L1 | 2 | valid 1 | yes | 124,262,712,389 | 3,266,282,779 | **120,996,429,610** |

The non-allocator count is identical in all 12 Toy16 runs and in all 4 L1
runs. The allocator count varies with the key and with the mode, because the
heap state on entry depends on the earlier, variable-time key generation and
public phase. For a given key it is the same for both valid ciphertexts.
See R6.

At L1 the heap grows during the secret phase. A first summary that did not
count `sbrk` and `brk` as allocator code showed a 174-instruction
difference between the two keys. All of it was in those two libc
heap-growth routines, so it is allocator work, not code of this crate.

An earlier run, `run-20261011T1827Z-sidechannel` on commit 19c641d12, gave
the same picture at Toy16, with 1,116,186,003 non-allocator instructions in
all 12 runs. The count is different here because the hardening commit adds
buffer wiping to the secret phase.

**A leak this check found and fixed.** The first version of the probe showed
counts that differed between keys and between ciphertexts in
`safegcd::modinv`, by 180 instructions in one comparison. Line-level
annotation traced it to `normalize`. There, LLVM had compiled the masked
add `r.v[i] + (m.v[i] & cond_add)` into a conditional that skipped the add
when the mask was zero. That leaks the sign of a secret-dependent
intermediate on every field inversion. The masks now pass through
`core::hint::black_box`, and the counts above were taken after the fix.
The source had no branch at all: only a check of the compiled code can
catch this class of problem.

**Statistical timing test.** `examples/dudect_decaps.rs` times the secret
phase on one Toy16 key for two classes, interleaved in random order. Class 0
is one fixed ciphertext. Class 1 is a pool of 16 fresh ciphertexts, half of
them FO-rejected. It reports Welch's t, the dudect method. The result
with 1,000 samples per class, in the same verification run, after the
probes and alone on the machine:

```
samples/class 1000; median ns: fixed 122373850, random 120844155
t(raw) = +1.17
t(crop <= p50) = +0.14  (n = 473, 528)
t(crop <= p75) = +1.46  (n = 730, 771)
t(crop <= p90) = +2.36  (n = 899, 902)
```

|t| stays under the usual 4.5 threshold; the largest value is 2.36. This is weaker evidence than
the instruction counts: it is one key, wall-clock time on a shared VM, and
it compares ciphertext classes rather than keys. It does add coverage of
effects that instruction counts miss, such as variable-latency
instructions.

What this check does and does not show:

* Equal instruction counts across secret keys and across accepted and
  rejected ciphertexts are strong evidence of secret-independent control
  flow.
* They do not show secret-independent memory addresses. The code has no
  secret-indexed tables. All selections are full-width masks over values
  that are all loaded.
* They do not cover other compilers, other targets, or other optimisation
  settings. Rerun the probe on the deployment build.

## 6. Residual issues and known leaks

* **R1. Key generation is variable time.** Rejection sampling on f, Cornacchia,
  Miller–Rabin, the torsion-basis search and `Omega::find` use `num-bigint`
  and early exits. The number of rejected candidates and the Cornacchia
  running time depend on f. Run key generation in a trusted environment with
  no co-resident attacker. A constant-time key generation, especially a
  constant-time norm-equation solver, is open work.
  Running the probe exposed a liveness bug in key generation, now fixed.
  Tonelli–Shanks in `ntheory::sqrt_mod_prime` accepted z = 0 as a
  "non-residue", after which it looped forever. This happened for some RNG
  streams on small primes q ≡ 1 mod 8 inside `two_squares`. The regression
  test is `tests/kani_unit.rs::sqrt_mod_prime_small_q_terminates`.
* **R2. Encapsulation and decapsulation have public variable-time parts.**
  The canonical-basis search and the point lifting depend on the curves A and
  C, which are public. This leaks nothing about s or the secret key, but an
  implementation must not mix secret data into them.
* **R3. Degenerate scalar edge case.** `smul_ct` with y-recovery returns the
  point at infinity when [k + 1]P = O, that is k ≡ −1 modulo the order of P.
  Key generation now rejects keys whose c′, λ_S or λ_T is ≡ −1 mod 2^E, or
  whose c or c′ is ≡ −1 mod 3^F3, so the secret-key scalars never hit it.
  One case remains: the ciphertext-dependent DLP value k used to recover
  ker φ̂, which hits it with probability about 3^−F3, below 2^−256 at L1.
  The code has no branch there: the result is wrong, the timing is the same,
  and the FO check rejects.
* **R4. Degenerate inputs to the (2,2)-chain.** If a formula hits a zero
  denominator, for example an intermediate divisor that is not in general
  position, the output is garbage, with no branch or early exit. This gives
  an attacker no oracle. Any ciphertext that is not an honest encryption is
  answered with H_R(z, ct), whether or not a degenerate input occurred. Only
  the failure of an *honest* ciphertext would be visible, as a mismatch with
  the encapsulated key. Its probability is heuristically O(E/p), and none
  occurred in the recorded runs: 3,000 trials at Toy16, plus Toy32, L1, L3
  and L5 trials, and the robustness test.
* **R5. Temporaries are only partly wiped.** `SecretKey` is
  `Zeroize + ZeroizeOnDrop`, and its `Debug` impl redacts the secrets. Secret
  scalars (`Scalar` is `ZeroizeOnDrop`) are wiped on drop, and so are the
  secret byte buffers: the encoding of s, the mask hash, the shared secret
  and the rejection key. Points and field elements of decapsulation are not
  wiped.
* **R6. Heap allocation in the secret phase.** `Vec`s are allocated with
  sizes fixed by the parameter set, so the allocator calls do not depend on
  secrets. glibc malloc's own instruction count does depend on heap state
  from earlier public work. Section 5 therefore reports allocator functions
  separately. A production implementation should use fixed stack buffers.
* **R7. Compiler and hardware.** See section 5 and the platform assumption in
  section 1.
* **R8. Power, EM and faults** are not addressed. The FO re-encryption check
  catches faults that corrupt the recovered s, but not faults that skip it.
