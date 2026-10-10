# Post-quantum scheme speed

One table, one unit: **kilocycles per operation**, median of repeated calls, with
a correctness column. Produced by `cargo bench --bench pqc_speed`.

Every "fast" row is a second implementation of the same scheme, held to the
readable one byte for byte by differential tests. The readable implementations
stay; nothing here replaces them. `pqc::ml_kem` is validated against the NIST
ACVP vectors, so the differential tests carry that validation onto
`pqc::fast::ml_kem`.

## Machine

All numbers in this document come from one machine, and are only comparable to
each other:

| | |
|---|---|
| CPU | Intel Xeon @ 2.10 GHz |
| features | AVX2, AVX-512F, BMI2, AES-NI, SHA-NI |
| cores used | 1 |
| toolchain | rustc 1.98.1, `--release` |
| cycle source | `rdtsc`, which counts at a fixed reference frequency |

`rdtsc` counts reference cycles, not core cycles, so on a machine whose clock
boosts these numbers understate the core-cycle count. They are consistent within
a run, which is what comparing two implementations needs. Do not compare them
against published figures from other hardware without saying so.

## ML-KEM

| scheme | operation | before | after | speedup |
|---|---|---:|---:|---:|
| ML-KEM-512 | keygen | 215.7 | 41.8 | 5.16× |
| ML-KEM-512 | encaps | 244.8 | 42.6 | 5.75× |
| ML-KEM-512 | decaps | 221.1 | 56.3 | 3.93× |
| ML-KEM-768 | keygen | 341.5 | 70.4 | 4.85× |
| ML-KEM-768 | encaps | 391.0 | 68.7 | 5.69× |
| ML-KEM-768 | decaps | 353.6 | 88.3 | 4.00× |
| ML-KEM-1024 | keygen | 476.2 | 110.5 | 4.31× |
| ML-KEM-1024 | encaps | 537.5 | 103.7 | 5.18× |
| ML-KEM-1024 | decaps | 492.0 | 151.8 | 3.24× |

Decapsulation gains least because it re-encrypts to check the ciphertext, so it
pays for an encapsulation as well as a decryption and cannot fall below the
encapsulation row.

Correct in every row, where correct means the shared secret round-trips **and**
the reference implementation decapsulates the fast implementation's ciphertext
to the same secret. A fast implementation that is wrong in a self-consistent way
cannot show up here as a speedup.

### Where the time was

Before optimising anything, key generation was attributed between hashing and
ring arithmetic, because the received wisdom is that Keccak dominates the
lattice schemes and acting on that without checking would have wasted the
effort:

| ML-KEM-768 keygen | kilocycles | share |
|---|---:|---:|
| SHAKE128, nine matrix polynomials | 27.4 | 8% |
| SHAKE256, six CBD samples | 13.1 | 4% |
| everything else | 296.1 | 88% |
| total | 336.6 | |

So hashing was not the problem here, and a faster Keccak would have bought at
most a few per cent. The arithmetic was.

### What each change bought

The two largest items were not modular reductions but setup work repeated per
call:

- **Twiddle factors computed once, at compile time.** `ntt` and `ntt_inv` each
  opened by calling `zetas()`, which rebuilds all 128 twiddles by
  square-and-multiply: about 1800 modular multiplications to set up a transform
  that performs 1024. Worse, `multiply_ntts` called `zeta_pow` inside its loop,
  so one pointwise product ran 128 modular exponentiations. Both are pure
  functions of nothing, and are now `const` tables.
- **Byte packing through a word accumulator.** This was the single largest
  change, worth about 2.5× on its own: with the twiddle tables fixed but the
  packing still bitwise, ML-KEM-768 keygen sat at 208 kilocycles, and the
  accumulator took it to 70. `ByteEncode_12` had been running 3072 loop
  iterations, each with a division and a modulo by 8, to write 384 bytes, and
  encode or decode is called once per polynomial in every operation.
- **Signed coefficients and Montgomery multiplication**, replacing `u32`
  arithmetic with `% q`. This is the change one expects to matter most and it
  matters least of the three: `q` is a compile-time constant, so the compiler
  was already turning `% q` into a multiply and a shift. The gain is that
  intermediate values need reducing less often, not that each reduction is
  cheaper.
- **Streaming rejection sampling.** `sample_ntt` cannot know how many SHAKE128
  bytes it needs. The reference asked for a fixed length and, when it ran short,
  re-ran the sponge from the beginning with a doubled request, repeating every
  permutation already done. The fast one squeezes another block from the same
  state. This is rare enough to be worth little in cycles; it is here because
  the alternative is quadratic in the unlucky case.
- **No heap allocation in the hot paths.** The reference returned `Vec` from
  `byte_encode`, `prf` and every polynomial operator, so one key generation made
  hundreds of short-lived allocations.

### What was not done

- ~~No SIMD.~~ Done in the second round below, which supersedes this item;
  the first-round table above is kept as its "before".
- No constant-time work. This went the other way if anything: see below.

### Second round: public-data reuse and SIMD

The "after" column above became this round's baseline. Profiled again, the
picture had changed: in ML-KEM-768 encapsulation Keccak was ~45% of
instructions, expanding `Â` ~42%, and hashing `ek` for `H(ek)` ~10%. That is
the split Sayed, Taha and Nijjer report on a Cortex-M7 in *System-Level
Optimization Beyond Cryptographic Kernels* (arXiv:2610.01960), and their
first lever, caching work that depends only on public data, was the first
thing tried here. Evidence, per-step tables, host and commands are in
[`research/pqc_ml_kem_system_level_20261005/`](../research/pqc_ml_kem_system_level_20261005/README.md).

Class: **engineering**. Same algorithm, same bytes: every row is held to the
ACVP-validated reference, and keys, ciphertexts and shared secrets are
unchanged.

Kilocycles, median of 1400 interleaved samples per cell, isolated, on a
2.1 GHz Xeon with AVX-512 (a different session from the first round's table, so
compare within this table only). "prepared" is the same operation against a
key prepared once with `MlKemPreparedEncapsKey` / `MlKemPreparedDecapsKey`,
the preparation itself excluded and listed separately.

| scheme | operation | main | this round | prepared | main / round | main / prepared | correct |
|---|---|---:|---:|---:|---:|---:|:-:|
| ML-KEM-512 | keygen | 46.7 | 17.1 | 18.7 | 2.73× | 2.50× | yes |
| ML-KEM-512 | encaps | 54.2 | 15.9 | 8.3 | 3.40× | 6.55× | yes |
| ML-KEM-512 | decaps | 72.0 | 21.3 | 16.9 | 3.38× | 4.25× | yes |
| ML-KEM-768 | keygen | 83.1 | 28.0 | 29.6 | 2.97× | 2.81× | yes |
| ML-KEM-768 | encaps | 96.1 | 25.9 | 9.8 | 3.72× | 9.81× | yes |
| ML-KEM-768 | decaps | 119.8 | 33.1 | 22.3 | 3.62× | 5.38× | yes |
| ML-KEM-1024 | keygen | 136.3 | 36.5 | 37.8 | 3.74× | 3.61× | yes |
| ML-KEM-1024 | encaps | 149.0 | 34.5 | 13.0 | 4.32× | 11.43× | yes |
| ML-KEM-1024 | decaps | 180.5 | 44.4 | 30.0 | 4.07× | 6.01× | yes |

Preparing a key costs 7.7 / 16.2 / 21.6 kilocycles for the encapsulation
side and 5.6 / 12.5 / 16.1 for the decapsulation side (512 / 768 / 1024),
less than one unprepared encapsulation, so the cache pays for itself on its
first use. Key generation that captures the cache
(`ml_kem_keygen_prepared_internal`) costs within noise of plain key
generation, because it computes all of it anyway. The A/A spread on this
host is ±2%.

#### What each step bought, ML-KEM-768

Each row is a separate interleaved run against `main`, so the baseline
column moves by its noise; read the ratios.

| step | keygen | encaps | decaps | encaps, prepared | decaps, prepared |
|---|---:|---:|---:|---:|---:|
| main | 84.4 | 94.4 | 117.7 | | |
| public-data reuse | 82.8 | 93.6 | 117.6 | 43.4 | 70.7 |
| + four-way Keccak | 56.4 | 62.8 | 85.8 | 38.7 | 66.2 |
| + AVX2 NTT, inverse NTT, pointwise product | 38.4 | 38.4 | 52.9 | 14.3 | 32.3 |
| + grouped packing, single sponges in one AVX-512 lane | 30.4 | 30.3 | 38.5 | 11.8 | 24.1 |
| + vector compression | 30.8 | 28.1 | 35.2 | 9.8 | 21.8 |
| + vector rejection sampling | 28.0 | 25.9 | 33.1 | 9.8 | 22.3 |

- **Public-data reuse.** A prepared key holds `Â`, `t̂`, `H(ek)` and, for
  decapsulation, `ŝ`. On encapsulation that removes the matrix expansion
  and the `H(ek)` hash, 2.18× by itself; on decapsulation it removes the
  matrix from re-encryption, 1.66×: 54% and 40% fewer cycles. The paper
  measures 41–62% and 38–58% for its cached-`A` knob alone on the M7, so
  the lever transfers. What does not transfer is its ranking: there, the
  remaining system knobs were worth at most ~2.5%; here, once the cache
  existed, SIMD was worth more than the cache.
- **Four-way Keccak** (`pqc::fast::keccak4`). The k² matrix entries and the
  2k (+1) noise PRFs are independent sponges, so they run four to a
  256-bit register. On this host four AVX-512VL permutations (`vprolq`,
  `vpternlogq`) take ~610 cycles against ~880 for one scalar permutation,
  and four AVX2 permutations take ~1,500. So a leftover single sponge goes
  through the vector unit with AVX-512 and through the scalar code with
  AVX2.
- **AVX2 ring arithmetic** (`pqc::fast::ml_kem::avx2`). The largest single
  step. It is *bit-identical* to the scalar code, not merely congruent:
  `vpmulhw`/`vpmullw` give exactly the scalar Montgomery result, and the
  scalar Barrett quotient `⌊(va + 2²⁵)/2²⁶⌋` equals
  `⌊(⌊va/2¹⁶⌋ + 2⁹)/2¹⁰⌋` by the nested-floor identity. So coefficients
  keep the standard order (no pq-crystals-style permuted NTT and repacking)
  and the vector and scalar paths are tested for array equality.
- **Grouped packing.** Eight `d`-bit values are exactly `d` bytes, so
  `ByteEncode_d`/`ByteDecode_d` now assemble a group in one register. The
  byte-at-a-time accumulator it replaced had become 18% of a prepared
  encapsulation. The same step routed the scalar sponges (`H`, `G`, `J`,
  SHA3-512) through one lane of the AVX-512 kernel, ~630 cycles against
  ~880 per permutation; this also speeds up `pqc::fast::ml_dsa`, which
  shares the sponge.
- **Vector compression.** `Compress_d` divides by q. Done four lanes per
  `vpmuludq` with the reciprocal `⌈2³⁵/q⌉`, which is exact for every
  numerator below 13.7 million; ML-KEM's stay below 2²³. The test checks
  every input at every width.
- **Vector rejection sampling.** 24 bytes, 16 candidates per step, compacted
  with a 256-entry `pshufb` table. Worth ~2 kilocycles on the paths that
  still expand `Â`; nothing on prepared keys, which do not.

#### Per hardware class

The code picks AVX-512, AVX2 or scalar at run time, with no `target-cpu`
flag, so one binary runs everywhere. Each class measured on this host with
the others disabled, interleaved against `main` (ML-KEM-768 shown; all sets
in the research directory):

| class | keygen | encaps | decaps | encaps, prepared | decaps, prepared |
|---|---:|---:|---:|---:|---:|
| AVX-512VL + AVX2 | 2.97× | 3.72× | 3.62× | 9.81× | 5.38× |
| AVX2 only | 2.25× | 2.51× | 2.66× | 7.84× | 4.44× |
| scalar only | 1.16× | 1.17× | 1.21× | 2.26× | 1.77× |

No class got slower. The scalar row is the path Arm64 and pre-AVX2 x86
take, but it was measured on x86 with SIMD disabled, not on Arm64, which
has not been measured. There is no NEON backend.

Instructions per operation from callgrind, which is deterministic and does
not depend on the VM's neighbours (valgrind does not emulate AVX-512, so
these are the AVX2 path): ML-KEM-768 encapsulation 527,337 → 204,964, and
56,993 prepared; decapsulation 644,420 → 236,201, and 141,976 prepared.

#### Using the prepared keys

```rust
use crypto_lib::pqc::fast::ml_kem as kem;
use crypto_lib::pqc::ml_kem::ML_KEM_768;

// Key owner: generate and keep the cache, which costs nothing extra.
let (ek, dk, prepared) = kem::ml_kem_keygen_prepared(&ML_KEM_768);
// A peer holding only the standard encapsulation key prepares it once.
let peer = kem::MlKemPreparedEncapsKey::new(&ML_KEM_768, &ek).unwrap();
let (ct, ss) = peer.encaps();
assert_eq!(prepared.decaps(&ct), Some(ss));
```

The cache owns a copy of the key it came from, so it cannot be paired with
another key through this API; `ŝ` and `z` are zeroised on drop. A cache is
12 KiB for ML-KEM-768 (6 / 20 KiB for 512 / 1024).

#### Not done, and not claimed

- ~~**No comparison with the pq-crystals AVX2 implementation.**~~ Since
  made on this host, see the next section.
- **`J(z ‖ c)` is still computed on every decapsulation.** Skipping it when
  the ciphertext checks out would be the cheapest remaining saving (nine
  permutations at ML-KEM-768) and would hand an attacker a timing oracle on
  ciphertext validity, the plaintext-checking oracle that chosen-ciphertext
  attacks on Kyber are built from. It stays.
- **No AVX-512 ring arithmetic, no NEON.** The NTT is AVX2 only.

### The reference: pq-crystals AVX2 on the same host

Until this row existed the only comparison was against our own `main`.
pq-crystals' AVX2 Kyber (FIPS 203 version, `3edd5af`), built with its own
`Makefile` (`-march=native`) on the same machine and run interleaved with
ours, is the boundary. Kilocycles; ratio is ours ÷ reference, so above 1 is
slower. Details, kernel-level breakdown and raw runs:
[`research/pqc_ml_kem_reference_20261005/`](../research/pqc_ml_kem_reference_20261005/README.md).

| operation | pq-crystals AVX2 | ours | ratio | ours, prepared | ratio |
|---|---:|---:|---:|---:|---:|
| ML-KEM-768 keygen | 18.10 | 28.1 | 1.55 | 29.5 | 1.63 |
| ML-KEM-768 encaps | 18.48 | 26.0 | 1.41 | 9.7 | 0.52 |
| ML-KEM-768 decaps | 21.20 | 33.2 | 1.57 | 21.8 | 1.03 |

So after this round we are still **1.35–1.63× behind** the reference
across all three parameter sets on the operations both implement. We are
ahead only with a prepared key, an API the reference does not have and
would gain from too. Per kernel, our NTT, inverse NTT and pointwise
product are 1.8–1.95× slower than its hand-scheduled assembly, and our
matrix expansion 1.25–1.31× slower; those explain about half of the
encapsulation gap, and the rest is not yet attributed. The next round's
success condition, stated in that note before any change: every unprepared
ML-KEM-768 operation within 1.10× of this reference.

#### Where the gap is (profile)

[`research/pqc_ml_kem_profile_20261005/`](../research/pqc_ml_kem_profile_20261005/README.md)
attributes it, by callgrind on the same AVX2 path for both and by per-phase
cycles. Two corrections to the table above come with it:

- **Accounting.** That table's ratios are one sample. Five interleaved
  rounds put ML-KEM-768 encapsulation at 1.37–1.63× the `-march=native`
  build and 1.0–1.3× an AVX2-only build of the same reference, and two
  builds of our own harness differed by ~12% on identical code. The stable
  figure is instructions: **1.40×** encaps, 1.45× decaps, 1.50× keygen.
- **Attribution.** The gap is mostly bookkeeping, not arithmetic. Our
  pack/unpack/compress/decompress and one-bit message kernels are 5–11×
  pq-crystals' per call; ring arithmetic 1.7–2.0×; zero-filling and copying
  polynomials ~23k instructions per encapsulation; byte-wise constant-time
  compare ~13k per decapsulation. Our matrix expansion is already faster
  than the reference's AVX2 build here, by ~1k cycles, thanks to the
  AVX-512 Keccak.

## ML-DSA

ML-DSA-65, median of four back-to-back runs. Absolute counts on this machine
wandered by up to 30% between runs because other builds were competing for
cores; the ratios are the stable quantity, because both rows of a pair are
measured back to back in the same process.

| operation | reference | fast | speedup | per-run range |
|---|---:|---:|---:|---|
| keygen | 368.0 | 289.2 | 1.27× | 1.24–1.29 |
| sign | 1220.7 | 565.4 | 2.16× | 1.99–2.31 |
| verify | 415.9 | 253.4 | 1.64× | 1.57–1.69 |

The `correct` column asserts byte equality of the public key, secret key and
signature against the reference, and that each implementation verifies the
other's signature. A self-consistent wrong answer cannot pass it.

An independent run on a quieter machine gave 292.4 → 208.6, 985.3 → 450.2 and
314.2 → 200.3, that is 1.40× / 2.19× / 1.57×. Signing and verification agree
with the table; key generation came out above the four-run range quoted there,
so treat the keygen ratio as roughly 1.3× with real spread rather than a
figure good to two digits. Every absolute number on this page was taken on a
machine with other builds running, and the spread is wider than any single run
suggests.

### What each change bought

Four configurations built in one process out of tree, so they are directly
comparable. A and B are the fast module with exactly one piece swapped back to
the reference's; all four produce byte-identical keys and signatures.

| configuration | keygen | sign | verify |
|---|---:|---:|---:|
| reference | 390.8 | 1254.5 | 419.4 |
| A: streaming sponge, reference `% q` arithmetic | 331.6 | 1084.4 | 375.9 |
| B: Montgomery arithmetic, reference one-shot XOF | 301.7 | 631.7 | 271.1 |
| fast (both) | 260.4 | 598.4 | 265.4 |

Isolating one knob at a time:

- **Montgomery arithmetic** (A to fast): −21% keygen, −45% sign, −29% verify.
  The large item by far.
- **Streaming sponge** (B to fast): −14% keygen, −5% sign, −2% verify. Real but
  modest, and it shrinks as the arithmetic gets faster.
- **Residual**, the fixed-array restructuring mixed with the non-additivity of
  the other two: −18 / −137 / −38 kilocycles. Not separated further, so it is
  reported as one figure rather than split on a guess.

### Why this is 2× and not 5×

The ML-DSA reference does **not** have the two pathologies the ML-KEM reference
had. Its zeta table is already a compile-time constant, its packing is
nibble-oriented rather than bit-at-a-time, and it already hoists the matrix
expansion and the hashed prefixes out of the rejection loop. No factor of ten
was available.

What remained was arithmetic, not hashing, which is the opposite of the usual
expectation for a lattice signature. The reference multiplies with
`(a as i64 * b as i64) % q` and reduces every addition with `rem_euclid`, so a
64-bit division-by-constant sequence runs on roughly 30,000 butterflies and
12,000 pointwise products per signing attempt. Replacing that is 45% of
signing; the XOF is 5%.

### A convention trap between the two modules

The reference's `inv_ntt` multiplies by n⁻¹ and is a true inverse, so
`inv_ntt(ntt(p)) == p`. The fast one is the pq-crystals `invntt_tomont`, which
also scales by R, cancelled by the R⁻¹ that the pointwise multiply introduces.
So in the fast module `inv_ntt(ntt(p)) == R·p mod q`, and only the composition
with a pointwise multiply is the identity. Anyone moving a helper between the
two files needs to know this.

### On comparing against published figures

pq-crystals publishes, for round-3 Dilithium3 on one core of a Core-i7 6600U:
C reference 544,232 / 2,348,703 / 522,267 and AVX2 256,403 / 529,106 / 179,424.
Those are **not** recorded as a row above, deliberately. They are round-3
Dilithium rather than final FIPS 204 ML-DSA-65, on different hardware, and the
C reference has not been built on this machine. A peer row would need a
same-machine build.

## Isogeny arithmetic (SQIsign verification core)

### What this is, and what it is not

The existing `pqc::sqisign` is a toy: p = 431, 37 supersingular j-invariants, a
precomputed graph, breadth-first search standing in for the KLPT algorithm, and
a signing function that ignores the secret key. Making *that* faster is
meaningless.

`pqc::fast::isogeny` is instead the layer where a real SQIsign spends its time:
field and isogeny arithmetic at real parameter sizes. **It is not SQIsign.** A
real implementation still needs quaternion order arithmetic, the Deuring
correspondence and KLPT for signing, which is tens of thousands of lines and
was deliberately not attempted. What is here is roughly the verification-side
arithmetic, which is self-contained and has a clean correctness specification.

### The prime

`p = 3·2^324 − 1`, 326 bits, six 64-bit limbs. This is the **genuine SQIsign
NIST level-1 prime**, not a stand-in. It was established rather than recalled:
the official reference implementation was cloned at a named commit and the
value read out of its parameter file, then cross-checked two ways before any
constant was written down, with the limb decomposition and the recomputed
`R mod p` and `R² mod p` matching the reference's own constants limb for limb.
It was also Miller-Rabin tested, and `p ≡ 3 mod 4` confirmed so that
`F_p[i]/(i²+1)` is a field.

Worth recording because I got this wrong when scoping the work: I said level-1
was about 254 bits. That was the 2023 round-1 submission. The current
dimension-2 SQIsign moved to 326, which is six limbs rather than four, so every
field multiplication costs roughly 2.25× what a 254-bit one would.

The shape pays for itself twice. Because `p ≡ −1 mod 2^64` the Montgomery
constant is 1, so the quotient digit needs no multiplication; and because
`m·p = 48m·2^320 − m`, the `−m` cancels the low limb exactly with no borrow, so
dividing by `2^64` is a pure limb shift. A reduction round costs six
multiplications and one small one instead of twelve.

### Measured

Cycles per operation, from the final tree:

| operation | cycles |
|---|---:|
| F_p multiply (6 limbs) | 61 |
| F_p inverse (addition chain, 324S + 11M) | 16,918 |
| F_p² multiply (Karatsuba, 3M) | 212 |
| F_p² square (2M) | 144 |
| F_p² inverse | 17,700 |
| xDBL (A24plus : C24) | 1,270 |
| xDBL (a24 normalised) | 1,025 |
| xADD | 1,420 |
| ladder [k]P, k ≈ 2^326 | 768,648 |
| 2-isogeny step (codomain + evaluation) | 1,335 |
| 2-isogeny evaluation only | 1,028 |
| 4-isogeny step (codomain + evaluation) | 2,479 |
| **2^324 chain, optimal strategy** | **3,055,572** |
| 2^324 chain, naive walker | 34,685,256 |

These are **dependency-chain latencies**, measured with each operation feeding
the next 128 deep, which is the pessimistic end. Throughput with instruction
level parallelism would be better. An independent run of an earlier build gave
lower absolute figures throughout, for the reason in the next paragraph.

A 326-bit field multiplication at 61 cycles is in the range a tuned
implementation should reach. Scaled to 254 bits, four limbs and so roughly
`(4/6)²` of the multiplications, that is about 27 cycles-equivalent, the "tens
of cycles" figure one wants, in plain Rust with no `unsafe`, no `mulx`/`adcx`
intrinsics and no assembly. The extension multiply at 212 is 3.5× the base
one, as Karatsuba's three products plus five modular additions should cost.
Hand-written ADX assembly would likely be 1.5× to 2× faster; that is an
estimate, not a measurement, since no such implementation was benchmarked
here.

**Read the sub-100-cycle rows with care.** At this size the numbers move
*between builds* far more than within them: two back-to-back runs of one binary
agree to 0.5%, but adding or removing unrelated code shifts a field-multiply
figure by a quarter, and the multiply has been seen anywhere from 48 to 61
cycles across builds of the same algorithm. Code layout and inlining decisions
dominate. This was established the hard way, by removing a candidate
optimisation and finding that the two benchmark rows, then running identical
code, still differed by 5%. Any comparison of two primitives at this scale that
rests on a difference under about 10% is measuring the linker, not the
algorithm. The chain and ladder rows, which are thousands of times larger, do
not have this problem.

**The optimal strategy is worth 11.4× over the naive walker** on a full-length
chain, 3.06 against 34.7 megacycles. That is the `Θ(n log n)` against `Θ(n²)`
gap at n = 162. That is the single largest structural win
in this module and the reason the naive walker is kept: it is the thing the
strategy is measured against, and the two are checked to produce the same
isogeny.

### A negative result worth keeping

A dedicated squaring using the usual off-diagonal trick, 21 limb
multiplications instead of 36, was written, differentially tested, and measured
at 51 cycles against the general multiply's 48. No faster, and within this
shared machine's run-to-run noise, so not reliably slower either. With the CIOS
reduction already integrated and cheap, the 12-limb doubling shift and the
separated reduction's carry loop appear to cost about what the 15 saved
multiplications save. Forty lines of carry plumbing that buys nothing
measurable was dropped, so `Fp::sqr` is `mul(a, a)`. This is a
not-measured-any-better result, not a proof that the trick cannot be made to
win with a more integrated reduction.

### Correctness

Seventeen tests. The load-bearing ones are the differential tests against
`num-bigint`, which reimplement the same operations schoolbook-style and
compare over thousands of seeded random inputs, for both the base field and the
quadratic extension. On top of those: field axioms, Legendre symbol and square
roots, the ladder against repeated addition and its multiplicativity, points
staying on the curve under doubling, a 4-isogeny equalling two 2-isogenies, the
optimal-strategy chain agreeing with the naive walker, a full-length chain of
degree 2^324, the supersingular group order, and the modular polynomial
recognising known pairs.

Two checks are independent of the implementation's own formulas, which matters
because a differential test against your own algebra can agree on a shared
misconception:

- **The classical modular polynomial** `Φ₂(j(E), j(E')) = 0` on every 2-isogeny
  codomain, which knows nothing about Montgomery models or projective curve
  coefficients. Its coefficients were themselves checked against four classical
  pairs first, so a typo cannot make the test vacuous.
- **The 4-isogeny against two composed 2-isogenies**, compared on j-invariants,
  which is convention-free. This one earns its keep: the reference
  implementation's 4-isogeny is the Costello–Longa–Naehrig formula composed
  with `x ↦ −x`, the opposite sign convention from the 2-isogeny taken from
  Renes, and the test confirms the two are the same isogeny up to that
  isomorphism.

Four tests failed on the first pass and all four were worth having. One was a
bug in the *test*, which computed a modular inverse by Fermat's little theorem
on a composite modulus. Two surfaced a real precondition: the random kernel
sometimes hit `(0,0)`, the one 2-torsion point whose Montgomery 2-isogeny needs
a square root to name its codomain, which every x-only implementation excludes;
that exclusion is now documented and the generator rejects it. The fourth found
an API footgun in the new code itself: `PartialEq` had been derived on the
projective types, which is simply wrong for projective coordinates, and the two
chain walkers reach the same curve by different representatives. The derive is
gone, replaced by an explicit same-curve test.

### What a real SQIsign still needs

So that nobody reads this as more than it is. On top of what is here: the
quaternion algebra and its maximal orders and ideals with lattice reduction;
KLPT equivalent-ideal search; the Deuring correspondence's ideal-to-isogeny
translation; the theta-model (2,2)-isogeny machinery that dimension-2 SQIsign
uses; and deterministic point compression, torsion-basis generation and the
Fiat–Shamir plumbing. The first three are the signer, which is essentially
absent here. That list is tens of thousands of lines.

There is also **no known-answer test**, because no test vector exists for this
layer in isolation and the reference C was not linked against to generate one.
Correctness is internal consistency plus the two independent checks above.
Cross-checking against the reference's own field and curve test binaries is the
obvious next strengthening.

## Security

**These implementations are not constant-time and must not be used for
anything real.** That is true of the whole library, and more so here: speed and
side-channel resistance pull in opposite directions, and only the first was the
goal. The rejection samplers branch on their output and the compressions are
written for throughput. See [`SECURITY.md`](../SECURITY.md).
