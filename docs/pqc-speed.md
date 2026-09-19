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

- No SIMD. The machine has AVX2 and AVX-512F and the NTT is the obvious
  candidate; a vectorised butterfly would be worth considerably more than
  anything above. It is not written, so it is not claimed.
- No constant-time work. This went the other way if anything: see below.

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

The prime is `p = 3·2^324 − 1`, six 64-bit limbs, chosen for the right shape
(`p + 1` divisible by a large power of two) rather than copied from the NIST
submission. The module documents that it is representative rather than the
submitted parameter.

### Measured

| operation | kilocycles | per second |
|---|---:|---:|
| F_p multiply (6 limbs) | 0.040 | 52.9 M |
| F_p square (= multiply, see below) | 0.045 | 46.4 M |
| F_p inverse (addition chain) | 13.675 | 154 k |
| F_p² multiply (Karatsuba) | 0.167 | 12.6 M |
| F_p² square | 0.113 | 18.6 M |
| F_p² inverse | 13.902 | 151 k |
| xDBL (A24plus : C24) | 1.002 | 2.10 M |
| xDBL (a24 normalised) | 0.910 | 2.31 M |
| xADD | 1.262 | 1.66 M |
| ladder [k]P, k ≈ 2^326 | 923.9 | 2.3 k |
| 2-isogeny step (codomain + evaluation) | 1.102 | 1.91 M |
| 2-isogeny evaluation only | 0.846 | 2.48 M |
| 4-isogeny step (codomain + evaluation) | 2.208 | 951 k |
| 2^324 chain, optimal strategy | 2524.5 | 832 |
| 2^324 chain, naive walker | 28566.1 | 74 |

A 326-bit field multiplication at 40 cycles is in the range a tuned
implementation should reach, and the quadratic extension multiply at 167 is
about 4.2× that, consistent with Karatsuba's three base multiplications plus
reduction and carry overhead.

**The optimal strategy is worth 11.3× over the naive walker** on a full-length
chain, 2.52 against 28.57 megacycles. That is the single largest structural win
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

Four of these failed on the first pass, in the two- and four-isogeny formulas,
the strategy-versus-naive comparison and the prime constants. That is what
those tests are for.

## Security

**These implementations are not constant-time and must not be used for
anything real.** That is true of the whole library, and more so here: speed and
side-channel resistance pull in opposite directions, and only the first was the
goal. The rejection samplers branch on their output and the compressions are
written for throughput. See [`SECURITY.md`](../SECURITY.md).
