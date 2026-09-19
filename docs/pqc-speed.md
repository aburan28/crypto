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

## Security

**These implementations are not constant-time and must not be used for
anything real.** That is true of the whole library, and more so here: speed and
side-channel resistance pull in opposite directions, and only the first was the
goal. The rejection samplers branch on their output and the compressions are
written for throughput. See [`SECURITY.md`](../SECURITY.md).
