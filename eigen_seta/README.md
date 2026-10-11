# eigen_seta: Rust implementation of the Eigen-Séta KEM

This is a research implementation of **Eigen-Séta v2**, the secret-degree
variant. It is an isogeny-based KEM whose trapdoor is a secret endomorphism α
of a public curve. Decryption is a "lollipop" endomorphism decoded by a
(2^E, 2^E)-isogeny chain in dimension 2, built with Kani's lemma. The design
and the earlier Sage toy prototype are in `RESEARCH_ENDRING_PKE_STRATEGY.md`
and `research/eigen_seta_prototype_20261011/`. Section 4 corrects two
statements made there.

**Status.** Research code. It decapsulates correctly at the NIST level 1, 3
and 5 parameter sizes. Decapsulation is written to be constant time, and the
instruction-count check is in `SIDE_CHANNELS.md`. The security levels are
**conjectured** targets from the analysis below. No cryptanalysis by others
exists. Do not use this to protect data.

## 1. The scheme in one paragraph

Let E_D be the curve with CM by Z[(1+√−7)/2] (j = −3375) over F_p, with p =
2^(E+1)·3^F3·c − 1. Key generation draws a secret odd degree f ∈ [f_lo,
2 f_lo) and solves N(x) + p·N(y) = f(2^E − f) by Cornacchia for α ∈
End(E_D) of degree f(2^E − f). It then evaluates the Kani kernel
{([f]R, αR)} with a (2,2)-chain. The degree-f factor gives the public curve
E_A and the 2^E-torsion images. The public key is the curve plus x(λ1·P),
x(λ2·Q) and x(P − Q), where the eigenbasis scaling is part of the secret key.
Encapsulation picks s mod 3^F3, walks the 3^F3-isogeny with kernel P_D + s·Q_D
to E_C, and masks the torsion images with (a, b) = G(s, H(pk)).
Decapsulation runs the Kani chain on E_C × E_C with a fixed public step
schedule. It identifies ψ′ = α − c on E_C[3^F3], recovers ker φ̂ and then s by
pairings and a DLP in μ_{3^F3}. The Fujisaki–Okamoto transform with implicit
rejection gives the KEM.

## 2. Parameters, sizes and timings

Parameters come from `cargo run --release --example gen_params`. That Rust
port reproduces the former Sage script's output byte for byte.

| Set | target | p bits | limbs | E | F3 | d = 3^F3 | f bits | pk = ct bytes |
|---|---|---|---|---|---|---|---|---|
| Toy16 | none | 216 | 4 | 174 | 21 | 2^33.3 | 51 | 216 |
| Toy32 | none | 369 | 6 | 297 | 41 | 2^65.0 | 80 | 376 |
| L1 | 128 | 1341 | 21 | 1076 | 162 | 2^256.8 | 278 | 1344 |
| L3 | 192 | 1991 | 32 | 1594 | 243 | 2^385.1 | 409 | 1992 |
| L5 | 256 | 2630 | 42 | 2106 | 324 | 2^513.5 | 536 | 2632 |

Shared secret: 32 bytes. Public key and ciphertext are each 4 F_p² elements.

Wall time, one thread, measured in run `run-20261011T1813Z-bench` on an
Intel Xeon at 2.80 GHz (shared cloud VM, rustc 1.97.0, release build). The
records and the host description are in
`research/eigen_seta_rust_20261011/reports/run-20261011T1813Z-bench/`.
Each cell shows the range over the trials, with the trial count in
brackets. Every decapsulation returned the encapsulated key.

| Set | keygen (s) | encaps (s) | decaps (s) |
|---|---|---|---|
| Toy16 | 0.38 [1] | 0.005–0.008 [5] | 0.12–0.19 [5] |
| Toy32 | 2.3 [1] | 0.009–0.012 [3] | 0.45–0.59 [3] |
| L1 | 32.3 [1] | 0.24–0.26 [3] | 13.4–15.0 [3] |
| L3 | 118.6 [1] | 1.03 [2] | 51.8–53.9 [2] |
| L5 | 271.8 [1] | 2.17–2.20 [2] | 120.2–124.0 [2] |

Keygen is a single sample. Its rejection sampling makes the time highly
variable from key to key.

Decapsulation is slow because of how the chain is implemented. A callgrind
profile at Toy16 puts most of the time in the Richelot steps in Mumford
coordinates. Each step inverts a polynomial modulo a quartic for every point
it moves. That costs a few hundred F_p² multiplications per point, a figure
estimated from the formulas, not measured. Dimension-2 chains in theta
coordinates, as used by SQIsign2D and related work, avoid this polynomial
arithmetic. Porting the chain to theta coordinates is the main open
performance item. Its gain has not been measured here.

In the same run, the field layer measured as follows, in ns per operation.
Multiplication is the best of 15 runs and inversion the mean of 200.

| limbs | mul | inversion, Fermat | inversion, safegcd (used) |
|---|---|---|---|
| 4 (Toy16) | 54 | 18,403 | 4,741 |
| 6 (Toy32) | 72 | 67,657 | 10,316 |
| 21 (L1) | 611 | 1,426,022 | 35,608 |
| 32 (L3) | 1,516 | 6,708,390 | 109,506 |
| 42 (L5) | 2,449 | 13,453,052 | 96,327 |

## 3. Security, as conjectured

The full analysis, with proofs of correctness and an attack-cost table, is
the paper in `papers/eigen-seta/`. The design targets λ = 128, 192 and 256
classical bits, and the arguments are heuristic. No third party has analysed
the scheme.

* **Key recovery.** Knowing f gives a working secret key in polynomial time,
  so key recovery means finding f. The best attacks we know are these:
  * Exhaustive search over the admissible f. Keygen's public filters keep
    about 2^-8.1 to 2^-9.0 of the odd candidates, leaving about 2^268, 2^399
    and 2^525 values at L1, L3 and L5.
  * Finding any isogeny from E_D to E_A, after which lattice reduction yields
    the unique short isogeny of degree f. This costs at least p^(1/3) ≈ 2^447
    classically and p^(1/4) ≈ 2^335 quantumly at L1.
  * A meet in the middle on a smooth f, at cost at least √f ≈ 2^138 at L1.

  The structural assumption to scrutinise is that the eigenbasis, which
  reveals τ = εωε̂ on E_A[2^E] up to the unknown scalar f, gives no handle
  on f.
* **Message recovery, the binding constraint.** The ciphertext is a
  FESTA-type masked-torsion instance. Guessing the mask ratio modulo about
  √(3^F3) and then running Robert's dimension-8 attack costs about 2^128,
  2^192 and 2^256 trials classically. With Grover it costs the square root:
  2^64, 2^96 and 2^128 trials. Claw finding on the 3^F3-isogeny costs √d
  classically. Its degree exceeds that of the SIKE parameter sets assessed at
  the same categories.
* **NIST categories.** Categories 1, 3 and 5 are defined by AES key search,
  which is itself a 2^λ search classically and a 2^(λ/2) Grover search
  quantumly. The binding attack has the same shape, with a far more
  expensive oracle. We therefore conjecture categories 1, 3 and 5.

  This corrects the statement in the first version of this README, which
  said the quantum security fell below the category of the same name. The
  margin rests on the oracle cost. Any attack that needs fewer than about
  √d bits of mask information would break the claim.
* **Chosen-ciphertext security.** The FO transform with implicit rejection
  gives IND-CCA from one-wayness. The correctness error is heuristically
  O(E/p), and no decapsulation failure has been observed.
* **Earlier variants are broken.** The Sage prototype's `walk` variant leaks
  its trapdoor through a public degree. See section 4.

## 4. Corrections to the earlier merged material

* **The `walk` variant is insecure.** The README and the design note in
  `research/eigen_seta_prototype_20261011/` described it as a hidden CM curve
  behind a secret 5^k-walk. Its decryption key (s0, c, c′) is computed from
  Δ′ = 7·25^k alone, as `eigen_seta_toy.py` lines 133–155 show, and k is
  public. So anyone can compute the key and decrypt. An exploratory run, not
  archived, decrypted 5 of 5 toy ciphertexts this way. The `cm` variant was
  already labelled insecure. v2 fixes this by making the degree f secret.
* **Non-cyclic kernels.** The note (§5.4) reported ψ′ non-cyclic on
  E_C[3^f] in 1 of 8 toy trials at the 44-bit size, and its cause was never
  isolated. In v2 the case cannot occur. Proposition 3 of the paper shows
  that τ has no eigenline mod 3, because −7 is a non-square mod 3, so ψ maps
  E_C[3^F3] onto a cyclic group of order 3^F3.

## 5. Layout

| Module | Contents |
|---|---|
| `fp`, `safegcd`, `fp2`, `field` | Constant-time F_p (const-generic Montgomery, Montgomery-friendly reduction), Bernstein–Yang inversion, F_p² |
| `scalar` | Fixed-width scalars mod 2^E and mod 3^F3 |
| `mont`, `isog3`, `strategy` | x-only Montgomery ladders, y-recovery, 3-isogeny chains with optimal strategies |
| `weier`, `kani`, `chain` | Weierstrass points, Vélu, gluing, Richelot and splitting (Castryck–Decru formulas, MIT), the fixed-schedule decapsulation chain |
| `pairing` | Weil pairing of order 3^k, constant-time DLP in μ_{3^k} |
| `cm`, `ntheory` | Key-generation helpers (variable time): E_D, ω, torsion bases, Cornacchia, Miller–Rabin |
| `kem` | Key generation, encapsulation, decapsulation (public phase and constant-time core), FO transform |

## 6. Running

```
cargo test --release                        # KEM, KATs, robustness, chain, pairing, field tests
cargo test --release -- --ignored           # also Toy32 round trips and the L1 KAT
cargo run --release --example run -- L1 3   # keygen + 3 encaps/decaps at L1
cargo run --release --example bench_field   # F_p mul and inversion timings
cargo run --release --example ct_probe -- 1 1 valid L1   # constant-time probe (run under callgrind)
cargo run --release --example dudect_decaps -- 1000      # Welch t-test of decaps timing
cargo run --release --example f_filter_rate -- L1 200    # acceptance rate of keygen's public f filters
cargo run --release --example gen_params > src/params_generated.rs
```

`tests/kem.rs::kem_toy32` is `#[ignore]` because of its run time. Use
`cargo test --release -- --ignored` to run it.
