# eigen_seta: Rust implementation of the Eigen-Séta KEM and PKE

This is a research implementation of **Eigen-Séta v2**, the secret-degree
variant. It is an isogeny-based KEM whose trapdoor is a secret endomorphism α
of a public curve. Decryption is a "lollipop" endomorphism decoded by a
(2^E, 2^E)-isogeny chain in dimension 2, built with Kani's lemma. The design
and the earlier Sage toy prototype are in `RESEARCH_ENDRING_PKE_STRATEGY.md`
and `research/eigen_seta_prototype_20261011/`. Section 4 corrects two
statements made there.

**Status.** Research code. It decapsulates correctly at the nominal L1, L3
and L5 parameter sizes. Decapsulation is written to be constant time, and the
instruction-count check is in `SIDE_CHANNELS.md`. The security levels are
**conjectured** targets from the analysis below. No cryptanalysis by others
exists. Do not use this to protect data.

`pke::EigenSetaPke` adds message encryption using a fresh KEM encapsulation
per message and ChaCha20-Poly1305 authenticated encryption. Its wire format
is the fixed-size KEM ciphertext, a 12-byte nonce, then authenticated
ciphertext and 16-byte tag. The caller may supply associated data; decryption
rejects a changed context or ciphertext. This makes the in-memory API
end-to-end for plaintexts but does not upgrade the conjectural security of
the underlying KEM.

The API takes a `PublicKey`, plaintext, optional associated data and a
cryptographic RNG; the receiver supplies the same associated data:

```rust
use eigen_seta::kem::EigenSeta;
use eigen_seta::params::L1;
use eigen_seta::pke::EigenSetaPke;
use rand_core::OsRng;

let mut rng = OsRng;
let (pk, sk) = EigenSeta::<L1, 21>::keygen(&mut rng);
let ct = EigenSetaPke::<L1, 21>::encrypt(&pk, b"hello", b"context", &mut rng).unwrap();
let pt = EigenSetaPke::<L1, 21>::decrypt(&sk, &ct, b"context").unwrap();
assert_eq!(pt, b"hello");
```

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

The design targets λ = 128, 192 and 256 classical bits. The arguments are
heuristic.

* **Key recovery** needs α or f. The degree f is drawn from [f_lo, 2 f_lo),
  a range of about 2^277, 2^408 and 2^535 values at L1, L3 and L5. Every
  known way to use the public torsion images first needs f. Endomorphism-ring computation for E_A is the generic attack. Its
  cost is above 2^λ at these p sizes.
* **Message recovery** is a FESTA-type problem: recover a 3^F3-isogeny from
  masked torsion images. Guessing the masks or a meet-in-the-middle on the
  3^F3-isogeny costs about 2^λ classically. Grover-type search costs about
  2^(λ/2) quantum. That is the same exposure as FESTA, so the quantum
  security is not at the NIST level of the same name.
* **Earlier variants are broken.** The Sage prototype's `walk` variant leaks
  its trapdoor through a public degree. See section 4.

## 4. Corrections to the earlier merged material

* **The `walk` variant is insecure.** The README and the design note in
  `research/eigen_seta_prototype_20261011/` described it as a hidden CM curve
  behind a secret 5^k-walk. Its decryption key (s0, c, c′) depends only on
  Δ′ = 7·25^k, and k is public. An outsider recomputes the key from public
  data. The Sage experiment recovered s for 5 of 5 ciphertexts using only the
  public key and k. The `cm` variant was already labelled insecure. v2 fixes
  this by making the degree f secret.
* **Non-cyclic kernels.** The note (§5.4) reported ψ′ non-cyclic on
  E_C[3^f] in 1 of 8 toy trials at the 44-bit size, and the toy README
  listed it as a corner case. Its cause was never isolated, so the
  observation is unverified. It does not occur in v2: decapsulation here
  requires a full-order image of E_C[3^F3] under ψ′, and it has never failed
  at Toy16, Toy32, L1, L3 or L5.

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
cargo test --release                       # Toy16 KEM, chain, pairing, field tests
cargo test --release -- --ignored          # Toy32 KEM round trip
cargo run --release --example run -- L1 3  # keygen + 3 encaps/decaps at L1
cargo run --release --example pke_run -- L1 # encrypt/decrypt a message at L1
cargo run --release --example bench_field  # F_p mul and inversion timings
cargo run --release --example gen_params > src/params_generated.rs
```

`tests/kem.rs::kem_toy32` is `#[ignore]` because of its run time. Use
`cargo test --release -- --ignored` to run it.

Untrusted public keys are checked for a nonsingular curve, a consistent
differential point, and an independent full-order 2-power basis before
encapsulation. The public 3-power basis search is bounded so malformed
curve inputs cannot make it loop forever; malformed ciphertexts use
implicit rejection. `.github/workflows/eigen-seta.yml` runs formatting,
Clippy, the release tests including Toy32, and L1 KEM/PKE end-to-end smoke
tests.
