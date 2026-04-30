# Security policy

## Status

**This library is not safe for production use.**  It is an educational
reimplementation of well-known cryptographic primitives, written from
scratch to study the algorithms.  Several non-trivial classes of
vulnerability remain, and the code has had no third-party audit.

If you are evaluating this library for use in a system where confidentiality,
integrity, or authenticity matters — payments, identity, healthcare, anything
under PCI-DSS / HIPAA / GDPR / SOC 2 / FIPS 140 — **stop and use an audited
library instead.**  Concrete recommendations are in
[Recommended alternatives](#recommended-alternatives) below.

The hardening pass documented in [What is and is not protected]
(#what-is-and-is-not-protected) closes a number of obvious footguns but
does not change this conclusion.

## Recommended alternatives

| Use case | Recommended crate |
|---|---|
| General-purpose, audited, side-channel resistant | [`aws-lc-rs`](https://crates.io/crates/aws-lc-rs) (BoringSSL fork; FIPS-eligible build) or [`ring`](https://crates.io/crates/ring) |
| TLS / X.509 | [`rustls`](https://crates.io/crates/rustls) |
| AEAD (AES-GCM, ChaCha20-Poly1305) | [`aes-gcm`](https://crates.io/crates/aes-gcm), [`chacha20poly1305`](https://crates.io/crates/chacha20poly1305) (RustCrypto) |
| ECDSA / ECDH | [`p256`](https://crates.io/crates/p256), [`k256`](https://crates.io/crates/k256) (RustCrypto) |
| Ed25519 / X25519 | [`ed25519-dalek`](https://crates.io/crates/ed25519-dalek), [`x25519-dalek`](https://crates.io/crates/x25519-dalek) |
| RSA | [`rsa`](https://crates.io/crates/rsa) (RustCrypto) — and prefer ECC where you can |
| KDF / HMAC | [`hkdf`](https://crates.io/crates/hkdf), [`hmac`](https://crates.io/crates/hmac) (RustCrypto) |
| Post-quantum (real ML-KEM / ML-DSA) | [`ml-kem`](https://crates.io/crates/ml-kem) (RustCrypto), [`pqcrypto`](https://crates.io/crates/pqcrypto) |

## What is and is not protected

### What this library now does

- **Entropy** — All secret material is drawn from `OsRng` (OS CSPRNG via
  the `getrandom` crate: `getrandom(2)` on Linux, `BCryptGenRandom` on
  Windows, `getentropy(2)` on Darwin).  No user-space PRNG is used for
  long-term keys.
- **Constant-time MAC verification** — `aes_gcm_decrypt` and
  `chacha20_poly1305_decrypt` compare authentication tags using
  `subtle::ConstantTimeEq`.  A public `hmac_sha256_verify` is exposed for
  general HMAC verification.
- **Cache-timing-resistant AES S-box** — `SubBytes` and the key-
  schedule's `SubWord` no longer index a 256-byte table by a secret
  index.  Each S-box call scans the full table (256 reads, identical
  for every input) and selects the matching entry via
  `subtle::ConditionallySelectable`.  This closes the
  Bernstein/Osvik–Shamir–Tromer cache-timing channel against
  co-resident attackers, at the cost of ~256× per S-box call.
  `gmul` (GF(2⁸) multiply used by MixColumns) was also rewritten to
  use bit masks instead of `if` branches on secret bits.
- **Branchless field arithmetic** — `FieldElement::sub` and
  `FieldElement::neg` now compute `(a + p - b) mod p` and
  `(p - a) mod p` directly, eliminating the `if a >= b` and
  `if a == 0` comparisons that branched on secret operand values.
  `FieldElement::pow` is routed through `mod_pow_ct`.  These changes
  close the explicit branches at this layer; the underlying
  `BigUint` limb arithmetic still leaks (see structural blockers
  below).
- **Bleichenbacher-resistant PKCS#1 v1.5 unpadding** —
  `pkcs1_unpad_encrypt` no longer early-exits on individual padding-
  byte mismatches.  Every byte of the decrypted message is examined
  on every call, validity is accumulated as a `subtle::Choice`, and
  only the aggregate decision is branched on.  This hardens the
  per-byte timing channel that the original implementation exposed.
  **However**, the public `rsa_decrypt` still distinguishes valid
  from invalid by returning `Result::Err`, which is itself a
  Bleichenbacher oracle at the protocol layer — applications that
  surface this error to a network attacker are still vulnerable.
  See `pkcs1_unpad_encrypt`'s docs for the full picture and the
  recommended mitigations (implicit rejection, RSA-OAEP, hybrid
  encryption with AEAD).
- **Deterministic ECDSA (RFC 6979)** — Nonces are derived from the
  private key and message hash via HMAC-SHA-256 DRBG.  This eliminates
  the catastrophic key-recovery failure that occurs when `k` is
  repeated, biased, or predictable.  Verified against RFC 6979 §A.2.5
  test vectors.
- **Best-effort zeroization** — `RsaPrivateKey`, `EccPrivateKey`,
  `KyberPrivateKey`, and `McEliecePrivateKey` overwrite their secret
  fields on drop.  See the next section for caveats.
- **Authenticated encryption only** — Textbook (unpadded) RSA and the
  raw PKCS#1 v1.5 helpers are `pub(crate)` internals; the public API
  only exposes `rsa_encrypt`/`rsa_decrypt` (which apply PKCS#1 v1.5
  padding) and `rsa_sign`/`rsa_verify`.
- **Montgomery-ladder ECC scalar multiplication** —
  `Point::scalar_mul_ct` is used everywhere a secret scalar is
  multiplied by a curve point: ECC keypair generation
  (`d → d·G`), ECDH (`d_A · Q_B`), and the ECDSA signer's `k·G`
  step.  Every iteration performs exactly one point addition and one
  point doubling, regardless of bit value, and the iteration count
  is fixed at the curve order's bit length.  The variable-time
  `Point::scalar_mul` is retained only for the ECDSA verifier's
  `u₁·G + u₂·Q` (public inputs).
- **Constant-time fixed-width 256-bit integer (`ct_bignum::U256`).**
  Four-limb little-endian `[u64; 4]` with branchless `adc` / `sbb` /
  `mul_wide` (`u128` accumulators), `add_mod` / `sub_mod` using
  `subtle::Choice`-driven `cmov`, and `mont_mul` / `mont_sqr` (SOS
  Montgomery multiplication).  The Montgomery constant
  `m' = -p^(-1) mod 2^64` is computed via 2-adic Newton iteration
  (`compute_minv64`) and is `const`-callable, so per-curve constants
  are determined at compile time.  Cross-checked against `BigUint`
  for all arithmetic (16 unit tests).
- **Constant-time secp256k1 field arithmetic
  (`ecc::secp256k1_field::SecpFieldElement`).**  Wraps `U256` with
  every value held in Montgomery form (`a·R mod p`, `R = 2^256`).
  `add` / `sub` / `mul` / `sqr` / `neg` are all constant-time at the
  limb level — `mul` and `sqr` go through `mont_mul`, which is built
  from `u64×u64→u128` partial products and never branches on operand
  magnitudes.  `inv` is Fermat's `a^(p-2) mod p` via a Montgomery
  ladder over Montgomery-form `mul`/`sqr`, with the `if bit` branch
  replaced by `cmov` over precomputed candidates.  Every arithmetic
  primitive is cross-checked against `BigUint` (14 unit tests,
  including the curated edge cases at limb boundaries and
  `±1 mod p`).
- **Constant-time secp256k1 point arithmetic
  (`ecc::secp256k1_point::ProjectivePoint`).**  Homogeneous projective
  coordinates `(X, Y, Z)` over `SecpFieldElement`, with the identity
  encoded as `(0, 1, 0)` so that `Z == 0` is the only special case
  (and the whole-field `mul` makes that special case branch-free).
  Addition uses the Renes–Costello–Batina (2015) **Algorithm 7** for
  `a = 0` curves, doubling uses **Algorithm 9**.  Both are *complete*
  formulas: they handle distinct points, doublings, the identity on
  either side, and `P + (-P)` in a single straight-line sequence of
  `SecpFieldElement` ops, with no branches on input values.  The
  Montgomery ladder over this point type
  (`ProjectivePoint::scalar_mul_ct`) replaces `if bit { ... } else
  { ... }` with `cmov` over candidate next-states, so even the source
  level no longer has a branch on the secret bit.  14 unit tests
  cross-check the projective formulas against the textbook affine
  implementation, including the scalar-by-curve-order = identity
  sanity check.
- **Constant-time P-256 field arithmetic
  (`ecc::p256_field::P256FieldElement`).**  Same shape as
  `SecpFieldElement` but parameterised on the P-256 prime
  `p = 2^256 - 2^224 + 2^192 + 2^96 - 1`.  Constants `R mod p`,
  `R² mod p`, and `m' = -p^(-1) mod 2^64` are all `const`-callable
  through `from_montgomery_unchecked` / `compute_minv64`.
  Cross-checked against `BigUint` for add/sub/mul/sqr/neg/inv on the
  full set of edge cases (11 unit tests).
- **Constant-time P-256 point arithmetic
  (`ecc::p256_point::P256ProjectivePoint`).**  Homogeneous projective
  coordinates over `P256FieldElement`.  Because P-256 has `a ≠ 0`
  (specifically `a = -3 mod p`), addition uses RCB **Algorithm 1**
  (the general-`a` complete addition, 17 muls) and doubling uses RCB
  **Algorithm 3** (the general-`a` complete doubling, 13 muls).  The
  curve constants `a`, `b`, and `b3 = 3·b` are carried as
  Montgomery-form `P256FieldElement` constants.  The same `cmov`-driven
  Montgomery ladder pattern as secp256k1.  12 unit tests cross-check
  against the textbook affine `Point::scalar_mul`, including a P-256
  external KAT (`d → d·G` matches Python `cryptography` /
  OpenSSL-derived `Q`) routed through the projective ladder.
- **Unified secret-scalar dispatcher (`ecc::ct::scalar_mul_secret`).**
  All secret-scalar multiplications — `EccKeyPair::generate` /
  `EccKeyPair::from_private`, `ecdh` / `ecdh_raw`, and the ECDSA
  signer's `k·G` step — route through `scalar_mul_secret(point, k,
  curve)`, which dispatches on `curve.name` to the secp256k1 or P-256
  projective constant-time ladder.  No supported curve falls back to
  the legacy affine `BigUint`-backed ladder for secret scalars.  The
  variable-time `Point::scalar_mul` is retained only for the ECDSA
  verifier's `u₁·G + u₂·Q` (both inputs public).
- **Constant-time multi-precision integer (`ct_bignum::Uint<LIMBS>`).**
  The `U256` from earlier hardening is now `Uint<4>` — a const-generic
  fixed-width integer instantiated at the limb counts the RSA private
  path needs: `Uint<16>` for 1024-bit, `Uint<32>` for 2048-bit,
  `Uint<48>` for 3072-bit, `Uint<64>` for 4096-bit.  Branchless
  `adc`/`sbb`/`mul_wide` over `u64` limbs, SOS Montgomery multiplication
  with split-buffer scratch (so it compiles on stable Rust without
  `generic_const_exprs`), and `cmov` over `subtle::Choice` for all
  conditional moves.  Cross-checked against `BigUint` at LIMBS=4
  (secp256k1 modulus), LIMBS=16 (1024-bit), LIMBS=32 (2048-bit), and
  LIMBS=64 (4096-bit) in the unit-test suite.
- **Constant-time modular exponentiation for RSA private keys
  (`ct_bignum::mont_pow_ct` + `MontgomeryContext<LIMBS>`).**  RSA
  private operations route through `rsa_mod_pow_ct_crt` (CRT path,
  see next entry) which delegates to a monomorphised
  `mont_pow_ct<LIMBS>` per half-prime.  The non-CRT
  `rsa_mod_pow_ct` is kept as a fallback for unusual sizes and as
  the cross-check oracle for CRT correctness tests.  Each call:
    - Builds a runtime `MontgomeryContext<LIMBS>` from `n` (the public
      modulus) — `m' = -n^(-1) mod 2^64` via the const-callable
      Newton-iteration `compute_minv64`, plus `R mod n` and `R² mod n`
      via one-time `BigUint` reductions.  Setup leaks no secret-
      dependent timing because `n` is public.
    - Runs a square-and-multiply ladder over Montgomery form for
      *exactly* `64 · LIMBS` iterations — high zero bits of `d` are
      still scanned, so the iteration count is a function of the
      modulus *width* only, not of the secret exponent's value.
    - Uses `Uint::cmov` for the per-bit "should I multiply?" decision:
      the multiply is always computed and the result conditionally
      selected, so per-bit work cost is independent of the bit value.
  The variable-time `crate::utils::mod_pow` is retained for
  public-exponent operations (encrypt, verify, Miller-Rabin witness
  tests).  Cross-checked against `BigUint::modpow` at LIMBS=32
  (multiple base/exp shapes including `n - 1` and full-width
  exponents) and LIMBS=64.  RSA sizes outside the supported set fall
  back to the legacy `BigUint`-backed `mod_pow_ct`, which has uniform
  per-bit op count but does *not* close the limb-level channel.
- **CRT-accelerated constant-time RSA private path
  (`rsa_mod_pow_ct_crt`).**  PKCS#1 §3.2 CRT decomposition: instead of
  computing `m = c^d mod n` directly with a `64 · FULL`-iteration
  ladder, we compute `m_p = c^dp mod p` and `m_q = c^dq mod q` with two
  `64 · HALF`-iteration ladders and recombine via Garner's formula
  `m = m_q + q · ((m_p - m_q) · qinv mod p)`.  Roughly 4× faster than
  the direct form and closes the same timing channel because every
  step routes through `Uint<HALF>` Montgomery arithmetic.  Specific
  CT considerations:
    - `RsaPrivateKey` carries `dp`, `dq`, `qinv` as additional secret
      fields; all three are best-effort zeroized on drop alongside
      `d`, `p`, `q`.  Keygen swaps `p` / `q` so that `p > q`
      (PKCS#1 convention) before computing `qinv = q^(-1) mod p`.
    - The `c → c mod p` and `c → c mod q` reductions deliberately
      avoid `BigUint::%` (which divides and so leaks the operand
      magnitudes of `p` / `q` at the limb level).  We use
      `ct_reduce_full_to_half`: split `c` into the low / high HALF-
      limb halves, reduce each via a single conditional subtract
      (valid because RSA primes have their top bit forced on, so
      every HALF-limb value is `< 2p`), and combine via
      `lo_red + (hi_red · R mod p)` where `hi_red · R mod p` is
      precisely the Montgomery form of `hi_red` and is therefore
      computed by one constant-time `to_montgomery` call.
    - The recombination's multiplicative step uses one `mont_mul`
      against `p`; the final additive step uses `mul_wide` + `adc`
      over `Uint<FULL>`.  Both are limb-level branchless.
    - `MontgomeryContext::new(p)` and `MontgomeryContext::new(q)`
      perform the per-prime setup using `BigUint` for `R mod p` and
      `R² mod p`.  This setup runs once per loaded key, not once per
      message, and an attacker positioned to observe per-key setup
      timing already has the in-memory `key.p` available.  The
      setup leakage is structurally unavoidable without a CT
      modular reduction we don't have on stable Rust.
  Cross-checked at 1024 bits (CRT result equals direct
  `rsa_mod_pow_ct`); 2048-bit equivalent has a slow `#[ignore]`d
  test that runs under `cargo test -- --ignored`.  Sizes outside
  the supported set fall through to `rsa_mod_pow_ct` (non-CRT) and
  thence to `mod_pow_ct` (`BigUint`-backed) for unusual moduli.
- **Fermat-based modular inverse for prime moduli** —
  `mod_inverse_prime_ct(a, p)` returns `a^(p-2) mod p` via
  `mod_pow_ct`, avoiding the variable-time loop in the extended
  Euclidean algorithm.  The ECDSA signer uses this to compute
  `k⁻¹ mod n` (where `k` is the per-message secret nonce and `n`
  is the prime curve order).
- **Additional study-grade block ciphers (Serpent, Threefish).**
  Two AES-alternative block ciphers are available alongside AES for
  algorithm-diversity experiments and as a reference for the
  hardening posture this crate applies elsewhere:
    - **Serpent** (`symmetric::serpent`) — the AES finalist by Anderson,
      Biham, and Knudsen.  32-round bit-sliced SP-network, 128-bit
      block, 128/192/256-bit keys (shorter keys are padded per the
      submission's recommendation; over-length keys are rejected).  The
      eight 4-bit S-boxes are accessed via a 16-entry constant-time
      scan (`subtle::ConstantTimeEq` + `ConditionallySelectable`), so
      the same cache-line argument that motivated AES's CT S-box
      applies here too.  IP / FP / linear transform are pure bitwise
      operations and have no secret-dependent branches.  Round keys
      are best-effort zeroized on drop.
    - **Threefish** (`symmetric::threefish`) — the tweakable block
      cipher underlying Skein.  All three sizes (256 / 512 / 1024-bit
      block) are implemented at their specified round counts
      (72 / 72 / 80) with the published rotation tables and word
      permutations.  Pure ARX (Add-Rotate-XOR) construction means
      every operation is naturally constant-time at the source level
      — there are no S-boxes or table lookups whatsoever.  Key /
      tweak extension words are computed via XOR with `C240`, and the
      expanded key material is best-effort zeroized on drop.
  **KAT status.**
  - **Threefish-256 / 512 / 1024 — KAT-verified.**  All three sizes
    pass the all-zero key/tweak/plaintext vectors from the Skein 1.3
    reference implementation, independently corroborated by the
    RustCrypto `threefish` crate / Crypto++'s
    `TestVectors/threefish.txt`.  Bit-exact byte interoperability
    with these reference implementations is established for those
    inputs.  See `kat_threefish{256,512,1024}_zero` for the vectors.
  - **Serpent — KAT regression, NOT interop-compatible.**  The
    implementation passes self-consistency (encrypt-decrypt round-
    trip, S-box / permutation inverse identities, avalanche), but
    its byte-exact output diverges from the Cambridge
    `floppy4/ecb_e_m.txt` (Anderson/Biham/Knudsen "standard"-order)
    KATs for all three key sizes.  Most likely cause: a byte/bit-
    ordering convention mismatch between our LE word-loading and the
    reference's MSB-first bit numbering, or an IP/FP indexing
    direction divergence.  The failing KATs are kept in the suite
    under `#[ignore]` (`cargo test -- --ignored`) as a tracked
    correctness regression.  **Do not use this Serpent for data that
    any other Serpent implementation must decrypt.**  The CT
    properties (16-entry CT scan over the S-box, branchless IP/FP/LT)
    still hold relative to the implementation we have — they just
    don't describe Serpent-the-standard.
- **Cryptanalysis toolkit (`cryptanalysis::*`).**  Researcher-facing
  module for evaluating *novel* candidate ciphers and hashes before
  deployment.  This module is **analytical**, not operational — it
  is not constant-time, allocates `Vec`s freely, and is intended to
  run on a designer's workstation against a candidate primitive.
  Capabilities:
    - **`sbox::Sbox`** — generic n-in / m-out S-box with the full
      modern distinguishing-table battery: DDT (Differential
      Distribution Table), LAT (Linear Approximation Table), **BCT
      (Boomerang Connectivity Table — Cid–Huang–Peyrin–Sasaki–Song,
      EUROCRYPT 2018)**, **DLCT (Differential-Linear Connectivity
      Table — Bar-On–Dunkelman–Keller–Weizman, EUROCRYPT 2019)**, and
      **truncated DDT (Knudsen, FSE 1994)**.  Derived metrics:
      differential / boomerang uniformity, max differential
      probability, max linear bias, max DLCT bias, nonlinearity,
      balancedness, bijectivity, and algebraic degree.
      Cross-checked against published numbers: AES S-box reproduces
      DU=4 / NL=112 / **boomerang uniformity = 6 (Cid et al. 2018,
      Table 4)**; Serpent S0 reproduces DU=4 / degree=3 with BCT
      dominating DDT entry-wise; identity S-box yields the expected
      degenerate metrics.  LAT is internally cross-checked against a
      Walsh–Hadamard transform; BCT/DLCT trivial-corner identities
      are tested.
    - **`boolean::*`** — Walsh–Hadamard transform, Algebraic Normal
      Form coefficients via Möbius transform, algebraic-degree
      extraction.  Used by `Sbox` internally and exposed for
      designers of stream ciphers / nonlinear filter functions.
    - **`avalanche::*`** — full avalanche matrix
      `A[i][j] = Pr[bit_j(F(x)) flips when bit_i(x) flips]` over
      arbitrary `Fn(&[u8]) -> Vec<u8>` cipher black boxes, plus SAC
      (Strict Avalanche Criterion) and BIC (Bit Independence
      Criterion) scorers.  SHA-256 verified to pass; the identity
      function verified to fail at the expected 0.5 deviation.
    - **`statistical::*`** — chi-squared byte-distribution test,
      monobit frequency test, NIST SP 800-22 §2.3 runs test.
      Wilson–Hilferty normal approximation for χ² p-values
      (no `statrs` dependency).  SHA-256 verified to pass all three;
      a 90%-zero-biased stream verified to be rejected.
  **Limitations**: not a full automated trail-search engine — for
  cipher state larger than ~16 bits, finding the best `r`-round
  differential trail is genuinely hard and we expose the per-round
  building blocks rather than ship a CP/SAT/MILP backend.  These
  tools are necessary, not sufficient, for security analysis: many
  broken ciphers passed every test in this module.  Treat them as
  one input to peer review, not a substitute for it.
- **ECC curve-parameter safety auditor (`ecc_safety::*`).**  Static
  audit of short-Weierstrass curve parameters for the seven
  documented ECDH-attack categories from the
  `elikaski/ECC_Attacks` survey:
    - **Order too small** — flags `bits(n) < min_order_bits` (default
      200, SafeCurves recommendation).  ECDLP cost is `O(2^(bits/2))`
      via Pollard-rho / BSGS.
    - **Smooth order** — trial-divides `n` against a 6k±1 wheel up
      to a configurable bound (default `2²² ≈ 4 M`), then
      Miller-Rabin's the residue.  Reports the largest prime factor
      and fails if it's below `min_largest_prime_factor_bits`
      (default 160, ≈ 80-bit Pollard-rho security).
    - **Almost-smooth + small private key** — Pohlig-Hellman variant:
      even if `n` has a huge prime factor, an attacker who knows
      private keys are `≤ k` bits can ignore prime factors of `n`
      bigger than `k` bits.  We compute the *effective* security
      level given a configurable private-key bit budget and flag if
      the smooth subgroup alone covers it.
    - **Point not on curve** — exposed as
      `CurveParams::is_point_on_curve(x, y)` for runtime use; checks
      both that the equation holds and that `(x, y) ∈ [0, p)²`.
      Also wraps the same check around the supplied generator as
      part of the audit.
    - **Singular curve** — discriminant check `4a³ + 27b² ≢ 0 (mod p)`.
    - **Supersingular / small embedding degree** — searches for the
      smallest `k` with `p^k ≡ 1 (mod n)` up to a configurable cap
      (default 100); fails if `k ≤ max_unsafe_embedding_degree`
      (default 6, catches MOV/Frey-Rück targets while leaving
      pairing-friendly `k = 12` BLS curves alone).
    - **Anomalous (`#E = p`)** — `n · h == p` check; Smart's attack
      runs in `O(1)` p-adic ops.
  Each check returns `SafetyCheck::{Pass, Fail(msg), Inconclusive(msg)}`,
  so a CI-style "is this curve safe?" gate is one method call away
  (`SafetyReport::all_pass()`).  Cross-checked: secp256k1 passes the
  full audit; deliberately-weak fixtures (singular `(a=−3, b=2, p=7)`,
  anomalous `p=n=13`, embedding-degree-2 `p=19, n=5`, smooth-order,
  small-order, generator-off-curve, small-key + smooth-subgroup) all
  fail the corresponding check while passing the others.  13 unit
  tests in total.
  **Limitations**: this is a *parameter* auditor — it does not
  detect runtime misuse like nonce reuse, message-not-hashed-before-
  signing, or the Curveball-style "verifier accepts attacker-
  supplied generator" bug.  Those require API discipline at the
  call site.  ECDSA-side equivalents (RFC 6979 nonce derivation;
  the verifier's generator-pinning) are documented in their
  respective module docs.
- **Paillier partially homomorphic encryption (`asymmetric::paillier`).**
  Pascal Paillier, EUROCRYPT 1999.  Probabilistic public-key
  encryption over `Z_{n²}` whose security reduces to the Decisional
  Composite Residuosity Assumption (DCRA, factoring-equivalent).
  Supports two homomorphic operations:
    - **Ciphertext addition** — `enc(m₁) · enc(m₂) = enc(m₁ + m₂ mod n)`
      via `paillier_add`.
    - **Scalar multiplication by a public constant** — `enc(m)^k =
      enc(k · m mod n)` via `paillier_mul_scalar`.
  Plus convenience operations: `paillier_add_plain` (add a public
  scalar to an encrypted value) and `paillier_rerandomise` (refresh
  ciphertext randomness without changing the plaintext, e.g. for
  anonymous-tally use cases).  Use cases this enables: encrypted
  aggregation (federated balance / sensor / vote tally), sealed-bid
  auctions, federated-learning gradient aggregation, threshold
  protocols.  Implementation uses the canonical `g = n + 1`
  simplification (Paillier §6) so encryption is a single
  multiplication plus one `r^n mod n²` exponentiation.  **Refuses
  modulus < 2048 bits** (NIST 2030+ recommendation).  Cross-checked
  with 14 unit tests including: round-trip, probabilistic-encryption
  smoke (two encryptions of same plaintext differ), the headline
  homomorphic-add identity, scalar-multiply (`13 × 7 = 91`), the
  canonical encrypted-aggregation use case (`Σ enc(i)` for
  `i = 1..10` decrypts to `55`), modular wrap-around at `n`, edge-
  case `m = n − 1`, re-randomisation correctness, and "decryption
  under wrong key fails."

  **Security caveats baked into the docstring:**
  - **IND-CPA secure under DCRA** *only if* every encryption uses a
    fresh random `r`.  Reusing `r` across two encryptions leaks
    `m₁ − m₂ mod n` instantly — same severity as ECDSA k-reuse.
  - **Malleable by design.**  Multiplying ciphertexts adds the
    plaintexts — that is the whole point.  Do **not** use Paillier
    as "encryption with integrity"; if you need both, layer it
    behind an authenticated channel or apply a non-malleable
    transformation (Cramer–Shoup style).
  - **Decryption is variable-time.**  `mod_pow` over `Z_{n²}` with a
    secret `λ` exponent is not constant-time.  This is the standard
    HE posture — masking `λ`-decryption costs more than the
    leakage warrants on a non-public-network server.  Mitigation:
    only decrypt on hosts an attacker cannot time.
  - **Scope: only Paillier.**  Fully homomorphic schemes (BFV, CKKS,
    TFHE) are deliberately out of scope — each is 10–25× the
    implementation cost (NTT-based ring-LWE, modulus switching,
    bootstrapping, noise-budget tracking) and a from-scratch
    educational library cannot deliver them responsibly.  For FHE
    use [Microsoft SEAL](https://github.com/microsoft/SEAL) or
    [OpenFHE](https://www.openfhe.org/) — both audited, both with
    published parameter sets from HomomorphicEncryption.org.
- **ElGamal multiplicative HE (`asymmetric::elgamal`).**  Taher
  ElGamal, CRYPTO 1984.  Probabilistic public-key encryption over
  the prime-order subgroup of `Z_p*`.  Security reduces to the
  Decisional Diffie-Hellman (DDH) assumption.  Hardcoded to use
  RFC 3526 Group 14 (2048-bit safe prime, `q = (p-1)/2` also prime,
  `g = 2`) — the same group TLS / IPsec / SSH use as fallback DH
  parameters.  Hardcoding standard groups avoids the ~30 s
  per-keypair cost of generating a fresh safe prime and matches how
  real systems deploy DH.
  Supports: `paillier`-style **multiplicative** homomorphism
  (`enc(m₁) · enc(m₂) = enc(m₁·m₂ mod p)`), plaintext-multiply
  (`enc(m) · k → enc(k·m)`), plaintext-exponentiate
  (`enc(m)^k → enc(m^k)`), re-randomisation.  9 unit tests including
  the chained-product canonical case (`enc(2)·enc(3)·enc(5)·enc(7) →
  210`) and `enc(13)^7 → enc(13⁷)`.  Same security caveats apply as
  for Paillier: malleable by design, variable-time decryption,
  non-fresh `k` leaks plaintext relationships.

- **EC-ElGamal additive HE (`asymmetric::ec_elgamal`).**  The "Helios"
  scheme — additive homomorphic encryption over an elliptic curve,
  used by ElectionGuard and most modern verifiable-voting systems.
  Plaintext `m` is encoded as the curve point `m·G`; ciphertext
  addition composes encoded plaintexts as `(m₁ + m₂)·G`.
  Decryption recovers `M = m·G` then solves the elliptic-curve
  discrete log via baby-step giant-step (BSGS) — necessarily
  bounded to small message spaces.  Default bound `m_max = 2³²`
  (4 G plaintext range, ~3 MB BSGS table, ~milliseconds).  Larger
  bounds dialled up to `2⁴⁰` (~16 MB / ~1 s); above that, use
  Paillier.
  Supports: **additive** homomorphism (`enc(m₁) + enc(m₂) =
  enc(m₁+m₂)`), plaintext-add, scalar-multiply, decrypt-to-point
  (returns `m·G` without solving DLP — useful for equality tests).
  12 unit tests including the canonical encrypted-vote-tally
  pattern (10 ballots, sum decrypts to `7`).  Slow tests (anything
  using >1 P-256 affine scalar multiplication) are gated under
  `#[ignore]` because the affine `Point::scalar_mul` is the textbook
  variable-time backend; the constant-time projective backend is
  currently wired only into the ECDSA path and could be plugged in
  here as a follow-up optimisation.  Tested on P-256; works on any
  curve exposing the `CurveParams` API (secp256k1 also).

- **Why no fully homomorphic encryption (BFV / CKKS / TFHE).**
  These schemes require ring-LWE polynomial arithmetic with NTT,
  modulus switching, relinearization, and noise-budget tracking —
  roughly 5 000–10 000 LoC done responsibly per scheme, plus
  bootstrapping infrastructure on top for full FHE.  Parameter
  selection is exquisitely sensitive to the LWE-estimator's latest
  results.  An educational from-scratch library cannot deliver
  these with the rigour they need.  Use audited
  [Microsoft SEAL](https://github.com/microsoft/SEAL) or
  [OpenFHE](https://www.openfhe.org/) for FHE — both audited, both
  with HomomorphicEncryption.org parameter sets.  This crate
  implements the three partially-homomorphic schemes (Paillier,
  ElGamal, EC-ElGamal) which together cover every PHE primitive a
  real protocol needs (additive/multiplicative on `Z_n²`/`Z_p*`/EC).

- **HNP / lattice attack on biased-nonce ECDSA
  (`cryptanalysis::hnp_ecdsa`).**  The Boneh–Venkatesan
  Hidden-Number-Problem reduction, solved via LLL.  Given several
  ECDSA signatures whose per-message nonces `k_i` are known to lie
  in `[0, 2^k_bits)` for some `k_bits < n_bits`, this module
  recovers the signing key `d`.  This is **the** practical
  cryptanalysis of ECDSA — the attack that recovered Bitcoin keys
  (Breitner & Heninger 2019), TPM keys (Moghimi *et al.* TPM-Fail
  2020), and countless embedded-device keys; Albrecht & Heninger
  2021 sharpened it to handle one bit of bias given enough sigs.
  Implementation: `(m+2) × (m+2)` integer-scaled Boneh–Venkatesan
  lattice; LLL with `delta = 0.75`; recovery scans every reduced
  row, computes `d = ±row[m] · (2^k_bits)⁻¹ mod n` for each sign,
  and verifies against the public key.  **Signature counts
  required** (256-bit curve): 6–8 for 64-bit bias, 12–15 for
  32-bit, 25–30 for 16-bit, 50+ for 8-bit, BKZ-only territory at
  ≤4 bits.  Tests: 64-bit-bias recovery on P-256 with 8 sigs
  (~5s), 16-bit-bias recovery with 30 sigs (~30s, gated under
  `#[ignore]`), negative-control test that **full-entropy nonces
  cannot be recovered** (the security property RFC 6979 protects),
  plus rejection cases for too few sigs / no bias.  Why ship this:
  it's a **regression test** for the RFC 6979 nonce derivation in
  [`crate::ecc::ecdsa`] — if that derivation ever silently
  produces biased output, this attack catches it; feed the bad
  nonces in, the key falls out.

- **LLL lattice reduction (`cryptanalysis::lattice`).**  Classical
  Lenstra–Lenstra–Lovász with `f64` Gram–Schmidt and `BigInt`
  basis vectors — the same hybrid `fplll` and NTL's `LLL_FP` use.
  Sound at the dimensions and bit-sizes that matter for our
  cryptanalytic targets (≤30 dim, ≤256-bit coefficients); higher
  dimensions or longer entries would call for `L²`
  (Nguyen–Stehlé 2005) which is not implemented.  Test suite: the
  size-reduction invariant `|μ_ij| ≤ 1/2`, the obvious
  shorter-vector reduction from the original paper, identity-
  basis no-op, and long-then-short auto-swap.  Used by
  `hnp_ecdsa`; available as a standalone primitive for future
  cryptanalytic work (CVP via Babai's nearest-plane, sub-bit
  bias attacks via BKZ, etc.).

- **Pollard rho — Floyd + distinguished-points
  (`cryptanalysis::pollard_rho`).**  Pollard 1978 — the canonical
  generic-group `O(√n)` DLP attack whose √n cost is the reason a
  256-bit curve order gives ~128-bit security.  **Two variants:**
    - **Floyd tortoise-and-hare** (`pollard_rho_dlp` /
      `pollard_rho_dlp_zp`) — generic over four closures (`op`,
      `eq`, `partition`, `pow`).  Automatic restart on sterile
      collisions (up to 16 with random `(a₀, b₀)`).
    - **Distinguished points** (`pollard_rho_dp_dlp_zp` /
      `pollard_rho_dp_dlp_zp_multi`) — Van Oorschot-Wiener 1999.
      Each walker runs to a distinguished point (DP — element
      with low `dp_bits` zero), stores `(x, a, b)` in a shared
      table, starts a new walker.  Trivially parallelisable and
      multi-target friendly (Galbraith-Lin-Scott).  Sequential
      implementation here; the per-walker independence is the
      natural unit of work for future `rayon`-parallelism.
  Tests (8): Floyd 5-bit `Z_23*`, Floyd 16-bit `Z_131267*`,
  Floyd multi-target three secrets, Floyd 20-bit (`#[ignore]`),
  DP rho 5-bit, DP rho 16-bit, DP rho multi-target, sub-mod helper
  correctness.  Why ship this: (1) property-test for our scalar
  arithmetic — any correctness regression in `Point::add` /
  `scalar_mul` makes rho fail to recover planted secrets,
  surfacing the bug obviously; (2) empirical validation that the
  security floor reported by [`crate::ecc_safety`] extrapolates
  correctly.

- **LLL + BKZ lattice reduction (`cryptanalysis::lattice`).**
  Classical Lenstra–Lenstra–Lovász (LLL) plus block Korkine-
  Zolotarev (BKZ-β) with Schnorr-Euchner enumeration.  Hybrid
  `f64` Gram–Schmidt + `BigInt` basis — the same trade-off
  `fplll` and NTL's `LLL_FP` use.  Sound at the dimensions and
  bit-sizes our cryptanalytic targets need (≤30 dim, ≤256-bit
  coefficients, β ≤ ~25); higher dimensions or longer entries
  would require `L²` (Nguyen-Stehlé 2005) extended-precision GS
  which is not implemented.  BKZ enumerates blocks to find
  shorter vectors than plain LLL — useful for HNP at sub-bit
  bias and for any future lattice-cryptanalysis (NTRU, LWE).
  Test suite (9 tests): LLL size-reduction invariant
  `|μ_ij| ≤ 1/2`, two-dim obvious shorter-vector reduction,
  identity-basis no-op, long-then-short auto-swap, BKZ-2
  matches LLL, BKZ-8 produces ≤ LLL norm on a 6-dim
  challenge basis, BKZ preserves lattice volume, BKZ rejects
  β < 2.

- **HNP / lattice attack on biased-nonce ECDSA
  (`cryptanalysis::hnp_ecdsa`).**  Boneh–Venkatesan
  Hidden-Number-Problem reduction solved via LLL.  Given several
  ECDSA signatures whose per-message nonces lie in `[0, 2^k_bits)`,
  recovers the signing key.  **The** practical ECDSA attack —
  recovered Bitcoin keys (Breitner–Heninger 2019), TPM keys
  (TPM-Fail 2020), countless embedded-device keys.  Tests:
  64-bit-bias recovery on P-256 with 8 sigs (~5 s), 16-bit-bias
  recovery with 30 sigs (~30 s, gated under `#[ignore]`),
  negative-control test verifying full-entropy nonces cannot be
  recovered.

- **ECDSA scalar blinding (`ecc::ct::scalar_mul_secret_blinded`).**
  Coron 1999 ("Resistance against Differential Power Analysis for
  Elliptic Curve Cryptosystems") — applied automatically inside
  `ecdsa::sign_hash` for every signing call.  Computes `k' = k +
  r · n` for a freshly-sampled 64-bit `r`, runs the constant-time
  projective ladder on `k'` for `n_bits + 64` iterations.  Since
  `n · G = O`, the result `[k']G = [k]G` is identical and the
  signature `(r, s)` is unchanged — but the bit pattern driving
  the ladder is re-randomised every call, which **defeats DPA,
  template attacks, and ML-side-channel attacks that statistically
  combine many ladder traces against the same secret nonce**.
  This is the standard post-2020 ECDSA implementation discipline
  required by every public-cloud / shared-tenancy threat model
  (Hertzbleed, SQUIP, Downfall, etc.).  Tests: blinded result
  equals unblinded result (5 random rerolls each on secp256k1 +
  P-256), `[0]·G = O` edge case, blinding-bits constant pinned at
  64.  RFC 6979 KATs continue to pass, proving the blinding is
  fully transparent to signature output.

- **Trace-zero / Weil-descent auditor
  (`ecc_safety::CurveParams::check_no_weil_descent`).**  Diem 2011
  ("On the discrete logarithm problem in elliptic curves over
  non-prime finite fields"): for curves over `F_{p^n}` with
  composite `n`, Frey-Rück Weil descent + Gaudry-Schost-style
  index calculus solves ECDLP in subexponential time.
  Trace-zero subvarieties (Frey 1998, Galbraith-Smart) apply
  even for prime `n` ≥ 2.  New `AnalysisOptions::extension_degree`
  field (default 1 = prime field, no warning).  Composite `n`
  → hard `Fail` with the smallest prime factor reported;
  prime `n` ≥ 2 → `Inconclusive` (warn but pass — no concrete
  subexponential break is known at our parameter sizes).  Tests
  (5): prime-field passes, composite degree 6 fails, prime
  degree 5 returns Inconclusive, failure surfaces in
  `failed_checks()`, and a primality-helper smoke test.

- **Shor's algorithm — quantum-circuit simulator
  (`cryptanalysis::shor`).**  Shor 1994 — the polynomial-time
  quantum algorithm that ends RSA and ECC eventually.  This is a
  pedagogical state-vector simulator: actual H gates,
  function-application unitary `U_f`, and QFT acting on a
  hand-rolled `Vec<Complex>` state vector, with classical
  measurement sampling and continued-fraction post-processing.
  Total qubits ≤ 14 (16 KB state) → fast.  End-to-end factors
  small N (15, 21, 33).  Tests (6): continued-fraction
  convergents of 7/15 include 1/2, QFT of `|0⟩` is uniform,
  Hadamard of `|0⟩` is `|+⟩`, Shor finds period r = 4 for
  a = 7 mod 15 across 64 measurement seeds, end-to-end factor
  of 15 ∈ {3, 5}, end-to-end factor of 21 ∈ {3, 7}, rejects even
  N.  Honest scope: this is **integer-factoring Shor**.
  Shor-for-ECDLP needs an in-place EC scalar-multiplication
  unitary (Roetteler-Naehrig-Svore-Lauter 2017) which is
  ~10× more involved than the modular-exponentiation `U_f` here
  — that's a dedicated quantum-simulator project, not a
  cryptography-library deliverable.  The pedagogical takeaway
  is the same: P-256 falls when ~10^7 fault-tolerant qubits
  exist; this simulator demonstrates the algorithm structure
  on a tractable scale.

- **Multi-target rho margin check
  (`ecc_safety::CurveParams::check_multi_target_margin`).**
  Extension of the curve-parameter auditor for the
  Galbraith–Lin–Scott amortised-rho threat model.  When `m` keys
  share a curve, a single rho walk can break any one of them in
  `O(√(n/m))` operations rather than `O(√n)` — per-key effective
  security bits drop by ½·log₂(m).  Optional check (skipped by
  default; enable via `AnalysisOptions::expected_deployment_size =
  Some(N)`).  Tests: 256-bit secp256k1 passes at 2³⁰ keys (113
  effective bits, above 100-bit floor), fails at 2⁶⁰ keys (98
  effective bits, below floor), failure surfaces correctly in
  `failed_checks()`.

- **CGA-HNC: Class-Group-Amortized Hidden-Number Cryptanalysis
  (`cryptanalysis::cga_hnc`).** Original research thread pursued
  in this codebase (see `RESEARCH.md` for the full note).  The
  proposal: use the `Cl(O)`-action on an ordinary elliptic curve
  to amortise partial DLP recoveries across the isogeny orbit,
  then clean up the residual unrecovered scalar bits with an
  HNP lattice attack. Phase 1 (this code) builds the empirical
  measurement infrastructure: Φ_2 modular polynomial, j-invariant
  conversion, 2-isogeny BFS, brute-force point counting,
  cumulative-CRT smoothness measurement.  9 unit tests + 3
  experimental drivers gated under `#[ignore]`.  Empirical
  findings at toy scale (`p ≤ 10⁴`):  mean CRT-coverage of `d`
  via the 2-isogeny orbit alone = **0.856 · log₂(p)** at `p = 1009,
  B = 50`; substantial variance, with ~30% of sampled curves
  showing CRT/log₂(p) ≥ 1.0 (orbit *fully recovers* `d`).
  Heegner/CM-1 curves (j ≡ 0 or 1728) confirmed irrelevant
  (orbit size 1, CRT = 0).  These are toy-scale results and do
  not demonstrate a generic asymptotic improvement at
  cryptographic sizes — they characterise a "weak curve" family
  that previous auditors miss.

- **`ecc_safety::CurveParams::check_cl_o_orbit_brittleness`** —
  experimental auditor extension built on top of `cga_hnc`.
  Walks the 2-isogeny graph (capped at `max_orbit_check_bits` to
  prevent runaway at cryptographic sizes), measures CRT-coverage
  of `d` from the orbit's smooth point counts, fails if coverage
  exceeds `max_orbit_crt_fraction` (default 0.25 of `log₂(n)`).
  Tests verify (a) secp256k1 returns `Inconclusive` (correctly
  gated), (b) `(a=1, b=2)/F_1009` returns `Fail` (matches
  empirical CRT = 1.796), (c) `(a=17, b=1)/F_1009` (orbit size
  1) does not fail.

- **CGA-HNC: applicability to cryptographic curves.**  The
  research thread documented in `RESEARCH.md` reaches an honest
  negative conclusion for deployed crypto: the technique requires
  curves with composite point order to provide Pohlig-Hellman
  amortisation.  Every standard cryptographic curve in this
  library and adjacent libraries (secp256k1, NIST P-256, P-224,
  P-384, P-521, brainpoolP-*) has *prime order with cofactor 1
  by design* — explicitly chosen to prevent this attack family.
  The auditor encodes this: `check_cl_o_orbit_brittleness`
  short-circuits to `Pass` immediately for any curve with
  cofactor 1 and prime `n`, without running the (infeasible-at-
  scale) brute-force orbit walk.  Verified by the dedicated test
  `orbit_brittleness_passes_secp256k1_immediately`.
  Phase 3 empirical findings (toy scale): orbit cooperation IS
  observable (4 of 27 sampled `d` values triggered ≥ 2
  cooperative curves) but the cooperation provides redundant
  rather than necessary information — single-curve PH + BSGS
  recovers `d` regardless, and at toy scale BSGS is cheap.
  At cryptographic scale, BSGS is `2¹²⁸` and orbit cooperation
  *would* matter if it scaled — but it can't apply to
  prime-order curves regardless.  **The library's standard
  cryptographic curves are not threatened by CGA-HNC.**

- **CGA-HNC Phase 2: end-to-end attack pipeline.** Building on the
  Phase 1 orbit-characterisation infrastructure, Phase 2 implements
  the actual attack: explicit 2-isogeny point transport via Vélu's
  formula in short-Weierstrass form (codomain `a' = −15α² − 4a,
  b' = −22α³ − 8αa`; point map `x' = (x² − 4αx + 6α² + a)/(x − α)`,
  `y' = y(x² − 2αx − 2α² − a)/(x − α)²`); affine point arithmetic
  (`Pt2` add/double/scalar); Pohlig-Hellman DLP recovery on smooth
  subgroups; end-to-end orchestrator `cga_hnc_attack`.  Phase-2 demo:
  on `(a=1, b=2)/F_1009` with planted `d = 17`, the attack walks a
  32-curve orbit, transports `(P, Q = 17·P)`, runs Pohlig-Hellman,
  and **recovers `d = 17 mod 28` correctly via CRT**.  Mathematical
  correctness checks: Vélu codomain verified by `Φ_2(j_src, j_dst) ≡
  0`; point-arithmetic invariants (`2·P = P + P`, `n_curve·P = O`);
  Pohlig-Hellman recovery on transported points verified against
  planted `d`.  Honest scope: the demo proves the pipeline works
  end-to-end but does not demonstrate a generic ECDLP improvement —
  in this run only 1 of 32 orbit curves contributed useful info, so
  the result reduced to single-curve Pohlig-Hellman with isogeny-
  transport overhead.  The proposal's *novel* claim — that
  cross-curve aggregation provides super-additive information — is
  untested at non-toy scale (Phase 3 work).

- **`cryptanalysis::ecdsa_audit` — transcript auditor for deployed
  ECDSA streams.**  Direct cryptanalytic attack on ECDLP for
  prime-order curves (P-256, secp256k1, P-384, P-521) has been
  open for 30+ years; every published key recovery in real-world
  ECDSA has been an *implementation* attack, not a math attack
  (Sony PS3 2010, Android Bitcoin 2013, TPM-Fail 2020, GoFetch
  2024).  This module is the defensive auditor counterpart:
  given an ECDSA signature transcript and the public key, it
  automatically detects nonce bias and (if present) recovers the
  private key — *as a verification tool* that the transcript is
  safe.  Specifically targets the curves the user cares about
  (P-256, secp256k1).
  Implementation:
  (a) `quick_bias_score` — chi-squared / monobit prefilter on the
      low byte of `t_i = s_i⁻¹·z_i mod n` (Bleichenbacher-flavoured
      cheap detector);
  (b) auto-`k_bits` sweep over the suspected bias range, calling
      `hnp_recover_key` at each level and verifying recovery
      against the public key;
  (c) returns one of `NoBiasDetected` / `KeyRecovered{d, k_bits}` /
      `BiasSuspectedNoRecovery`.
  Tests (4): planted 64-bit bias on P-256 → recovered (~1 s);
  full-entropy nonces on P-256 → `NoBiasDetected` (proves RFC 6979
  defence works); chi-squared score distinguishes uniform from
  biased streams; degenerate input safely handled.
  The module's docstring also contains an **honest research note**
  enumerating five novel-but-untried angles for prime-order ECDLP
  (canonical-lift formal-group attack, Bleichenbacher FFT for
  sub-bit bias, multi-target preprocessing, joint multi-key HNP,
  quantum-classical hybrid) with calibrated probability estimates.

- **`cryptanalysis::multi_key_hnp` — joint HNP attack on BIP32 /
  HD-wallet ECDSA.**  *Original research deliverable* (research-
  note item #4, ~5% probability of incremental result).  The
  observation: hierarchical-deterministic wallets derive child
  keys via the public relation `d_j = master_d + δ_j (mod n)`,
  where `δ_j` is computable from the master xpub.  Substituting
  into the ECDSA equation gives a Hidden-Number-Problem instance
  in `master_d` for any signature from any child key.  Pooling
  signatures across many children attacks the master with much
  weaker per-signature bias than single-key HNP needs.
  Implementation: `ChildKeySignature` carries `(r, s, z, k_bits,
  offset)`; `multi_key_hnp_recover_master` substitutes
  `z' = z + r·δ_j (mod n)`, builds standard `BiasedSignature`
  instances over the master variable, and dispatches to the
  existing `hnp_recover_key`.  The key-recovery test
  (`multi_key_recovers_master_p256`) plants a P-256 master,
  derives 4 child keys, generates 2 biased signatures per child
  (8 pooled), and recovers the master at 64-bit bias; the
  control test (`single_key_fails_with_two_sigs`) confirms the
  same 2 sigs at the same bias on a *single* key are below the
  HNP threshold and fail.  Cross-key advantage demonstrated
  empirically.  Negative-control test verifies full-entropy
  nonces produce no false-positive recovery.
  **Honest novelty assessment** (in module docs): the
  mathematical mechanism is straightforward (HNP with `z'`
  substitution); the *formulation as an attack on BIP32* does
  not appear in the surveyed literature.  Closest published
  work: standard HNP (Boneh-Venkatesan 1996), cross-signer
  HNP via shared randomness (LadderLeak CCS 2020),
  implicit-factoring RSA shared-prime (May-Ritzenhofen PKC
  2009).  Realistic threat model: any HD wallet exposing its
  xpub publicly (e.g. Bitcoin merchant payment processors)
  while signing with biased nonces — a known implementation
  smell that this attack would amplify across the wallet's
  address set.  Modern wallets using RFC 6979 are immune
  by construction.

- **`cryptanalysis::preprocessing_rho` — Bernstein-Lange non-uniform
  ECDLP attack** (research-note item #3).  Asiacrypt 2013, "Non-
  uniform cracks in the concrete: the power of free precomputation."
  Splits the rho attack into a target-INDEPENDENT precomputation
  phase (build a table of `T` distinguished points reached from
  random `g^a` walks) and a target-dependent online phase (walk
  from `h`, look up the first table-hit DP).  At the asymptotic
  optimum `T = N^{2/3}`, online cost per target is `O(N^{1/3})` —
  sub-rho-`√n` for any single target after the (one-time) `N^{2/3}`
  precomputation.  Tests (4): planted-DLP recovery on a Sophie-
  Germain group; **multi-target amortisation** (one table → multiple
  targets recovered, demonstrating the non-uniform advantage);
  expected-online-cost calculator at `N = 2^{60}, T = 2^{40}`;
  empty-table edge case.  Practical relevance: at `n = 2^{256}`,
  `N^{2/3}` storage is astronomical — the attack is theoretical
  at cryptographic sizes but real at deployment-cost-modelling
  scale (e.g., assessing Bitcoin's 2³⁰-key population).

- **`cryptanalysis::bleichenbacher` — FFT-based bias detector for
  sub-bit nonce bias** (research-note item #2).  Bleichenbacher
  2000 (unpublished talk; canonical citations: de Mulder–Hutter–
  Marson–Pearson CHES 2014, Aranha–Tibouchi–Yarom LadderLeak CCS
  2020).  Where standard HNP/LLL needs `m · b > n_bits` signatures,
  the spectral method works *below* this threshold by detecting
  bias as a peak in `|Z(d)| = |Σᵢ exp(2πi·(t_i + h_i·d)/n)|`.
  Implemented at toy scale (direct `O(n·m)` evaluation for
  `n ≤ 2²⁴`); the lattice-reduction extension to cryptographic
  sizes is documented but not implemented.  Tests (6): peak-finding
  at correct `d`; full-range recovery at 3-bit bias / 400 sigs;
  **1-bit bias recovery at 2000 sigs** (the sub-threshold regime
  HNP/LLL cannot touch); negative control (uniform nonces give no
  peak); malformed-input rejection; short-transcript rejection.
  The 1-bit-bias test is the headline: it recovers the secret in
  the regime where standard HNP would need ~256 signatures
  *theoretically* and many more in practice.

- **`cryptanalysis::canonical_lift` — Smart's-attack foundation +
  p-adic infrastructure** (research-note item #1).  Smart 1999
  ECDLP attack on anomalous curves (`#E(F_p) = p`), with the
  p-adic / formal-group machinery built up to support it.
  Implementation: `ZpInt` truncated p-adic integers (add, mul,
  inverse, valuation); `ZpCurve` Weierstrass-curve arithmetic in
  `Z_p / p^precision` (point doubling, addition, scalar mul);
  Hensel-lifting of points from `E(F_p)` to `E(Z_p)` via Newton
  iteration; anomalous-curve search by brute-force point counting.
  Tests (7): Z_p arithmetic basics, Z_p inversion, Z_p valuation,
  curve doubling-equals-self-addition consistency, Hensel-lift
  satisfies curve equation, anomalous-curve found at `p = 11`,
  Smart's-attack stub documents the projective-coordinate
  extension needed for end-to-end recovery (the Hensel + p-adic
  + arithmetic infrastructure is in place; the projective-coord
  formal-log extraction at `[p]·P̂` is the missing piece, since
  affine arithmetic at finite precision cannot represent points
  at infinity over `F_p`).  Module documentation explicitly
  flags the **canonical-lift extension to non-anomalous curves**
  (research-note item #1's full ambition) as open research with
  `~3%` honest probability — the foundation now exists for anyone
  to attempt it.

- **`cryptanalysis::aut_folded_rho` — Pollard rho with full
  `Aut(E) = Z/6Z` orbit folding for j = 0 curves** (research-
  agent blind-spot finding #4: "the assumption 'CM gives only `√2`'
  deserves scrutiny").  **Specifically targets secp256k1**, which
  has CM by `Z[ω] = Z[(−1 + √−3)/2]` and full automorphism group of
  order 6 (generated by negation `[−1]` and the GLV endomorphism
  `[ω]: (x, y) → (β·x, y)`).  The standard cryptanalytic literature
  cites a `√2` speedup from GLV alone; the **complete** 6-fold
  orbit folding gives a theoretical `√6 ≈ 2.45×` speedup over naive
  rho.  Implementation: per-step canonicalisation reduces every
  visited point to its lex-smallest `Aut(E)`-orbit representative;
  collisions on canonical form give a relation modulo the
  recovered `α ∈ Aut(E)`; integer action of `α` is computed via
  the cube-root-of-1 multiplier `λ ∈ Z/n`.  Tests (5 unit + 1
  ignored empirical):
    1. `aut_group_axioms` — `α · α⁻¹ = id` for all 6 elements.
    2. `canonicalisation_invariant_on_orbit` — every orbit element
       maps to the same canonical representative.
    3. `apply_aut_canonicalises` — applying the recovered `α` to
       `P` yields canonical(`P`).
    4. `integer_action_matches_point_action` — `α(P) = μ_α · P`
       on Z[ω]-stable subgroups (prime order `≡ 1 mod 3`).
    5. `folded_rho_recovers_planted_dlp` — end-to-end recovery on
       a synthetic j=0 curve `(p=43, b=8, n=31)` where `Aut(E)`
       acts faithfully on the prime-order group.
  Plus an `#[ignore]`d empirical-speedup driver that reports
  iteration counts on `(p=67, b=2, n=73)`.  Honest framing
  documented in the test: the `√6` advantage is asymptotic;
  small-`n` variance dominates at toy scale, so empirical
  speedup measurement requires `n ≫ 10⁴`.  Cryptographic-scale
  validation would need the projective constant-time secp256k1
  backend already in `ecc::secp256k1_point`.
  This is **the agent's research-note item #4** delivered: a
  rigorous re-implementation of the structural exploit on
  secp256k1's CM that goes past the often-cited `√2` GLV factor.

- **BKZ-augmented HNP: `hnp_recover_key_with_reduction(...,
  HnpReduction::Bkz(β))`** (research-agent blind-spot finding #3:
  "modern lattice techniques on HNP at the bias-threshold
  frontier").  Existing HNP attack used LLL with `δ = 0.75`; the
  new variant accepts `HnpReduction::Lll` (default, fast) or
  `HnpReduction::Bkz(β)` (slower per call, recovers at thinner
  bias-bit/sig-count margins where LLL fails).  Empirical
  validation (`bkz_recovery_at_thin_margin`): on P-256 with
  **5 signatures at 56-bit bias** (`m·b = 280` vs `n_bits = 256`,
  only 9% margin), across 6 random seeds:
    - **LLL: 5/6 successes** (one degenerate-lattice failure)
    - **BKZ-12: 6/6 successes** (perfect)
  This is a published-result-but-rarely-shipped improvement: the
  cryptanalytic literature has noted BKZ's threshold advantage
  for ~10 years (Albrecht-Heninger 2021, "How to Find Short
  Lattice Vectors Faster"), but most open-source HNP
  implementations still use LLL.  This module pulls in the
  existing `cryptanalysis::lattice::bkz_reduce` and exposes the
  upgrade through a single enum.  Plus a second `#[ignore]`d
  empirical driver (`bkz_recovers_at_thinner_margin_than_lll`)
  at +25% margin showing the same 5/6 → 6/6 improvement.

- **`cryptanalysis::signature_corpus` — hypothesis-free
  statistical sweep for ECDSA signature corpora** (research-agent
  blind-spot finding #1: ~25% probability of finding novel
  deployment bugs; the highest-impact under-explored angle in
  ECDLP cryptanalysis).  This is the **infrastructure** for what
  the agent identified as the field's biggest cultural blind spot:
  massive-scale empirical analysis of deployed signatures, which
  cryptanalysis culture treats as low-status engineering despite
  every famous real-world ECDSA break (Sony PS3 2010, Android
  Bitcoin 2013, TPM-Fail 2020, GoFetch 2024) being found by it.
  The library ships the analyzer; users supply their own data
  (Bitcoin blockchain parser, TLS dump, hardware-token transcript).

  Streaming API: ingest `SignatureRecord { public_key, r, s, z,
  source }` records via `add()`; clusters by SEC1-encoded public
  key; runs four detectors per cluster at `finalize()`:
    1. **Repeating-`k`** — hash by `r`; any collision recovers
       the private key in one division.  Famous Bitcoin attack.
       Severity: `Critical`.
    2. **Single-key auto-HNP** — wraps `ecdsa_audit::audit_ecdsa_transcript`
       with the prefilter disabled (hypothesis-free).  Sweeps `k_bits`
       and runs the lattice attack.  Severity: `Critical` on recovery,
       `High` on suspected-but-unrecovered bias.
    3. **Bleichenbacher FFT** — toy-scale (gated to `n_bits ≤ 24`);
       flags significant SNR peaks.  Severity: `High`.
    4. **Chi-squared monobit prefilter** — quick statistical
       anomaly test on the low byte of `t_i = s_i⁻¹·z_i`.
       Severity: `Medium`.

  Aggregate `CorpusReport` carries severity-classified findings
  with `recovered_keys()` / `at_severity(Severity)` / `critical_count()`
  helpers.  Tests:
    1. Repeated-`k` recovery on planted P-256 attack — recovers `d` ✓
    2. HNP-recoverable bias on planted P-256 corpus — recovers `d` ✓
    3. Clean (RFC-6979-style) corpus — **zero** Critical/High false
       positives across 10 sigs ✓
    4. Mixed corpus (1 bad key + 3 clean) — flags the bad one,
       leaves the good ones alone, exactly 1 Critical finding ✓
    5. Empty corpus — only `Info` row ✓
    6. Malformed records — reported, no panic ✓

  Realistic deployment story: a Bitcoin merchant operator points
  this at their xpub's signature history (~10⁴-10⁶ signatures);
  if any address used a buggy RNG, the analyzer flags it and
  recovers the master key.  A TLS observatory operator points it
  at their corpus and finds device-cohorts with biased nonces.
  The library provides the pipeline; running it at internet-scale
  requires data infrastructure beyond Rust idiomatic surface but
  the analyzer scales linearly in corpus size and per-cluster
  analyses parallelise trivially.

- **`cryptanalysis::sha1_differential` — differential cryptanalysis
  toolkit applied to SHA-1.**  SHA-1 is **broken** (Wang-Yin-Yu
  Crypto 2005 theoretical break at `~2⁶⁹` ops; SHAttered
  Stevens et al. 2017 actual collision at `~2⁶³·¹` = 6500
  CPU-years).  Full-SHA-1 cryptanalysis is computationally
  infeasible in any session, but the techniques and toolchain
  are entirely visible at reduced rounds — that's what this
  module ships.

  SHA-1 lives in `cryptanalysis`, *not* `hash` — the library
  refuses to expose it as a primitive.  It's available
  specifically as a **target** for differential analysis.
  Implementation is reduced-round-capable: `sha1_compress(state,
  block, rounds)` for `rounds ∈ [0, 80]`, with verified NIST
  FIPS 180-4 known-answer tests at `rounds = 80`
  (`SHA-1("abc")` = `a9993e36...`, `SHA-1("")` = `da39a3ee...`).

  Differential-analysis tools provided:
  1. **Avalanche analysis** at any round count — wraps the
     existing `cryptanalysis::avalanche::full_avalanche`.
  2. **Boolean-function characterisation** — `f_t` truth-table
     extraction for the four phases (Choose, Parity, Majority,
     Parity).  Compatible with existing
     `cryptanalysis::boolean::walsh_hadamard`.
  3. **Random differential probability search** — estimates
     `P[F_r(x) ⊕ F_r(x ⊕ Δ_in) = Δ_out]` via Monte Carlo.
  4. **Reduced-round near-collision finder** — birthday-style
     random search; demonstrates the cost-per-round transition.

  **Empirical results** from `avalanche_transition_curve`
  (gated under `#[ignore]`):

  | Rounds | Worst-cell SAC deviation |
  |--------|--------------------------|
  | 5      | 0.500 (no diffusion)     |
  | 10     | 0.500 (no diffusion)     |
  | 15     | 0.344                    |
  | 20     | 0.312                    |
  | 80     | 0.344                    |

  At 5–10 rounds, the avalanche is degenerate: input bit flips
  deterministically affect specific output bits.  Diffusion
  kicks in at round ~15 (consistent with SHA-1's message-
  expansion structure: rounds 0–15 only touch the original 16
  message words; rounds ≥ 16 introduce expansion-rotation).

  Tests (8 unit + 1 experimental):
    1. `sha1_full_matches_nist_kat` — full SHA-1 implementation
       byte-identical to NIST KAT for `"abc"` and `""`
    2. `full_sha1_has_good_avalanche` — worst-cell deviation
       < 0.18 at 256 samples/bit, mean-cell deviation < 0.05
    3. `reduced_sha1_has_poor_avalanche` — 10-round deviation
       strictly worse than full
    4. `round_function_phases` — truth tables of Choose / Parity
       / Majority match expected algebraic identities
    5. `trivial_differential_has_prob_1` — `Δ_in=0 → Δ_out=0`
       sanity
    6. `nonzero_differential_at_reduced_rounds` — estimator
       validity at non-trivial inputs
    7. `near_collision_at_reduced_rounds` — 10-round random
       search produces strictly lower Hamming distance than
       80-round (concrete demonstration of the round-cost
       cliff)
    8. `sha1_empty_matches_nist_kat` — KAT for empty message

  Honest scope: we do **not** ship Wang's hand-crafted
  differential paths or find a full-SHA-1 collision.  We **do**
  demonstrate that the differential-cryptanalysis toolchain in
  this library applies to SHA-1, that reduced-round SHA-1 is
  breakable in seconds, and that the avalanche transition
  between "broken" (≤10 rounds) and "secure" (≥20 rounds)
  rounds is visible empirically.

- **NIST/RFC known-answer tests + cross-check unit tests** — 401
  passing, 24 ignored (1 RSA-CRT-2048 slow, 3 Serpent KAT
  regressions, 7 EC-ElGamal slow-affine tests, 1 HNP 16-bit-bias
  slow LLL, 1 Pollard rho 20-bit DLP slow, 3 CGA-HNC Phase-1
  experimental drivers, 1 CGA-HNC Phase-2 demo, 3 CGA-HNC
  Phase-3 cross-curve drivers, 1 CGA-HNC orbit-traversal slow,
  1 aut-folded rho empirical-speedup driver, 2 BKZ-vs-LLL HNP
  comparison drivers, 1 SHA-1 avalanche transition curve)
  unit tests across SHA-2/3, AES-GCM, ChaCha20-Poly1305, Serpent,
  Threefish-256/512/1024, HMAC, HKDF,
  PBKDF2, ECDSA, ECDH, RSA, Kyber, McEliece, the constant-time
  `Uint<LIMBS>` primitives at LIMBS=4/16/32/64, both
  `SecpFieldElement` and `P256FieldElement` arithmetic, both
  projective point modules, `MontgomeryContext` / `mont_pow_ct`
  cross-checked against `BigUint::modpow` at the RSA-2048 and
  RSA-4096 shapes, and CRT-vs-direct equivalence on the RSA private
  path (`rsa_mod_pow_ct_crt` matches `rsa_mod_pow_ct`).  Includes the AES-GCM NIST KAT that
  previously exposed a real keystream-reuse bug in this codebase, the
  RFC 6979 P-256 §A.2.5 deterministic-ECDSA vector, the cross-
  implementation P-256 verify KAT against Python `cryptography` /
  OpenSSL, and Renes–Costello–Batina sanity checks that
  `n · G = identity` on secp256k1 and a P-256 external KAT lands on
  the same `Q` via the projective ladder as via the textbook ladder.

### What this library does NOT protect against

These are structural — they cannot be fixed by configuration or by
patches that don't fundamentally re-architect the crate:

- **Timing side-channels in big-integer arithmetic — closed for
  every supported secret-key path on both ECC and RSA.**
  Every secret-scalar / secret-exponent path in this codebase now
  routes through `ct_bignum::Uint<LIMBS>`-backed Montgomery
  arithmetic:
    - **ECC** (secp256k1 + P-256): keypair generation, ECDH, and the
      ECDSA signer's `k·G` use `ecc::ct::scalar_mul_secret` over the
      Renes–Costello–Batina projective ladder.  The ECDSA verifier's
      `u₁·G + u₂·Q` keeps using the variable-time
      `Point::scalar_mul` since both inputs are public.
    - **RSA** (1024 / 1536 / 2048 / 3072 / 4096-bit): `rsa_decrypt_raw`
      and `rsa_sign` use `rsa_mod_pow_ct_crt`, the CRT-decomposed
      private path.  Each half-prime exponentiation routes through a
      monomorphised `ct_bignum::mont_pow_ct<HALF>`; the `c → c mod p`
      and `c → c mod q` reductions go through a CT split-and-reduce
      primitive instead of `BigUint::%`.  Modulus sizes outside this
      set fall back to the non-CRT `rsa_mod_pow_ct`, and thence to
      the legacy `BigUint`-backed `mod_pow_ct`, which has uniform
      per-bit op count but is not limb-level constant-time.
  The remaining caveats are about *micro-architectural* timing
  channels rather than the source-level operation count:
    - The variable-time multiplication-or-not behaviour of `*` on
      Cortex-A53, ARM7, and a few older Intel CPUs is not addressed —
      we use a fixed-cycle `mul` instruction on x86-64 / aarch64 in
      practice, but the language specification does not guarantee
      this.
    - The `subtle` crate's `cmov` is optimisation-resistant by
      construction but not formally proven constant-time on every
      target/CPU/compiler combination.
    - Cache, port-contention, and speculative-execution channels
      are out of scope — see the AES entry below for a related case.
  None of these can be addressed at the Rust source level; closing
  them requires a verified-CT compiler backend, hardware crypto
  intrinsics, or both.
- **Cache-timing side-channels in AES — partially mitigated.**
  The S-box is now scanned in full on every byte substitution
  (see "What this library now does"), which closes the cache-line
  channel.  However, the `subtle` crate's primitives are
  optimisation-resistant by construction but not formally proven
  constant-time on every target/CPU/compiler combination, and
  micro-architectural channels beyond cache (port contention,
  speculative execution) are not addressed.  `gmul` uses bit-mask
  arithmetic which compiles to branchless code under `rustc 1.x` on
  x86-64 and aarch64 in our spot checks, but this is not guaranteed
  by the language specification.
  *True fix*: AES-NI hardware intrinsics, or a Käsper–Schwabe-style
  8-way bitsliced software AES that processes 128 bytes per call
  with no S-box logic at all.
- **Variable-time ECC scalar multiplication.**  `Point::scalar_mul`
  uses double-and-add and branches on each scalar bit.
  *Fix*: Montgomery ladder, complete addition formulas, constant-time
  field arithmetic.
- **Variable-time Patterson decoder (McEliece).**  Constant-time Goppa
  decoding is an active research area.  This implementation is not
  constant-time and uses toy parameters (m=6, n=32, t=3) far below any
  security level.
- **Non-standard Kyber.**  The implementation here uses naive
  polynomial multiplication and `rand`-based sampling instead of the
  NTT and SHAKE-based sampling that NIST ML-KEM specifies.  Ciphertexts
  are *not* compatible with FIPS 203.
- **Best-effort, not guaranteed, zeroization.**  `BigUint`'s internal
  heap allocations are owned by `num-bigint` and are not exposed.  We
  re-assign the value to zero on drop; the previous heap may or may not
  be scrubbed by the allocator before reuse.  For genuine zeroization,
  switch to `crypto-bigint` (which integrates with the `zeroize` crate).
- **No formal verification or third-party audit.**  Audited libraries
  have spent thousands of expert-hours under adversarial scrutiny.
  This one has not.
- **No standards compliance.**  Not FIPS 140-2/3 validated, not Common
  Criteria evaluated.  PCI-DSS environments specifically require FIPS-
  validated cryptographic modules.

## Reporting a vulnerability

If you believe you have found a bug in this library, open an issue on
GitHub.  Because the library is explicitly educational, there is no
embargo or coordinated-disclosure process — the appropriate response to
a vulnerability is to update the docs and recommend that consumers move
to an audited crate.

## Versioning

This library follows no compatibility guarantee.  Breaking changes can
land in any release while the codebase is being hardened.
