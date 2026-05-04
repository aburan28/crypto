# Deferred items — what's next on the roadmap

The library now ships ~430 unit tests across the cryptanalysis,
ECC, symmetric, asymmetric, hash, KDF, and PQC families.  Several
high-value items remain deferred from the most recent gap-audit.
This note records what they are, *why* they were deferred, and
what each would take to ship responsibly.

---

## Tier 1 (medium effort, high impact)

### 1. **Ed25519 + X25519** — modern ECC

**Why deferred**: needs a Curve25519 field-arithmetic implementation
before Ed25519 sits on top.  The library's existing field machinery
is parametrised by a runtime `BigUint` modulus; Curve25519 wants a
fixed `2²⁵⁵ − 19` prime with constant-time field operations
specifically tuned for it (radix-25.5 representation, Mersenne-form
reduction).

**Estimated cost**:
- Curve25519 field (~600 LoC): radix-25.5 limbs, Mersenne-reduce,
  constant-time inversion via Fermat's little theorem
- Edwards-form points (~300 LoC): twisted-Edwards addition formulas,
  Montgomery-ladder for X25519
- Ed25519 sign/verify (~400 LoC): SHA-512-based deterministic nonce,
  the EdDSA encoding/decoding rules
- Tests: RFC 8032 test vectors

**Total**: ~1500 LoC, 1-2 days of focused work.

### 2. **Argon2** + **scrypt** — modern password hashing

**Why deferred**: Argon2 needs the **BlaMka** function (a custom
permutation built on Blake2b's G-function), parallel-lane processing,
and careful memory-block management.  scrypt needs **BlockMix** and
**ROMix** with the Salsa20/8 core.  Both are intricate; getting them
right requires careful KAT validation.

**Estimated cost**:
- Argon2 (~600 LoC): BlaMka + Blake2b core + parallel lanes + KAT validation
- scrypt (~400 LoC): Salsa20/8 + BlockMix + ROMix + KAT validation
- BLAKE2b dependency (~250 LoC): if not already present

**Total**: ~1250 LoC, 1-2 days.

### 3. **ML-KEM (Kyber)** — NIST-standardised PQC KEM

**Why deferred**: needs ring-LWE polynomial arithmetic over
`Z_q[x]/(x²⁵⁶ + 1)` with Number-Theoretic Transform (NTT) for fast
multiplication.  The NTT alone is ~400 LoC of careful modular code;
the rest of Kyber (CBD sampling, compress/decompress, encryption,
KEM transformation) sits on top.

**Estimated cost**:
- NTT + polynomial ring (~500 LoC)
- Kyber-512/768/1024 parameter sets + sampling + encrypt/decrypt
  (~800 LoC)
- Fujisaki-Okamoto transform for KEM (~200 LoC)
- Known-answer tests against the FIPS 203 KAT vectors

**Total**: ~2000 LoC, 2-3 days.

The library already has a "simplified Kyber" in `src/pqc/`;
upgrading it to actual FIPS 203 ML-KEM is the deferred work.

---

## Tier 2 (large effort, very high impact)

### 4. **ML-DSA (Dilithium)** — NIST-standardised PQC signature

Same NTT infrastructure as ML-KEM, plus rejection-sampling
signature scheme.  ~2500 LoC.

### 5. **SLH-DSA (SPHINCS+)** — NIST-standardised stateless hash signature

Pure hash-based; no number theory.  Needs WOTS+, FORS, and a
Merkle-tree harness.  ~2000 LoC.  Easier than Dilithium but with
significantly larger signatures (~17 KB).

### 6. **FN-DSA (Falcon)** — NIST-standardised PQC signature (lattice)

Hardest of the four NIST PQC signatures: needs floating-point
fast Fourier transform, Gaussian sampling, and tower-of-fields
arithmetic.  ~4000 LoC.  Notoriously hard to implement correctly;
several early implementations had subtle floating-point bugs.

---

## Tier 3 (multi-week, separate project territory)

### 7. **zk-SNARKs (Groth16)**

The hardest item by far.  Requires:

- **BLS12-381** pairing-friendly curve (~3000 LoC):
  - `F_p` for the 381-bit base field
  - `F_{p²}` quadratic extension
  - `F_{p⁶}` cubic extension over `F_{p²}`
  - `F_{p¹²}` quadratic extension over `F_{p⁶}`
  - G1 (over `F_p`) and G2 (over `F_{p²}`) point arithmetic
  - Optimal Ate pairing (Miller loop + final exponentiation)
  - Frobenius endomorphism + GLV/GLS optimisations
- **R1CS constraint system** (~500 LoC): rank-1 quadratic constraints
  representation, witness assignment
- **QAP transformation** (~400 LoC): R1CS → quadratic arithmetic
  program via Lagrange interpolation
- **Trusted setup** (~300 LoC): structured reference string generation;
  toxic-waste handling
- **Groth16 prover** (~800 LoC): MSM-heavy; the bulk of cost
- **Groth16 verifier** (~400 LoC): three pairings
- **Tests against published circuits**: e.g. ZCash Sapling verification

**Total**: ~5500 LoC minimum for a *toy* but functional Groth16.
Production-grade adds maybe another 3000 LoC for batched MSM,
multi-exp tricks, and proper serialisation.

This is **multi-week** focused work, not session-scale.  Would
need its own repo / sub-crate.

### 8. **Bulletproofs** — range proofs without trusted setup

Needs Pedersen commitments + inner-product argument over Ristretto255
or BLS12-381 G1.  ~3000 LoC.  Less infrastructure than Groth16
because no pairings, but the inner-product argument is intricate.

### 9. **PLONK / Marlin / Halo2**

Universal-setup SNARK schemes.  Each is its own multi-thousand-LoC
project.  PLONK is conceptually cleaner than Groth16 but the
permutation-argument machinery is large.

### 10. **STARKs**

No trusted setup, post-quantum-secure.  Needs FRI commitment,
arithmetic intermediate representation (AIR), Reed-Solomon codes.
~5000 LoC.  Intricate but cleaner than SNARKs in many ways.

---

## Tier 4 (small but valuable)

### 11. **MD5 + cryptanalysis**

MD5 is 5x simpler than SHA-1 and the differential cryptanalysis
that breaks it is even cleaner pedagogically.  Wang et al.'s 2004
collision attack runs in seconds on a laptop today.  ~250 LoC for
MD5 + ~400 LoC for the simplified differential collision finder.

### 12. **DES / 3DES + linear cryptanalysis**

Matsui 1993, the original linear cryptanalysis paper, attacks DES
specifically.  Pedagogically perfect — would let the library ship
*both* differential cryptanalysis (already done for SHA-1) and
linear cryptanalysis.  ~400 LoC for DES + ~600 LoC for the linear
attack.

### 13. **CMAC + GMAC + Poly1305 standalone**

GMAC and Poly1305 already exist inside the AEAD modes; would just
need extraction into standalone interfaces.  CMAC needs ~150 LoC.

### 14. **AES-OCB / AES-SIV / AES-EAX**

Additional AEAD modes.  Each ~300 LoC.

### 15. **HMAC-DRBG / CTR-DRBG / Hash-DRBG** — NIST SP 800-90A

The library already has a special-cased HMAC-DRBG inside ECDSA
(for RFC 6979); it should be promoted to a standalone module.
~200 LoC.

### 16. **SipHash**

Keyed PRF used by Rust's `HashMap` for DoS resistance.  ~150 LoC.

### 17. **SHA-512 family**

Wider-state SHA-2.  Already conceptually present (the SHA-256
implementation is parametrised).  ~250 LoC.

---

## What this session shipped

| Item | Tests | LoC |
|---|---|---|
| RIPEMD-160 | 3 (incl. million-a + Bitcoin pubkey) | ~250 |
| Bitcoin HASH160 | (above) | (above) |
| DER codec for ECDSA | 6 | ~250 |
| AES-CBC + PKCS#7 | 4 | ~150 |
| Vaudenay 2002 padding-oracle attack | 2 (incl. multi-block recovery) | ~250 |
| Schnorr BIP-340 | 6 (incl. both official BIP test vectors) | ~400 |
| **Total** | **22** | **~1300** |

Net library state: **423 passing, 0 failed, 24 ignored**, +5 new
modules (`hash::ripemd160`, `encoding::der`, `symmetric::aes_cbc`,
`ecc::schnorr`, `encoding/`).

---

## Strategic recommendation

The next single highest-leverage addition would be **Ed25519 +
X25519** (Tier 1 #1) — it fills the most embarrassing gap (every
modern cryptographic deployment uses Curve25519) and pairs naturally
with the existing constant-time ECC infrastructure.  After that,
**Argon2** (Tier 1 #2) closes the password-hashing gap.

The PQC items (#3-#6) are valuable but would benefit from being
done together (shared NTT infrastructure for ML-KEM and ML-DSA;
shared hash-tree machinery for SLH-DSA).  Best as a focused multi-day
session.

zk-SNARKs (Tier 3) are a separate project: not "deferred" but
"out of scope for a single library."  If the user wants
verifiable computation, the right answer is to consume the
existing audited ZK ecosystems (`arkworks`, `bellman`, `halo2`)
rather than re-implement from scratch.

The Tier 4 items (#11-#17) are small enough to ship in a single
session each; they'd round out the educational story.
