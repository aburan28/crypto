# crypto — Comprehensive Educational Cryptography Library

A from-scratch implementation of the major cryptographic algorithm families in Rust.
Every algorithm is documented with the underlying mathematics, making this a study
companion as much as a working library.

> ## ⚠️ DO NOT USE IN PRODUCTION
>
> This library is **not safe** for any real use case where confidentiality,
> integrity, or authenticity matters.  It is an educational reimplementation
> intended for studying the algorithms.
>
> **What is wrong with it:**
> - **Not constant-time.** All big-integer arithmetic uses `num-bigint`,
>   which branches on secret values.  Every RSA, ECDSA, ECDH, and ECC
>   operation leaks via timing.  AES uses a lookup-table S-box vulnerable
>   to cache-timing attacks.
> - **No third-party audit.**  We've already found one real keystream-reuse
>   bug in AES-GCM during testing — there could easily be more.
> - **Non-standard Kyber.**  The implementation here is *not* compatible
>   with NIST ML-KEM (FIPS 203).
> - **Toy McEliece parameters** (m=6, n=32, t=3) — far below any security
>   level.
>
> See [`SECURITY.md`](./SECURITY.md) for the full list of structural
> limitations and recommended audited alternatives (`aws-lc-rs`, `ring`,
> RustCrypto, dalek, etc.).
>
> A hardening pass has been applied (OS-backed RNG, RFC 6979 deterministic
> ECDSA, constant-time MAC verification, best-effort key zeroization,
> textbook-RSA hidden from the public API).  This makes the codebase a
> better educational artifact.  **It does not make it production-safe.**

---

## Algorithms implemented

| Module | Algorithm | Standard | From scratch? |
|--------|-----------|----------|---------------|
| `ecc` | Finite field arithmetic | — | ✓ |
| `ecc` | Short-Weierstrass point add/double | SEC 1 | ✓ |
| `ecc` | Scalar multiplication (double-and-add) | SEC 1 | ✓ |
| `ecc` | ECDSA sign + verify | FIPS 186-4 | ✓ |
| `ecc` | ECDH key exchange | RFC 8422 | ✓ |
| `ecc` | secp256k1 (Bitcoin/Ethereum curve) | SEC 2 | — |
| `ecc` | P-256 (NIST / TLS curve) | FIPS 186-4 | — |
| `symmetric` | AES-128 / AES-256 block cipher | FIPS 197 | ✓ |
| `symmetric` | AES-CTR streaming mode | NIST SP 800-38A | ✓ |
| `symmetric` | AES-GCM authenticated encryption | NIST SP 800-38D | ✓ |
| `symmetric` | ChaCha20 stream cipher | RFC 8439 | ✓ |
| `symmetric` | Poly1305 MAC | RFC 8439 | ✓ |
| `symmetric` | ChaCha20-Poly1305 AEAD | RFC 8439 | ✓ |
| `asymmetric` | Miller-Rabin primality test | — | ✓ |
| `asymmetric` | RSA key generation | PKCS#1 | ✓ |
| `asymmetric` | RSA encrypt/decrypt (PKCS#1 v1.5) | RFC 8017 | ✓ |
| `asymmetric` | RSA sign/verify | PKCS#1 | ✓ |
| `hash` | SHA-256 / SHA-224 | FIPS 180-4 | ✓ |
| `hash` | SHA3-256 / SHA3-512 (Keccak) | FIPS 202 | ✓ |
| `hash` | SHAKE128 / SHAKE256 (XOF) | FIPS 202 | ✓ |
| `hash` | BLAKE3 | BLAKE3 spec | via crate |
| `kdf` | HMAC-SHA256 | RFC 2104 | ✓ |
| `kdf` | HKDF-SHA256 | RFC 5869 | ✓ |
| `kdf` | PBKDF2-HMAC-SHA256 | RFC 8018 | ✓ |
| `pqc` | Kyber / ML-KEM (simplified educational variant) | FIPS 203 | ✓ |

---

## Building and running

```bash
# Build
cargo build

# Run the full demo
cargo run -- demo

# Individual commands
cargo run -- hash sha256 "hello world"
cargo run -- ecc keygen
cargo run -- ecdsa sign "my message"
cargo run -- ecdh
cargo run -- aes
cargo run -- chacha
cargo run -- hkdf
cargo run -- pqc
cargo run -- rsa demo        # Note: slow — generates a 1024-bit RSA key

# Run all tests
cargo test
```

---

## Mathematics background

### Elliptic curves (ECC)

An elliptic curve over a prime field **F_p** is the set of points satisfying:

```
y² ≡ x³ + ax + b  (mod p)
```

together with a special "point at infinity" **O** that acts as the additive identity.

**Point addition** (P ≠ Q):
```
λ = (y₂ − y₁) / (x₂ − x₁)  mod p
x₃ = λ² − x₁ − x₂           mod p
y₃ = λ(x₁ − x₃) − y₁        mod p
```

**Point doubling** (P = Q):
```
λ = (3x₁² + a) / (2y₁)  mod p
```

**Scalar multiplication** k·P means adding P to itself k times, computed
efficiently with the double-and-add algorithm in O(log k) steps.

#### secp256k1 (Bitcoin's curve)
- a = 0, b = 7
- p ≈ 2²⁵⁶ (a 256-bit prime)
- Group order n ≈ 2²⁵⁶ (also prime)
- Private key: a random scalar d ∈ [1, n-1]
- Public key: Q = d·G

### ECDSA

Sign(d, z):  pick random k; r = x(k·G) mod n; s = k⁻¹(z + rd) mod n  
Verify(Q, z, r, s):  w = s⁻¹; check x(zw·G + rw·Q) mod n = r

### AES

AES operates on a 4×4 byte state through `Nr` rounds (10/14 for 128/256-bit keys).
Each round applies:
1. **SubBytes** — non-linear S-box substitution (inverse in GF(2⁸) + affine map)
2. **ShiftRows** — cyclic rotation of rows 0..3 by 0..3 positions
3. **MixColumns** — matrix multiplication over GF(2⁸) (skipped in the last round)
4. **AddRoundKey** — XOR with the current round key

**GCM** adds GHASH (polynomial multiplication over GF(2¹²⁸)) for authentication.

### ChaCha20-Poly1305

ChaCha20 is an ARX (Add-Rotate-XOR) stream cipher with a 32-byte key and 12-byte nonce.
The quarter-round operation mixes four 32-bit words in a 4×4 state matrix.
Poly1305 provides authentication over GF(2¹³⁰−5).

### RSA

1. Generate primes p, q using Miller-Rabin (probabilistic) primality testing.
2. n = p·q; λ(n) = lcm(p-1, q-1); e = 65537; d = e⁻¹ mod λ(n).
3. Encrypt: c = mᵉ mod n. Decrypt: m = cᵈ mod n.
4. Sign: s = H(m)ᵈ mod n. Verify: H(m) = sᵉ mod n.

### HKDF

Extract:  PRK = HMAC-SHA256(salt, IKM)  
Expand:   OKM = T(1) ‖ T(2) ‖ … where T(i) = HMAC-SHA256(PRK, T(i-1) ‖ info ‖ i)

### Kyber / ML-KEM (simplified)

Kyber is a Module-LWE-based KEM operating in the polynomial ring R_q = Z_q[x]/(xⁿ+1).
- Key gen: sample secret **s** and error **e**; public key t = A·**s** + **e**
- Encapsulate: sample r; ciphertext (u = Aᵀr + e₁, v = t·r + e₂ + μ)
- Decapsulate: recover μ = v − s·u, hash to shared secret

The simplified implementation here uses schoolbook polynomial multiplication
(O(n²)) instead of the NTT. Ciphertexts are not compatible with NIST FIPS 203.

---

## Project structure

```
src/
├── lib.rs              — library root
├── main.rs             — CLI demo (clap)
├── ecc/
│   ├── field.rs        — FieldElement over F_p
│   ├── point.rs        — Point (add, double, scalar_mul)
│   ├── curve.rs        — CurveParams (secp256k1, P-256)
│   ├── keys.rs         — EccKeyPair generation
│   ├── ecdsa.rs        — ECDSA sign + verify
│   └── ecdh.rs         — ECDH shared secret
├── symmetric/
│   ├── aes.rs          — AES-128/256, CTR, GCM
│   └── chacha20.rs     — ChaCha20 + Poly1305 + AEAD
├── asymmetric/
│   └── rsa.rs          — Miller-Rabin, RSA keygen, encrypt, sign
├── hash/
│   ├── sha256.rs       — SHA-256 / SHA-224 from scratch
│   ├── sha3.rs         — Keccak / SHA-3 / SHAKE from scratch
│   └── blake3.rs       — BLAKE3 via the blake3 crate
├── kdf/
│   ├── hkdf.rs         — HMAC-SHA256, HKDF extract + expand
│   └── pbkdf2.rs       — PBKDF2-HMAC-SHA256
├── pqc/
│   └── kyber.rs        — Simplified Kyber KEM
├── keys/
│   └── mod.rs          — Key management helpers
└── utils/
    ├── random.rs       — CSPRNG, random_scalar
    ├── encoding.rs     — hex, BigUint ↔ bytes
    └── mod.rs          — mod_pow, mod_inverse (ext. Euclidean)
```

---

## References

- FIPS 197 — AES: https://nvlpubs.nist.gov/nistpubs/FIPS/NIST.FIPS.197.pdf
- FIPS 180-4 — SHA: https://nvlpubs.nist.gov/nistpubs/FIPS/NIST.FIPS.180-4.pdf
- FIPS 202 — SHA-3 / Keccak: https://nvlpubs.nist.gov/nistpubs/FIPS/NIST.FIPS.202.pdf
- FIPS 203 — ML-KEM (Kyber): https://nvlpubs.nist.gov/nistpubs/FIPS/NIST.FIPS.203.pdf
- FIPS 186-4 — ECDSA: https://nvlpubs.nist.gov/nistpubs/FIPS/NIST.FIPS.186-4.pdf
- RFC 8439 — ChaCha20-Poly1305: https://www.rfc-editor.org/rfc/rfc8439
- RFC 5869 — HKDF: https://www.rfc-editor.org/rfc/rfc5869
- RFC 8018 — PBKDF2: https://www.rfc-editor.org/rfc/rfc8018
- SEC 2 — secp256k1: https://www.secg.org/sec2-v2.pdf
- BLAKE3 spec: https://github.com/BLAKE3-team/BLAKE3-specs
