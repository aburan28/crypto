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
- **NIST/RFC known-answer tests** — 119 unit tests across SHA-2/3,
  AES-GCM, ChaCha20-Poly1305, HMAC, HKDF, PBKDF2, ECDSA, ECDH, RSA,
  Kyber, and McEliece, including the AES-GCM NIST KAT that previously
  exposed a real keystream-reuse bug in this codebase.

### What this library does NOT protect against

These are structural — they cannot be fixed by configuration or by
patches that don't fundamentally re-architect the crate:

- **Timing side-channels in big-integer arithmetic.**  All RSA, ECDSA,
  ECDH, and ECC field operations route through `num-bigint`, which is
  not constant-time.  Modular exponentiation, modular inverse, and
  comparison all branch on the value of secret operands.  A network or
  co-resident attacker can recover keys from timing data.
  *Fix*: replace `num-bigint` with `crypto-bigint` (constant-time) and
  rewrite the affected modules.  Months of work.
- **Cache-timing side-channels in AES.**  The S-box is a 256-byte
  lookup table.  Bernstein (2005) and Osvik–Shamir–Tromer (2006) showed
  this is exploitable from a co-resident process.
  *Fix*: AES-NI hardware intrinsics, or a bitsliced software AES.
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
