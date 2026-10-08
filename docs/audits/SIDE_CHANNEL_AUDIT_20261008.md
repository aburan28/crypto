# Side-channel (constant-time) inventory of the `crypto` crate — 2026-10-08

Source-level review of every code path in `src/` that handles secret material, classifying each as
constant-time (CT) or variable-time (VT) with file and line references, plus the gaps in tooling and
documentation. A companion scope note for the `cryptanalysis` repository lives there as
`docs/SIDE_CHANNEL_SCOPE.md`. Baseline: `origin/main` at `4135dc700`.

This is a reading of the source, not a measurement. No timing or cache-trace experiments were run, and
the compiler's freedom to reintroduce branches from `subtle` selects was not checked (there is no
`black_box`, `write_volatile` or `#[inline(never)]` barrier anywhere in `src/`).

## 1. Scope: where secrets are

Only the library crate (`src/`) handles long-lived or per-message secrets: ECC keys and nonces, RSA,
ElGamal/Paillier, BLS, the ZK provers, symmetric keys, MAC/KDF keys and the PQC schemes. Everything else
in the repository (`src/cryptanalysis/*`, `src/bin/ic`, `examples/`, `ecc2k130/`, `gpu/*`, `research/*`,
`sage/`, `scripts/`) is search or benchmark code on public inputs with planted, seeded known answers;
constant time is irrelevant there and nothing below applies to it.

## 2. What is constant-time at source level

| Component | Evidence | Caveats |
|---|---|---|
| Fixed-width integers `ct_bignum::Uint` | `src/ct_bignum.rs:190-210` cmov add/sub; `275-374` Montgomery mul/sqr (branches only on loop indices); `525-556` `mont_pow_ct` runs all bits with cmov multiply | `from_biguint` (`:77`) stops at the byte length of its input |
| secp256k1 / P-256 field | `src/ecc/secp256k1_field.rs`, `src/ecc/p256_field.rs`: Montgomery form; inversion by Fermat ladder over a public exponent (`secp :150-172`, `P-256 :109-125`); `PartialEq` via `ct_eq` | `from_biguint` reduces with `BigUint %` (`secp :86-88`, `P-256 :67-69`) |
| secp256k1 / P-256 points | Renes–Costello–Batina complete formulas (`secp256k1_point.rs:132,181`; `p256_point.rs:126,181`); ladders `scalar_mul_ct` always run `order_bits` steps with cmov (`secp :258-274`, `P-256 :236-246`) | |
| Secret scalar multiplication dispatch | `src/ecc/ct.rs:39-51` `scalar_mul_secret`; `:79-105` Coron-blinded variant (`k + r·n`, 64-bit `r` from `OsRng`) | Doc comment (`:33-38`) says unknown curves panic; code falls through to the VT affine ladder (`:48`), see §3 |
| ECDSA nonce point, key generation, ECDH | `src/ecc/ecdsa.rs:72`; `src/ecc/keys.rs:55-80`; `src/ecc/ecdh.rs:38,60` all route through `scalar_mul_secret` | CT only on secp256k1 and P-256 |
| RFC 6979 nonce | `src/ecc/ecdsa.rs:104+` (HMAC-SHA-256 DRBG) | |
| RSA private operation | `src/asymmetric/rsa.rs:108-117` CRT over `mont_pow_ct` for 1024–4096-bit keys; CT PKCS#1 unpadding (`:618+`) | Other sizes fall back to BigUint `mod_pow_ct`; `rsa_decrypt` returns an `Err` on padding failure (`:671`), an oracle in the API shape |
| AES S-box, AEAD tags | `src/symmetric/aes.rs:98-115` full-table scan with `ct_eq`/`conditional_select`; tag compares with `ct_eq` in gcm, aes, chacha20, salsa20, aegis, ascon; XOR-accumulate-then-branch in ccm/eax/ocb3/siv/gcm_siv/kw | |
| POLYVAL (GCM-SIV) | `src/symmetric/modes/gcm_siv.rs:76-132` masked | |
| HMAC / HKDF / HMAC-DRBG verify | `src/hash/hmac.rs:63-83`, `src/kdf/hkdf.rs:25-38` `ct_eq` | |
| RNG | `src/utils/random.rs:25-58` `OsRng` with rejection sampling | |

## 3. Variable-time paths that touch secrets (ranked)

1. **Scalar arithmetic mod n in every signature scheme is `num-bigint`.** ECDSA `r·d`, `z + r·d`,
   `k⁻¹·(…)` (`src/ecc/ecdsa.rs:85-96`); Schnorr `s = k + e·d` (`src/ecc/schnorr.rs:186`) and the
   `y_is_odd` negations (`:143, :172`) with `neg_mod` branching on zero (`:87-91`). The `mod_pow_ct`
   ladder behind `mod_inverse_prime_ct` is itself a plain `if bit` over BigUint (`src/utils/mod.rs:72`).
   The CT point ladders are therefore wrapped in VT scalar arithmetic.
2. **Extended-Euclid inversion of a private key.** SM2 `mod_inverse(1 + d)` (`src/ecc/sm2.rs:151`);
   EC-KCDSA `mod_inverse(x)` (`src/ecc/ec_kcdsa.rs:85`). `mod_inverse` (`src/utils/mod.rs:108-135`) has
   data-dependent loop count and division.
3. **The generic "CT" ladder is not CT.** `Point::scalar_mul_ct` (`src/ecc/point.rs:128-152`) is an
   affine BigUint ladder with `if bit {…} else {…}` (`:143`) that starts from `Infinity`, so leading zero
   bits take a different `match` path. It is the fallback for every curve other than secp256k1/P-256
   (`src/ecc/ct.rs:48`) and is used directly by SM2 (`sm2.rs:138, 298`), GOST (`gost_3410_2012.rs:107,
   127, 196`), EC-KCDSA (`ec_kcdsa.rs:90, 127`) and EC-ElGamal (`ec_elgamal.rs:99-163, 371-372`).
4. **Ed25519 / Ed448 use plain double-and-add on the secret scalar and nonce.** `EdPoint::scalar_mul`
   (`src/ecc/ed25519.rs:187-197`, `src/ecc/ed448.rs:234-245`) branches on `scalar.bit(i)` and loops over
   `scalar.bits()`; called from signing at `ed25519.rs:280-299`, `ed448.rs:370-391`.
5. **X25519 / X448 ladder swap is a branch.** `if s == 1 { mem::swap … }` (`src/ecc/x25519.rs:92-98`,
   `src/ecc/x448.rs:118-150`), over BigUint field arithmetic with `modpow` inversion (`curve25519.rs:68-74`).
6. **VT `scalar_mul` / `mod_pow` on secrets elsewhere.** EC-ElGamal secret plaintext `m`
   (`ec_elgamal.rs:139, 343`); finite-field ElGamal keygen/encrypt/decrypt (`elgamal.rs:145, 178-179,
   194-195`); Paillier decrypt `mod_pow(c, λ)` (`paillier.rs:219`); BLS `sk·G1`, `sk·H(m)`
   (`bls12_381/signature.rs:85, 94`); Pedersen, Schnorr-ZKP, Chaum–Pedersen, Bulletproofs blinding and the
   KZG setup `tau` (`src/zk/*`, `kzg.rs:96-99`); CSIDH loops `|secret[i]|` times (`pqc/csidh.rs:323-330`).
7. **Key-dependent branches in symmetric modes.** GHASH `gcm_mult` branches on bits of the hash key
   (`src/symmetric/modes/gcm.rs:45, 55`; also `aes.rs:427, 439`); doubling steps branch on the top bit of
   key-derived values in `ocb3.rs:71`, `siv.rs:47`, `pmac.rs:56, 69, 79`.
8. **Binary-field and BLS field arithmetic.** `f2m.rs` multiplication branches on operand bits
   (`:265, :348`), so its Itoh–Tsujii inversion (`:389`) inherits that; `bls12_381/fq.rs` uses
   extended Euclid (`:74-76`) and `modpow` (`:79-81`).
9. **Secret hygiene.** The `zeroize` crate is declared (`Cargo.toml:34`) and never imported in `src/`;
   `Drop` impls call `BigUint::set_zero()` (`src/ecc/keys.rs:28-30`, `rsa.rs:400`, `ec_elgamal.rs:73`,
   `elgamal.rs:120`, `paillier.rs:103`), which resets the limb vector's length without overwriting
   memory. `EccPrivateKey` derives `Clone` and `Debug` (`keys.rs:22`).

## 4. Documentation drift

* `README.md:15-21` says all big-integer arithmetic is `num-bigint` and that AES uses a table S-box.
  Both are stale for secp256k1/P-256 and AES; both remain true for every scheme in §3.
* `SECURITY.md` ("What this library now does") states the ECC secret-scalar paths are closed for every
  supported secret-key path. That holds for key generation, ECDH, ECDSA and Schnorr on secp256k1 and
  P-256 only; the schemes in §3 items 2–6 are not covered and should be listed under "does NOT protect".
* `src/ecc/ct.rs:33-38` documents a panic that does not happen (fallback instead).

## 5. Tooling gap

There are no side-channel tests: no dudect-style timing distributions, no Valgrind/ctgrind secret
marking, no TVLA. CT code is validated for correctness against `BigUint` only. Rust has no
`black_box`/volatile barriers protecting the `subtle` selects from being compiled back into branches.

## 6. Recommended remediation order

1. Move all scalar arithmetic mod n in ECDSA/Schnorr onto `ct_bignum::Uint` (add, mul, Fermat inverse
   with a public exponent), and replace `mod_pow_ct`'s `if bit` with a `subtle` select. This closes the
   only remaining VT step in the two schemes already claimed CT.
2. Make `scalar_mul_secret` do what its doc says (reject unsupported curves) or give the affine ladder a
   real cswap; remove `mod_inverse` from SM2/EC-KCDSA in favour of a Fermat inverse.
3. Port Ed25519/X25519 to a fixed-width field (`Uint<4>`) with a masked cswap; until then, label them VT.
4. Replace the GHASH/OCB/SIV/PMAC branches with masks (`gcm_siv.rs` already has the pattern).
5. Use `zeroize` for real (on `Uint`, byte buffers; `BigUint` cannot be zeroized reliably), drop
   `Clone`/`Debug` on key types or route `Debug` through a redacting impl.
6. Add a timing harness (`dudect-bencher` or a two-class Welch t-test over `Instant` samples) for the CT
   paths in §2 and a CI job that runs it; add `core::hint::black_box` at the cmov sites.
7. Correct `README.md` and `SECURITY.md` to match §2/§3.
