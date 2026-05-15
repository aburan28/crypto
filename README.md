# crypto — Comprehensive Educational Cryptography Library

A from-scratch implementation of the major cryptographic algorithm families in Rust,
plus a substantial **cryptanalysis suite** for studying how those algorithms break.
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
> - **No third-party audit.**
> - **Non-standard Kyber.**  The implementation here is *not* compatible
>   with NIST ML-KEM (FIPS 203).
> - **Toy McEliece parameters** (m=6, n=32, t=3) — far below any security
>   level.
>
> See [`SECURITY.md`](./SECURITY.md) for the full list of structural
> limitations and recommended audited alternatives (`aws-lc-rs`, `ring`,
> RustCrypto, dalek, etc.).

---

## Quick tour

```bash
# Build
cargo build --release

# Run the protocol demos
cargo run --release -- demo

# Hash a string
cargo run --release -- hash sha256 "hello world"

# Auto-run every applicable cryptanalysis attack against a named cipher
cargo run --release -- cryptanalysis auto --cipher toyspn-2r
cargo run --release -- cryptanalysis auto --cipher aes-2r

# Discover what attacks/ciphers are registered
cargo run --release -- cryptanalysis list-ciphers

# Run the full empirical research bench
cargo run --release -- cryptanalysis bench

# Tests (all 843 of them)
cargo test --release
```

---

## Algorithm coverage

### Symmetric primitives

| Module                    | Algorithm                            | Standard / Origin           |
|---------------------------|--------------------------------------|-----------------------------|
| `symmetric::aes`          | AES-128 / AES-256 + CTR + GCM        | FIPS 197 / SP 800-38A,D     |
| `symmetric::aes_cbc`      | AES-CBC                              | NIST SP 800-38A             |
| `symmetric::chacha20`     | ChaCha20 + Poly1305 + AEAD           | RFC 8439                    |
| `symmetric::serpent`      | Serpent block cipher                 | Anderson–Biham–Knudsen 1998 |
| `symmetric::threefish`    | Threefish (Skein tweakable cipher)   | Skein NIST submission        |
| `symmetric::sm4`          | SM4 Chinese national block cipher    | GB/T 32907                  |
| `symmetric::gost_magma`   | GOST R 34.12-2015 Magma (64-bit)     | GOST R 34.12-2015           |
| `symmetric::kuznyechik`   | GOST R 34.12-2015 Kuznyechik (128-bit) | GOST R 34.12-2015           |
| `symmetric::cmac`         | CMAC-AES MAC                          | NIST SP 800-38B             |
| `symmetric::des`          | DES (Data Encryption Standard, broken — teaching) | FIPS PUB 46-3 |
| `symmetric::des3`         | Triple-DES (EDE3, deprecated by NIST 2017)   | NIST SP 800-67  |
| `symmetric::rc4`          | RC4 stream cipher (broken — teaching) | Rivest 1987 / RFC 6229 vectors |
| `symmetric::rc5`          | RC5-32/12/16 block cipher             | Rivest, FSE 1994            |

### Block-cipher modes (generic over `BlockCipher<N>` trait)

The `symmetric::modes` submodule hosts modes that work against any
[`BlockCipher<N>`] implementor — AES, SM4, Serpent, Kuznyechik, Magma.

| Mode                            | Type                           | Standard               |
|---------------------------------|--------------------------------|------------------------|
| `modes::ecb` (PKCS#7 padded)    | Confidentiality (broken — teaching only) | NIST SP 800-38A §6.1 |
| `modes::cfb`                    | Self-synchronising stream      | NIST SP 800-38A §6.3   |
| `modes::ofb`                    | Pre-computable keystream       | NIST SP 800-38A §6.4   |
| `modes::ccm`                    | **AEAD** (Counter + CBC-MAC)   | NIST SP 800-38C / RFC 3610 |
| `modes::siv`                    | **AEAD, misuse-resistant** (deterministic) | RFC 5297            |
| `modes::eax`                    | **AEAD** (three-CMAC construction)         | Bellare-Rogaway-Wagner FSE 2004 |
| `modes::ocb3`                   | **AEAD** (single-pass, offset-codebook)    | RFC 7253 |
| `modes::gcm_siv` (+ POLYVAL)    | **AEAD, nonce-misuse-resistant** (GCM variant) | RFC 8452 |
| `modes::kw`                     | Deterministic key wrap         | NIST SP 800-38F / RFC 3394 |
| `modes::xts`                    | Disk encryption (ciphertext stealing) | IEEE P1619 / NIST SP 800-38E |

### Hashes / XOFs

| Module           | Algorithm                  | Standard                |
|------------------|----------------------------|-------------------------|
| `hash::md4`      | MD4 (broken — teaching/target)  | RFC 1320            |
| `hash::md5`      | MD5 (broken — teaching/target)  | RFC 1321            |
| `hash::sha256`   | SHA-224 / SHA-256          | FIPS 180-4              |
| `hash::sha512`   | SHA-384 / SHA-512          | FIPS 180-4              |
| `hash::sha3`     | SHA3-256/512, SHAKE128/256 | FIPS 202 (Keccak)       |
| `hash::blake3`   | BLAKE3                     | BLAKE3 spec (via crate) |
| `hash::ripemd160`| RIPEMD-160                 | ISO/IEC 10118-3         |
| `hash::siphash`  | SipHash-2-4                | Aumasson–Bernstein 2012 |
| `hash::sm3`      | SM3 Chinese national hash  | GB/T 32905              |
| `hash::streebog` | GOST R 34.11-2012 (256/512)| GOST R 34.11-2012       |

### Hash-function cryptanalysis

Generic attacks over the `HashFunction` / `MerkleDamgardHash` trait:

| Attack                                  | Module                     | Notes                                |
|----------------------------------------|----------------------------|--------------------------------------|
| **Length-extension attack**            | `cryptanalysis::hash_attacks::length_extension_attack` | Works against any MD-Damgård hash used as `H(secret \|\| msg)` MAC.  Tested on MD4, MD5, SHA-1. |
| **Birthday-collision search**          | `cryptanalysis::hash_attacks::birthday_collision_search` | Generic `√(2^b)` collision finder over a truncated hash. |
| **Joux multicollision** (CRYPTO 2004)  | `cryptanalysis::hash_attacks::joux_multicollision` | Chain `t` collisions → `2^t` equivalent messages at cost `t · 2^(b/2)`. |
| **Differential bias**                  | `cryptanalysis::hash_attacks::differential_bias` | Per-output-byte bias measurement on a chosen input differential. |
| **MD5 Wang differential**              | `cryptanalysis::md5_differential` | Reduced-round near-collision finder. |
| **SHA-1 differential**                 | `cryptanalysis::sha1_differential` | Round-function table + near-collision search. |

CLI:

```bash
# Auto-run every applicable hash attack
crypto cryptanalysis hash-auto --hash md5
crypto cryptanalysis hash-auto --hash md4
crypto cryptanalysis hash-auto --hash sha1

# Length-extension demo (interactive)
crypto cryptanalysis length-extension --hash sha1 \
    --secret-len 16 --message-hex "&role=guest" \
    --suffix-hex "&role=admin" --digest-hex <hex>

# AES related-key attack: key-schedule diffusion + avalanche + local collision
crypto cryptanalysis aes-related-key --key-bits 256 --rounds 4
crypto cryptanalysis aes-related-key --key-bits 128 --rounds 4 --avalanche-trials 1024

# AES visual demos — ASCII art for every major attack
crypto cryptanalysis aes-visual-demo --demo all              # all four below
crypto cryptanalysis aes-visual-demo --demo truncated-diff   # 4-round trail grid
crypto cryptanalysis aes-visual-demo --demo integral         # 3-round Square balance
crypto cryptanalysis aes-visual-demo --demo dfa              # Piret-Quisquater fault diff
crypto cryptanalysis aes-visual-demo --demo key-schedule     # AES-128 vs AES-256 diffusion

# Visual demos for every NON-AES attack module
crypto cryptanalysis visual-all --target all                 # full bundle (Markdown)
crypto cryptanalysis visual-all --target sbox-ddt            # DDT heat-maps (AES + Serpent + PRESENT)
crypto cryptanalysis visual-all --target sbox-lat            # LAT heat-maps
crypto cryptanalysis visual-all --target walsh               # Walsh-Hadamard signed spectrum bars
crypto cryptanalysis visual-all --target pollard-rho         # ρ-iteration trajectory plot
crypto cryptanalysis visual-all --target hnp                 # HNP recovery probability curve
crypto cryptanalysis visual-all --target bleichenbacher      # nonce-bias histogram
crypto cryptanalysis visual-all --target length-extension    # MD-Damgård chain diagram
crypto cryptanalysis visual-all --target joux                # multicollision tree
crypto cryptanalysis visual-all --target j0-twists           # twist factor-size bars
crypto cryptanalysis visual-all --target birthday            # cumulative collision-probability curve
```

### Library-wide visual demos

Beyond the cryptanalysis-side visualizations, every major module
ships ASCII-art demos:

```bash
# Symmetric ciphers
crypto visual symmetric --target ecb-penguin   # leak structure through ECB
crypto visual symmetric --target modes         # CBC/CFB/OFB/CTR/GCM/CCM/SIV/XTS/KW dataflow
crypto visual symmetric --target aes-trace     # 11 AES-128 state grids round-by-round
crypto visual symmetric --target chacha        # ChaCha20 4×4 word state evolution
crypto visual symmetric --target avalanche     # AES avalanche heat-map at 1/2/3/5/10 rounds

# Hash functions
crypto visual hash --target md                 # Merkle-Damgård chain (MD5/SHA-2/RIPEMD)
crypto visual hash --target sponge             # SHA-3 / SHAKE sponge construction
crypto visual hash --target blake3             # BLAKE3 binary-tree mode
crypto visual hash --target sha256-trace       # SHA-256 working-variable trace
crypto visual hash --target avalanche          # MD5 vs SHA-256 single-bit avalanche

# Elliptic curves
crypto visual ecc --target curve               # scatter plot of points on toy F_p curve
crypto visual ecc --target ops                 # point-add and point-double geometric diagrams
crypto visual ecc --target ladder              # double-and-add scalar-mul trace
crypto visual ecc --target ecdsa               # ECDSA sign/verify dataflow

# Everything at once
crypto visual all
```

### Visual output samples

The `aes-visual-demo` subcommand renders ASCII diagrams.  For example,
`--demo truncated-diff` shows the classical wide-trail AES diffusion:

```
        input      R1 SR      R2 SR      R3 SR      R4 SR
  row 0    ▓ · · ·    ▓ · · ·    ▓ · · ·    ▓ ▓ ▓ ▓    ▓ ▓ ▓ ▓
  row 1    · · · ·    · · · ·    · · · ▓    ▓ ▓ ▓ ▓    · · · ·
  row 2    · · · ·    · · · ·    · · ▓ ·    ▓ ▓ ▓ ▓    · · · ·
  row 3    · · · ·    · · · ·    · ▓ · ·    ▓ ▓ ▓ ▓    · · · ·

**Active S-boxes per round**

round  0 │█                              1
round  1 │███████                        4
round  2 │██████████████████████████████ 16
round  3 │███████                        4

**Total active S-boxes** N = 25, Pr ≤ 2⁻¹⁵⁰
```

The same machinery powers every other attack's test output — run any
`#[ignore]`-flagged `demo_*` test with `--nocapture` to see it.

### Elliptic-curve cryptography

| Module                | Algorithm / Curve                                           | Standard               |
|-----------------------|-------------------------------------------------------------|------------------------|
| `ecc::field` / `point`| Generic short-Weierstrass arithmetic                        | SEC 1                  |
| `ecc::curve`          | secp256k1, P-256                                            | SEC 2, FIPS 186-4      |
| `ecc::curve_zoo`      | Brainpool {160,192,224,256,320,384,512} r1, FRP256v1, NIST P-{192,224,384,521}, SECG k-curves, GOST tc26 paramsets | RFC 5639 / FIPS 186-4 / SEC 2 / GOST |
| `ecc::p256_field`     | Optimised P-256 field arithmetic (Solinas reduction)        | FIPS 186-4             |
| `ecc::secp256k1_*`    | Optimised secp256k1 field + Jacobian coordinates            | SEC 2                  |
| `ecc::ct`             | Constant-time scalar multiplication primitives              | —                      |
| `ecc::curve25519`     | Curve25519 prime field                                      | RFC 7748               |
| `ecc::x25519`         | X25519 ECDH                                                 | RFC 7748               |
| `ecc::ed25519`        | Ed25519 EdDSA signatures                                    | RFC 8032               |
| `ecc::ecdsa`          | ECDSA with RFC 6979 deterministic nonces                    | FIPS 186-4 + RFC 6979  |
| `ecc::barrett_ecdsa`  | ECDSA with Barrett-reduced fields                           | —                      |
| `ecc::ecdh`           | ECDH key exchange                                           | RFC 8422               |
| `ecc::schnorr`        | Schnorr signatures (BIP-340 flavour)                        | BIP-340                |
| `ecc::sm2`            | SM2 signature + encryption + ZA                             | GB/T 32918             |
| `ecc::gost_3410_2012` | GOST R 34.10-2012 signature (256/512)                       | GOST R 34.10-2012      |
| `ecc::ec_kcdsa`       | Korean EC-KCDSA (Q = x⁻¹·G convention)                       | TTAS.KO-12.0011        |
| `binary_ecc`          | Binary-field ECC over F_{2^m}                                | NIST B-curves          |

### Pairings

| Module        | Algorithm                                | Notes                       |
|---------------|------------------------------------------|-----------------------------|
| `bls12_381`   | BLS12-381 pairing-friendly curve         | Optimal Ate pairing, Fq/Fq2/Fq6/Fq12 towers |

### Asymmetric / public-key

| Module                      | Algorithm                          | Standard         |
|-----------------------------|------------------------------------|------------------|
| `asymmetric::rsa`           | RSA keygen (Miller-Rabin), enc, sig, CRT-CT decrypt | PKCS#1 / RFC 8017 |
| `asymmetric::paillier`      | Paillier additively-homomorphic encryption | Paillier 1999 |
| `asymmetric::elgamal`       | ElGamal encryption (multiplicative group) | ElGamal 1985 |
| `asymmetric::ec_elgamal`    | EC-ElGamal (point-on-curve homomorphic) | — |

### Key derivation / MAC

| Module          | Algorithm                       | Standard   |
|-----------------|---------------------------------|------------|
| `kdf::hkdf`     | HMAC-SHA256, HKDF extract+expand| RFC 5869   |
| `kdf::pbkdf2`   | PBKDF2-HMAC-SHA256              | RFC 8018   |
| `kdf::hmac_drbg`| HMAC-DRBG deterministic RNG     | NIST SP 800-90A |

### Post-quantum (educational)

| Module                       | Algorithm                                | Notes                              |
|------------------------------|------------------------------------------|------------------------------------|
| `pqc::kyber`                 | Kyber / ML-KEM (simplified)              | Not FIPS-203 compatible            |
| `pqc::frodo`                 | FrodoKEM (LWE-based)                     | Conservative simplified            |
| `pqc::ntru`                  | NTRU lattice cryptosystem                | NTRU encrypt 1996                  |
| `pqc::ntru_prime`            | NTRU Prime (Streamlined / NTRU LPRime)   | NTRU Prime spec                    |
| `pqc::mceliece`              | Educational binary-Goppa McEliece (toy params) | Original 1978 |
| `pqc::classic_mceliece`      | Classic McEliece (NIST round-3) skeleton | NIST round-3 spec                  |
| `pqc::hqc`                   | HQC (code-based KEM)                     | NIST round-3 / round-4             |
| `pqc::bike`                  | BIKE (QC-MDPC code-based)                | NIST round-3                       |
| `pqc::csidh`                 | CSIDH (isogeny-based) — group action     | Castryck–Lange–Martindale–Panny–Renes 2018 |
| `pqc::x_wing`                | X-Wing hybrid (X25519 + ML-KEM)          | draft-connolly-cfrg-xwing-kem      |

### Zero-knowledge / commitments

| Module                | Algorithm                                       |
|-----------------------|-------------------------------------------------|
| `zk::schnorr_zkp`     | Schnorr identification proof                    |
| `zk::pedersen`        | Pedersen commitments                            |
| `zk::chaum_pedersen`  | Chaum-Pedersen equality-of-discrete-logs proof  |
| `zk::merkle`          | Merkle tree commitments                         |
| `zk::polynomial`      | Polynomial arithmetic over a field              |
| `zk::kzg`             | KZG polynomial commitments                      |
| `zk::bulletproofs`    | Bulletproofs range proofs                       |
| `zk::groth16`         | Groth16 zk-SNARK skeleton                       |
| `zk::plonk`           | PLONK zk-SNARK skeleton                         |
| `zk::stark` / `stark_v2` | STARK proof system skeleton                   |

---

## Cryptanalysis suite

This is the part of the repo that's grown the most.  Every implemented algorithm
above has at least one attack module that *targets it*, plus a generic framework
that makes adding new attacks cheap.

### Core analytical primitives

| Module                                   | What it gives you |
|------------------------------------------|-------------------|
| `cryptanalysis::sbox`                    | DDT, LAT, BCT (Cid et al. 2018), DLCT (Bar-On et al. 2019), truncated DDT, differential / linear uniformity, nonlinearity, algebraic degree, boomerang uniformity |
| `cryptanalysis::boolean`                 | Walsh-Hadamard, ANF, algebraic degree |
| `cryptanalysis::avalanche`               | Full avalanche matrix, SAC, BIC |
| `cryptanalysis::statistical`             | Chi-squared, monobit, runs test |
| `cryptanalysis::lattice`                 | LLL, BKZ reduction |

### Symmetric-cipher attacks

| Module                                      | Attack                              |
|---------------------------------------------|-------------------------------------|
| `cryptanalysis::boomerang`                  | **Generic boomerang distinguisher + rectangle + sandwich + branch-and-bound trail search** (Wagner FSE 1999 → Cid et al. EUROCRYPT 2018) |
| `cryptanalysis::aes::reduced`               | Reduced-round AES-128 with configurable `Nr` |
| `cryptanalysis::aes::square`                | Square / integral attack (Daemen-Knudsen-Rijmen 1997) |
| `cryptanalysis::aes::differential`          | Differential cryptanalysis on 2-round AES |
| `cryptanalysis::aes::linear`                | Matsui-style linear cryptanalysis |
| `cryptanalysis::aes::impossible`            | Biham-Biryukov-Shamir impossible differential |
| `cryptanalysis::aes::mixture`               | Grassi-Rechberger-Rønjom mixture quadruples |
| `cryptanalysis::aes::boomerang`             | AES-specific boomerang + BCT |
| `cryptanalysis::aes::yoyo`                  | Rønjom-Bardeh-Helleseth yoyo distinguisher |
| `cryptanalysis::aes::mitm`                  | Demirci-Selçuk MITM on 2-round AES |
| `cryptanalysis::aes::small_scale`           | Cid-Murphy-Robshaw SR(n,2,2,4) variant |
| `cryptanalysis::aes::algebraic`             | Multivariate quadratic equations |
| `cryptanalysis::aes::milp`                  | MILP encoding of differential trail search |
| `cryptanalysis::aes::quantum_grover`        | Grover-style quantum key-recovery cost model |
| `cryptanalysis::aes::biclique`              | Biclique attack skeleton |
| `cryptanalysis::aes::related_key`           | **Related-key attack framework** — key-schedule difference propagation, related-key avalanche, related-key boomerang, Biryukov-Khovratovich local-collision demo |
| `cryptanalysis::aes::truncated_diff`        | **Truncated differential cryptanalysis** (Knudsen FSE 1994) — active-byte trail propagation, MixColumns branch-number filter, minimum-active-S-box bounds |
| `cryptanalysis::aes::higher_order`          | **Higher-order differential** (Lai 1994, Knudsen 1995) — d-th derivatives, 3-round Square integral distinguisher, algebraic-degree probe |
| `cryptanalysis::aes::dfa`                   | **DFA / Piret-Quisquater fault attack** (FDTC 2003) — single-byte fault before round 9 reduces last-round-key entropy by ~29 bits per fault |
| `cryptanalysis::aes::visualize`             | **ASCII visualization helpers** — state grids, active-byte patterns, trail diagrams, heat maps, bar charts, boomerang quartet diagrams, recovery progress bars |
| `cryptanalysis::md5_differential`           | MD5 differential collision (Wang et al. 2004) |
| `cryptanalysis::sha1_differential`          | SHA-1 differential analysis (Wang-Yin-Yu 2005) |

### Discrete-log / ECDLP attacks

| Module                                       | Attack                                                  |
|----------------------------------------------|---------------------------------------------------------|
| `cryptanalysis::pollard_rho`                 | Pollard ρ for DLP / ECDLP, multi-shard, distinguished-points |
| `cryptanalysis::preprocessing_rho`           | Bernstein-Lange precomputation rho                       |
| `cryptanalysis::ml_rho_walks`                | Pollard ρ walks under learned partition functions       |
| `cryptanalysis::aut_folded_rho`              | Automorphism-folded rho (CM curves)                      |
| `cryptanalysis::ec_index_calculus`           | Semaev S₃ index calculus on prime-field curves          |
| `cryptanalysis::ec_index_calculus_j0`        | ζ-orbit-reduced IC on j=0 curves + Eisenstein-smooth FB |
| `cryptanalysis::j0_twists`                   | 6-twist enumeration on j=0 curves + smoothness flagging |
| `cryptanalysis::canonical_lift`              | Smart attack on anomalous curves (canonical lifting)    |
| `cryptanalysis::cm_canonical_lift`           | CM-curve canonical lift + p-adic logarithm              |
| `cryptanalysis::nonanom_formal_log`          | Non-anomalous formal-group logarithm                    |
| `cryptanalysis::mazur_tate_sigma`            | Mazur-Tate σ function for anomalous curves              |
| `cryptanalysis::coleman_integration`         | Coleman integration on hyperelliptic curves             |
| `cryptanalysis::ai_schoof`                   | AI-assisted Schoof's algorithm                          |
| `cryptanalysis::ecm`                         | Lenstra elliptic-curve method (ECM) factoring           |
| `cryptanalysis::modular_polynomial`          | Modular polynomials Φ_ℓ(X, Y)                            |
| `cryptanalysis::hilbert_class_poly`          | Hilbert class polynomial computation                    |
| `cryptanalysis::fght_snfs`                   | FGHT toy SNFS pipeline for trapdoored F_p* DLP          |
| `cryptanalysis::shor`                        | Shor's quantum factoring + order-finding                |
| `cryptanalysis::quantum_estimator`           | Resource estimator for quantum attacks                  |

### Signature / nonce attacks

| Module                                      | Attack                                              |
|---------------------------------------------|-----------------------------------------------------|
| `cryptanalysis::hnp_ecdsa`                  | Hidden Number Problem on ECDSA (LLL / BKZ)          |
| `cryptanalysis::multi_key_hnp`              | Multi-key HNP master-key recovery                   |
| `cryptanalysis::bleichenbacher`             | Bleichenbacher direct bias attack (Mironov / Tibouchi) |
| `cryptanalysis::ecdsa_audit`                | ECDSA transcript audit + quick-bias score           |
| `cryptanalysis::signature_corpus`           | Corpus-level reuse / weak-RNG / bias detection      |
| `cryptanalysis::orbit_homology`             | Orbit-homology bias detector                        |
| `cryptanalysis::cga_hnc`                    | Continued-fraction / HNC variant                    |
| `cryptanalysis::b_seed_profile`             | Behaviour-seed profile attack                       |

### Structural / curve-specific

| Module                                | Attack                                                  |
|---------------------------------------|---------------------------------------------------------|
| `cryptanalysis::p256_attacks`         | P-256 specific attack vectors                            |
| `cryptanalysis::p256_structural`      | P-256 structural / lattice attacks                       |
| `cryptanalysis::p256_isogeny_cover`   | Isogeny-cover speculation targeting P-256                |
| `cryptanalysis::solinas_correlations` | Solinas-prime reduction correlations                     |

### Binary-curve Weil descent (GHS / Hess pipeline)

A full Weil-descent attack pipeline on binary-field elliptic curves.  Lifts
ECDLP on `E/F_{2^n}` to HCDLP on a hyperelliptic Jacobian over `F_{2^{n/m}}`,
then runs index calculus on the higher-genus Jacobian.

| Module                                | Role                                                     |
|---------------------------------------|----------------------------------------------------------|
| `cryptanalysis::ghs_descent`          | Frey-Rück / GHS cover-curve construction                 |
| `cryptanalysis::ghs_full_attack`      | End-to-end orchestrator producing structured `AttackReport` |
| `cryptanalysis::ec_trapdoor`          | EC trapdoor / weak-curve detection                       |
| `cryptanalysis::binary_isogeny`       | Vélu-style isogeny computation in characteristic 2       |
| `binary_ecc::hyperelliptic`           | Genus-`g` hyperelliptic curve arithmetic over F_{2^m}    |
| `binary_ecc::poly_f2m`                | Polynomial arithmetic over F_{2^m}                        |
| `prime_hyperelliptic`                 | Hyperelliptic curves over prime fields (Coleman / odd-char descent) |

End-to-end m=1 ECDLP recovery runs on F_{2^6}; m=2 type-II symbolic descent
demonstrated on F_{2^6} with (n=3, ℓ=2).  See `examples/ghs_attack_demo.rs`
for three runnable scenarios.

### Auto-attack framework + CLI

| Module                                | Purpose                                                 |
|---------------------------------------|---------------------------------------------------------|
| `cryptanalysis::cipher_registry`      | Named-cipher catalog (`toyspn-{1..4}r`, `aes-{1..4}r`, `toyspn-present-2r`) with `BlockCipher` adapter |
| `cryptanalysis::auto_attack`          | Discovers and runs every applicable attack against a named cipher; emits Markdown report |
| `cryptanalysis::research_bench`       | Falsifiable-hypothesis bench with `linear_fit` / `log_log_fit`, polynomial and exponential scaling regimes |

---

## Cryptanalysis CLI

```bash
# Discover what's registered
crypto cryptanalysis list-ciphers

# Run every applicable attack against a named cipher
crypto cryptanalysis auto --cipher toyspn-2r
crypto cryptanalysis auto --cipher aes-2r --boomerang-pairs 4096

# Drive individual attacks
crypto cryptanalysis sbox      --cipher aes-2r
crypto cryptanalysis boomerang --cipher toyspn-3r --pairs 65536
crypto cryptanalysis boomerang --cipher toyspn-2r --alpha 0100 --delta 0100
crypto cryptanalysis rectangle --cipher toyspn-1r --pool 2048

# Run the full research bench
crypto cryptanalysis bench
```

Auto-attack discovery logic (`auto --cipher <name>`):

| Section                  | Condition                          |
|--------------------------|------------------------------------|
| S-box DDT/LAT/BCT/DLCT   | cipher exposes a single S-box      |
| Boomerang distinguisher  | always                             |
| Rectangle attack         | always                             |
| Differential trail search| cipher is 4×4-bit SPN (`entry.trail_searchable`); skipped for AES (8-bit S-box outside the current SPN model — see `aes/milp.rs` for the MILP route) |

Sample output (`auto --cipher toyspn-2r`):

```
# Auto-cryptanalysis: `toyspn-2r`
**Block size**: 2 bytes (16 bits)
**Rounds**: 2
**Canonical (α, δ)**: 0x0100 / 0x0100

## S-box analysis
| differential uniformity | 4    |
| max DP                  | 0.25 |
| boomerang uniformity    | 16   |

## Boomerang distinguisher
- right quartets: 4099 / 16384
- empirical (pq)²: 2.50e-1
- random baseline: 1.53e-5
- distinguishes from random (10×): ✓

## Rectangle attack
- right rectangles: 2 (pool 1024)

## Differential trail search
| rank | input α | output ω | probability | log₂ p |
| 1    | 0x0001  | 0xccc0   | 7.8e-3      | -7.00  |

**Total**: 4 sections in 30 ms.
```

---

## Research bench

A continuously-updated empirical study of cryptanalytic-attack complexity.  Each
attack is registered as a falsifiable [`Hypothesis`] with a theoretical complexity
claim; the bench runs it at multiple problem sizes, fits the measured
exponent by log-log regression, and prints a verdict.

Two scaling regimes are supported:

- **Polynomial** (`time ∝ n^α`): default; ECDLP attacks, IC variants, twist enum.
- **Exponential** (`time ∝ 2^(c·rounds)`): boomerang attack on r-round ToySpn,
  fits bits-per-round directly via `fit_boomerang_decay`.

Run:

```bash
cargo test --lib --release bench_demo -- --ignored --nocapture
# or
cargo run --release -- cryptanalysis bench
```

Live measurements: see [`docs/RESEARCH_BENCH_LOG.md`](./docs/RESEARCH_BENCH_LOG.md).

Headline empirical findings (most recent):

| Hypothesis                          | Theoretical  | Measured | R²    | Verdict        |
|-------------------------------------|-------------:|---------:|------:|----------------|
| Pollard ρ ECDLP (7–20 bits)         | α = 0.500    | 0.807    | 0.996 | ≈ approx       |
| Index calculus (generic)            | α = 1.500    | 0.691    | 0.970 | ✗ toy regime   |
| IC on j=0 (ζ-orbit reduced)         | α = 1.500    | 0.752    | 0.895 | ✗ toy regime   |
| Eisenstein-smooth IC (j=0)          | α = 1.500    | 0.825    | 0.992 | ✗ falsified speculation |
| j=0 twist enumeration               | α = 1.000    | 1.035    | 0.999 | ✓ matches      |
| Boomerang decay (ToySpn r=1..4)     | ~8 bits/round| 5.05     | 0.990 | empirical signature |

See also [`docs/ECDLP_ATTACK_MATRIX.md`](./docs/ECDLP_ATTACK_MATRIX.md) for the
attack/curve-family applicability matrix.

---

## CLI demo for non-cryptanalysis features

```bash
# Hash
crypto hash sha256 "hello world"
crypto hash sha3   "hello world"
crypto hash blake3 "hello world"

# ECC / ECDSA / ECDH
crypto ecc keygen
crypto ecdsa sign   "my message"
crypto ecdsa verify "my message" <pubkey-hex> <r-hex> <s-hex>
crypto ecdh

# RSA (slow — 1024-bit keygen)
crypto rsa keygen
crypto rsa demo

# Symmetric AEAD
crypto aes
crypto chacha

# Key derivation
crypto hkdf

# Post-quantum
crypto pqc

# Run all of the above demos at once
crypto demo
```

---

## Mathematics background

### Elliptic curves

An elliptic curve over a prime field **F_p** is the set of points satisfying:

```
y² ≡ x³ + ax + b  (mod p)
```

together with a special "point at infinity" **O** that acts as the additive identity.

**Point addition** (P ≠ Q):  λ = (y₂ − y₁) / (x₂ − x₁) mod p ; x₃ = λ² − x₁ − x₂ ; y₃ = λ(x₁ − x₃) − y₁.

**Point doubling** (P = Q):  λ = (3x₁² + a) / (2y₁) mod p.

**Scalar multiplication** k·P uses double-and-add in `O(log k)` group operations.

#### Curve families we support

- **Bitcoin/Ethereum**: secp256k1 (a=0, b=7, p ≈ 2²⁵⁶).
- **NIST / TLS**: P-192, P-224, P-256, P-384, P-521.
- **Brainpool** (RFC 5639): P{160,192,224,256,320,384,512}r1.
- **French ANSSI**: FRP256v1.
- **Chinese**: SM2 (`a = −3 mod p`).
- **Russian GOST**: tc26-256 and tc26-512 paramsets.
- **Korean**: EC-KCDSA convention `Q = x⁻¹·G`.
- **CFRG / RFC 7748**: Curve25519 (Montgomery) + Ed25519 (Edwards).

### AES

AES operates on a 4×4 byte state through `Nr` rounds (10/14 for AES-128/256).
Each round applies:

1. **SubBytes** — non-linear S-box (inverse in GF(2⁸) + affine map).
2. **ShiftRows** — cyclic rotation of rows 0..3 by 0..3 positions.
3. **MixColumns** — matrix mul over GF(2⁸) (omitted in the final round).
4. **AddRoundKey** — XOR with the round key.

**GCM** adds GHASH (polynomial mul over GF(2¹²⁸)) for authentication.

The cryptanalysis side studies all four operations: the S-box DDT/LAT/BCT
(`cryptanalysis::sbox`), the MixColumns branch number = 5 (`aes::differential`),
the Square/integral structure (`aes::square`), etc.

### Boomerang attack

Split `E = E₁ ∘ E₀`.  Find differential α →_{E₀} β with prob p, γ →_{E₁} δ with
prob q.  For random P: encrypt (P, P⊕α) to (C, C⊕β?); XOR δ onto both; decrypt.
Claim: P' ⊕ P'' = α with probability ≈ (pq)² (Wagner FSE 1999).

The **Boomerang Connectivity Table** (Cid–Huang–Peyrin–Sasaki–Song EUROCRYPT 2018)
corrects this naive estimate — at the boundary the two trails interact through
the S-box and the BCT entry replaces `(pq)²` with `(p²·q²·BCT[β][γ]/2^n)`.  We
implement BCT (per-S-box) and the **sandwich attack** (Dunkelman-Keller-Shamir
CRYPTO 2010) that uses BCT directly.

### Pollard ρ for DLP

Random walk on the cyclic group; collision via Floyd's tortoise-and-hare.
Expected `√(πn/2)` group operations to recover the discrete log.

### Index calculus (ECDLP, prime fields)

Semaev's summation polynomials `S_n(x_1, …, x_n)` vanish iff there exist
y-coords making (x_i, y_i) on the curve sum to zero.  We use `S_3` for
2-decompositions; cost `O(p · m)` per relation with factor-base size `m`,
optimum at `m = √p` giving `O(p^{3/2})` — asymptotically worse than ρ on prime
curves, but the bench studies the empirical exponent.

### Hidden Number Problem

ECDSA with biased nonces `k = k_low + 2^B · k_high` leaks `k_high` through
linear relations; LLL on the right lattice recovers the secret in polynomial
time given enough biased signatures.  See `hnp_ecdsa`, `multi_key_hnp`,
`bleichenbacher`.

---

## Project structure

```
src/
├── lib.rs / main.rs           — library root + CLI
├── ecc/                       — ECC core + curve zoo + ECDSA/ECDH/Schnorr/SM2/GOST/EC-KCDSA
├── binary_ecc/                — Binary-field ECC + polynomial-over-F_{2^m} + hyperelliptic Jacobians
├── prime_hyperelliptic/       — Odd-char hyperelliptic curves (for Coleman / descent)
├── bls12_381/                 — Pairing-friendly BLS12-381 curve + pairing
├── symmetric/                 — AES, ChaCha, Serpent, Threefish, SM4, Magma, Kuznyechik
├── hash/                      — SHA-2/3, BLAKE3, RIPEMD-160, SipHash, SM3, Streebog
├── asymmetric/                — RSA, Paillier, ElGamal, EC-ElGamal
├── pqc/                       — Kyber, Frodo, NTRU/NTRU-Prime, McEliece, HQC, BIKE, CSIDH, X-Wing
├── kdf/                       — HKDF, PBKDF2, HMAC-DRBG
├── zk/                        — Schnorr, Pedersen, Chaum-Pedersen, Merkle, KZG, Bulletproofs, Groth16, PLONK, STARK
├── cryptanalysis/             — see Cryptanalysis suite section above
│   ├── aes/                   — Reduced-round + Square + linear + diff + impossible + mixture + boomerang + yoyo + MITM + small-scale + algebraic + MILP + Grover + biclique
│   ├── boomerang.rs           — Generic distinguisher + rectangle + sandwich + trail search
│   ├── ghs_descent.rs         — Binary-curve Weil descent (cover construction)
│   ├── ghs_full_attack.rs     — GHS attack orchestrator
│   ├── binary_isogeny.rs      — Char-2 Vélu isogenies
│   ├── ec_trapdoor.rs         — Weak-curve / trapdoor detection
│   ├── cipher_registry.rs     — Named-cipher catalog
│   ├── auto_attack.rs         — Auto-discovery + dispatch
│   ├── research_bench.rs      — Falsifiable-hypothesis bench
│   └── …45+ other attack modules
├── examples/
│   └── ghs_attack_demo.rs     — Three runnable GHS scenarios
├── ct_bignum.rs               — Constant-time big-integer experiments
├── ecc_safety.rs              — ECC parameter-safety auditor
└── utils/                     — Modular arithmetic, encoding, randomness
```

---

## Documentation

- [`SECURITY.md`](./SECURITY.md) — structural limitations + recommended alternatives.
- [`RESEARCH.md`](./RESEARCH.md) — research notes.
- [`RESEARCH_P256.md`](./RESEARCH_P256.md) — P-256 specific structural studies.
- [`DEFERRED.md`](./DEFERRED.md) — known gaps + deferred work.
- [`docs/ECDLP_ATTACK_MATRIX.md`](./docs/ECDLP_ATTACK_MATRIX.md) — ECDLP attack taxonomy.
- [`docs/RESEARCH_BENCH_LOG.md`](./docs/RESEARCH_BENCH_LOG.md) — live empirical bench measurements.

---

## References

### Standards

- FIPS 197 — AES
- FIPS 180-4 — SHA-2
- FIPS 202 — SHA-3 / Keccak
- FIPS 203 — ML-KEM (Kyber)
- FIPS 186-4 — ECDSA + NIST curves
- RFC 8439 — ChaCha20-Poly1305
- RFC 7748 — X25519 / X448
- RFC 8032 — EdDSA / Ed25519
- RFC 5869 — HKDF
- RFC 8018 — PBKDF2
- RFC 5639 — Brainpool curves
- RFC 6979 — Deterministic ECDSA nonces
- SEC 1 / SEC 2 — Elliptic curve standards (Certicom)
- GB/T 32905 / 32907 / 32918 — SM3 / SM4 / SM2
- GOST R 34.10-2012, 34.11-2012, 34.12-2015 — Russian GOST suite

### Cryptanalysis literature

- **D. Wagner**, *The Boomerang Attack*, FSE 1999.
- **E. Biham, O. Dunkelman, N. Keller**, *The Rectangle Attack*, EUROCRYPT 2001.
- **O. Dunkelman, N. Keller, A. Shamir**, *A Practical-Time Related-Key Attack
  on KASUMI*, CRYPTO 2010.
- **C. Cid, T. Huang, T. Peyrin, Y. Sasaki, L. Song**, *Boomerang Connectivity
  Table*, EUROCRYPT 2018.
- **A. Bar-On, O. Dunkelman, N. Keller, A. Weizman**, *DLCT: A New Tool for
  Differential-Linear Cryptanalysis*, EUROCRYPT 2019.
- **J. Daemen, L. Knudsen, V. Rijmen**, *The Block Cipher Square*, FSE 1997.
- **E. Biham, A. Shamir**, *Differential Cryptanalysis of DES-like
  Cryptosystems*, CRYPTO 1990.
- **M. Matsui**, *Linear Cryptanalysis Method for DES Cipher*, EUROCRYPT 1993.
- **I. Semaev**, *Summation polynomials and the discrete logarithm problem
  on elliptic curves*, ePrint 2004/031.
- **C. Diem**, *On the discrete logarithm problem in elliptic curves*,
  Compositio Math. 147 (2011).
- **N. Smart**, *The discrete logarithm problem on elliptic curves of trace
  one*, J. Cryptology 12 (1999).
- **P. Q. Nguyen, I. Shparlinski**, *The insecurity of the digital signature
  algorithm with partially known nonces*, J. Cryptology 15 (2002).
- **D. Bleichenbacher**, *On the generation of one-time keys in DL signature
  schemes*, IEEE P1363 submission (2000).
- **P. W. Shor**, *Polynomial-time algorithms for prime factorization and
  discrete logarithms on a quantum computer*, SIAM 1997.
- **A. K. Lenstra, H. W. Lenstra Jr., L. Lovász**, *Factoring polynomials
  with rational coefficients*, Math. Ann. 261 (1982).
- **Fried, Gaudry, Heninger, Thomé**, *A kilobit hidden SNFS discrete
  logarithm computation*, EUROCRYPT 2017.
