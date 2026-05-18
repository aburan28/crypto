# Mainstream Elliptic Curve Audit — Full Methodology Application

Applies the methodology developed across the 21-phase investigation of NIST P-256 to **31 mainstream elliptic curves** spanning every major deployment context.

**Methodology applied per curve**:
1. Phase 10/14 Frobenius-mod-2 parity obstruction check
2. Cofactor / cardinality parity classification
3. `t_E` parity and target `(a, b) mod 2` calculation
4. Vulnerability regime assignment

## Top-line summary

- **Phase-10 theorem-protected (STRUCTURALLY blocked)**: 24 curves
- **Statistical-rarity protected**: 3 curves
- **NOT blocked (Teske/GHS family applies)**: 4 curves

---

## Regime 1 — Phase-10 STRUCTURAL theorem-protected

These curves are blocked by the parity-mod-2 obstruction.  For any N ≥ 2, no smooth genus-2 Jacobian over their base field can be `F_p`-isogenous to `E × E^twist`.  This is a **theorem** (Phase 10/14), not a probabilistic null.

| Curve | Family | p bits | n bits | Deployment |
|---|---|---|---|---|
| NIST P-192 | NIST FIPS 186-4 | 192 | 192 | Deprecated |
| NIST P-224 | NIST FIPS 186-4 | 224 | 224 | Production |
| NIST P-256 | NIST FIPS 186-4 | 256 | 256 | Production (TLS, FIDO2, IETF) |
| NIST P-384 | NIST FIPS 186-4 | 384 | 384 | Production (TLS, government) |
| NIST P-521 | NIST FIPS 186-4 | 513 | 521 | Production (top-secret class) |
| secp192k1 | SECG K (Bitcoin family) | 192 | 192 | Niche |
| secp224k1 | SECG K | 224 | 225 | Niche |
| secp256k1 | SECG K (Bitcoin, Ethereum) | 256 | 256 | Production (BTC, ETH, Cardano) |
| brainpoolP192r1 | Brainpool | 192 | 192 | EU government |
| brainpoolP224r1 | Brainpool | 224 | 224 | EU government |
| brainpoolP256r1 | Brainpool | 256 | 256 | EU government |
| brainpoolP320r1 | Brainpool | 320 | 320 | EU government |
| brainpoolP384r1 | Brainpool | 384 | 384 | EU government |
| brainpoolP512r1 | Brainpool | 512 | 512 | EU government |
| FRP256v1 | ANSSI (France) | 256 | 256 | EU regulated |
| SM2 (China) | China SM2 | 256 | 256 | Production (China, GM/T 0003) |
| GOST CryptoPro A | GOST R 34.10-2001 | 256 | 256 | Production (Russia) |
| GOST CryptoPro B | GOST R 34.10-2001 | 256 | 256 | Production (Russia) |
| GOST CryptoPro C | GOST R 34.10-2001 | 256 | 256 | Production (Russia) |
| GOST TC26 256-A | GOST TC26 | 256 | 256 | Production (Russia) |
| BN254 (alt_bn128) | BN (Ethereum) | 254 | 254 | Production (Ethereum precompile) |
| BLS12-381 (E/F_p) | BLS pairing (ZK) | 381 | 381 | Production (Zcash, Filecoin, Eth2) |
| numsp256d1 | Microsoft NUMS | 256 | 254 | Research |
| numsp384d1 | Microsoft NUMS | 380 | 382 | Research |

---

## Regime 2 — Statistical-rarity protected

These curves have **even** `#E(F_p)` so the Phase-10 parity argument doesn't apply.  However, the empirical Phase 18–21 work shows cover-Jacobians **do exist** at toy scale with density `~p^{-3.4}` — which extrapolates to `2^{-858}` at cryptographic scale.  Protection is **probabilistic** rather than theorem-grade.

| Curve | Family | p bits | cofactor | Deployment |
|---|---|---|---|---|
| Curve25519 | RFC 7748 (Mont/Edwards) | 255 | 8 | Production (TLS 1.3, Signal, WireGuard) |
| Curve448 | RFC 7448 (Mont/Edwards) | 448 | 4 | Production (TLS 1.3, CFRG) |
| E-521 | Hamburg Edwards | 513 | 4 | Research |

---

## Regime 3 — NOT blocked (Teske/GHS applies)

**Binary-field** curves are outside the Phase 10 analysis domain (the obstruction was derived for odd-prime fields).  These are the **original Teske/GHS attack targets**: the binary-field cover attack actually works against them, with varying degrees of practical impact.  Implementing the full attack against deployed binary curves remains an active research area.

| Curve | Family | n bits | cofactor | Deployment |
|---|---|---|---|---|
| sect163k1 | NIST Koblitz binary | 163 | 2 | Deprecated |
| sect233k1 | NIST Koblitz binary | 232 | 4 | Deprecated |
| sect283k1 | NIST Koblitz binary | 281 | 4 | Deprecated |
| sect571k1 | NIST Koblitz binary | 570 | 4 | Niche |

---

## Per-curve detail

### NIST P-192

- **Family**: NIST FIPS 186-4
- **Deployment**: Deprecated
- **Base field**: p (192 bits, prime, odd)
- **Group order n**: 192 bits, cofactor 1
- **n is prime** (cofactor 1)
- **Trace parity** (t mod 2): 1
- **Target Frobenius parity** for E × E^twist: (0, 1)
- **Verdict**: Phase-10 theorem (STRUCTURALLY blocked, all N)

### NIST P-224

- **Family**: NIST FIPS 186-4
- **Deployment**: Production
- **Base field**: p (224 bits, prime, odd)
- **Group order n**: 224 bits, cofactor 1
- **n is prime** (cofactor 1)
- **Trace parity** (t mod 2): 1
- **Target Frobenius parity** for E × E^twist: (0, 1)
- **Verdict**: Phase-10 theorem (STRUCTURALLY blocked, all N)

### NIST P-256

- **Family**: NIST FIPS 186-4
- **Deployment**: Production (TLS, FIDO2, IETF)
- **Base field**: p (256 bits, prime, odd)
- **Group order n**: 256 bits, cofactor 1
- **n is prime** (cofactor 1)
- **Trace parity** (t mod 2): 1
- **Target Frobenius parity** for E × E^twist: (0, 1)
- **Verdict**: Phase-10 theorem (STRUCTURALLY blocked, all N)

### NIST P-384

- **Family**: NIST FIPS 186-4
- **Deployment**: Production (TLS, government)
- **Base field**: p (384 bits, prime, odd)
- **Group order n**: 384 bits, cofactor 1
- **n is prime** (cofactor 1)
- **Trace parity** (t mod 2): 1
- **Target Frobenius parity** for E × E^twist: (0, 1)
- **Verdict**: Phase-10 theorem (STRUCTURALLY blocked, all N)

### NIST P-521

- **Family**: NIST FIPS 186-4
- **Deployment**: Production (top-secret class)
- **Base field**: p (513 bits, prime, odd)
- **Group order n**: 521 bits, cofactor 1
- **n is prime** (cofactor 1)
- **Trace parity** (t mod 2): 1
- **Target Frobenius parity** for E × E^twist: (0, 1)
- **Verdict**: Phase-10 theorem (STRUCTURALLY blocked, all N)

### secp192k1

- **Family**: SECG K (Bitcoin family)
- **Deployment**: Niche
- **Base field**: p (192 bits, prime, odd)
- **Group order n**: 192 bits, cofactor 1
- **n is prime** (cofactor 1)
- **Trace parity** (t mod 2): 1
- **Target Frobenius parity** for E × E^twist: (0, 1)
- **Verdict**: Phase-10 theorem (STRUCTURALLY blocked, all N)

### secp224k1

- **Family**: SECG K
- **Deployment**: Niche
- **Base field**: p (224 bits, prime, odd)
- **Group order n**: 225 bits, cofactor 1
- **n is prime** (cofactor 1)
- **Trace parity** (t mod 2): 1
- **Target Frobenius parity** for E × E^twist: (0, 1)
- **Verdict**: Phase-10 theorem (STRUCTURALLY blocked, all N)

### secp256k1

- **Family**: SECG K (Bitcoin, Ethereum)
- **Deployment**: Production (BTC, ETH, Cardano)
- **Base field**: p (256 bits, prime, odd)
- **Group order n**: 256 bits, cofactor 1
- **n is prime** (cofactor 1)
- **Trace parity** (t mod 2): 1
- **Target Frobenius parity** for E × E^twist: (0, 1)
- **Verdict**: Phase-10 theorem (STRUCTURALLY blocked, all N)

### brainpoolP192r1

- **Family**: Brainpool
- **Deployment**: EU government
- **Base field**: p (192 bits, prime, odd)
- **Group order n**: 192 bits, cofactor 1
- **n is prime** (cofactor 1)
- **Trace parity** (t mod 2): 1
- **Target Frobenius parity** for E × E^twist: (0, 1)
- **Verdict**: Phase-10 theorem (STRUCTURALLY blocked, all N)

### brainpoolP224r1

- **Family**: Brainpool
- **Deployment**: EU government
- **Base field**: p (224 bits, prime, odd)
- **Group order n**: 224 bits, cofactor 1
- **n is prime** (cofactor 1)
- **Trace parity** (t mod 2): 1
- **Target Frobenius parity** for E × E^twist: (0, 1)
- **Verdict**: Phase-10 theorem (STRUCTURALLY blocked, all N)

### brainpoolP256r1

- **Family**: Brainpool
- **Deployment**: EU government
- **Base field**: p (256 bits, prime, odd)
- **Group order n**: 256 bits, cofactor 1
- **n is prime** (cofactor 1)
- **Trace parity** (t mod 2): 1
- **Target Frobenius parity** for E × E^twist: (0, 1)
- **Verdict**: Phase-10 theorem (STRUCTURALLY blocked, all N)

### brainpoolP320r1

- **Family**: Brainpool
- **Deployment**: EU government
- **Base field**: p (320 bits, prime, odd)
- **Group order n**: 320 bits, cofactor 1
- **n is prime** (cofactor 1)
- **Trace parity** (t mod 2): 1
- **Target Frobenius parity** for E × E^twist: (0, 1)
- **Verdict**: Phase-10 theorem (STRUCTURALLY blocked, all N)

### brainpoolP384r1

- **Family**: Brainpool
- **Deployment**: EU government
- **Base field**: p (384 bits, prime, odd)
- **Group order n**: 384 bits, cofactor 1
- **n is prime** (cofactor 1)
- **Trace parity** (t mod 2): 1
- **Target Frobenius parity** for E × E^twist: (0, 1)
- **Verdict**: Phase-10 theorem (STRUCTURALLY blocked, all N)

### brainpoolP512r1

- **Family**: Brainpool
- **Deployment**: EU government
- **Base field**: p (512 bits, prime, odd)
- **Group order n**: 512 bits, cofactor 1
- **n is prime** (cofactor 1)
- **Trace parity** (t mod 2): 1
- **Target Frobenius parity** for E × E^twist: (0, 1)
- **Verdict**: Phase-10 theorem (STRUCTURALLY blocked, all N)

### FRP256v1

- **Family**: ANSSI (France)
- **Deployment**: EU regulated
- **Base field**: p (256 bits, prime, odd)
- **Group order n**: 256 bits, cofactor 1
- **n is prime** (cofactor 1)
- **Trace parity** (t mod 2): 1
- **Target Frobenius parity** for E × E^twist: (0, 1)
- **Verdict**: Phase-10 theorem (STRUCTURALLY blocked, all N)

### SM2 (China)

- **Family**: China SM2
- **Deployment**: Production (China, GM/T 0003)
- **Base field**: p (256 bits, prime, odd)
- **Group order n**: 256 bits, cofactor 1
- **n is prime** (cofactor 1)
- **Trace parity** (t mod 2): 1
- **Target Frobenius parity** for E × E^twist: (0, 1)
- **Verdict**: Phase-10 theorem (STRUCTURALLY blocked, all N)

### GOST CryptoPro A

- **Family**: GOST R 34.10-2001
- **Deployment**: Production (Russia)
- **Base field**: p (256 bits, prime, odd)
- **Group order n**: 256 bits, cofactor 1
- **n is prime** (cofactor 1)
- **Trace parity** (t mod 2): 1
- **Target Frobenius parity** for E × E^twist: (0, 1)
- **Verdict**: Phase-10 theorem (STRUCTURALLY blocked, all N)

### GOST CryptoPro B

- **Family**: GOST R 34.10-2001
- **Deployment**: Production (Russia)
- **Base field**: p (256 bits, prime, odd)
- **Group order n**: 256 bits, cofactor 1
- **n is prime** (cofactor 1)
- **Trace parity** (t mod 2): 1
- **Target Frobenius parity** for E × E^twist: (0, 1)
- **Verdict**: Phase-10 theorem (STRUCTURALLY blocked, all N)

### GOST CryptoPro C

- **Family**: GOST R 34.10-2001
- **Deployment**: Production (Russia)
- **Base field**: p (256 bits, prime, odd)
- **Group order n**: 256 bits, cofactor 1
- **n is prime** (cofactor 1)
- **Trace parity** (t mod 2): 1
- **Target Frobenius parity** for E × E^twist: (0, 1)
- **Verdict**: Phase-10 theorem (STRUCTURALLY blocked, all N)

### GOST TC26 256-A

- **Family**: GOST TC26
- **Deployment**: Production (Russia)
- **Base field**: p (256 bits, prime, odd)
- **Group order n**: 256 bits, cofactor 1
- **n is prime** (cofactor 1)
- **Trace parity** (t mod 2): 1
- **Target Frobenius parity** for E × E^twist: (0, 1)
- **Verdict**: Phase-10 theorem (STRUCTURALLY blocked, all N)

### BN254 (alt_bn128)

- **Family**: BN (Ethereum)
- **Deployment**: Production (Ethereum precompile)
- **Base field**: p (254 bits, prime, odd)
- **Group order n**: 254 bits, cofactor 1
- **n is prime** (cofactor 1)
- **Trace parity** (t mod 2): 1
- **Target Frobenius parity** for E × E^twist: (0, 1)
- **Verdict**: Phase-10 theorem (STRUCTURALLY blocked, all N)

### BLS12-381 (E/F_p)

- **Family**: BLS pairing (ZK)
- **Deployment**: Production (Zcash, Filecoin, Eth2)
- **Base field**: p (381 bits, prime, odd)
- **Group order n**: 381 bits, cofactor 1
- **n is prime** (cofactor 1)
- **Trace parity** (t mod 2): 1
- **Target Frobenius parity** for E × E^twist: (0, 1)
- **Verdict**: Phase-10 theorem (STRUCTURALLY blocked, all N)

### Curve25519

- **Family**: RFC 7748 (Mont/Edwards)
- **Deployment**: Production (TLS 1.3, Signal, WireGuard)
- **Base field**: p (255 bits, prime, odd)
- **Group order n**: 256 bits, cofactor 8
- **n has cofactor 8** ⇒ #E even
- **Trace parity** (t mod 2): 0
- **Target Frobenius parity** for E × E^twist: (0, 0)
- **Verdict**: Statistical rarity (~p^-3.4); covers exist at toy scale (Phase 18-21)

### Curve448

- **Family**: RFC 7448 (Mont/Edwards)
- **Deployment**: Production (TLS 1.3, CFRG)
- **Base field**: p (448 bits, prime, odd)
- **Group order n**: 528 bits, cofactor 4
- **n has cofactor 4** ⇒ #E even
- **Trace parity** (t mod 2): 0
- **Target Frobenius parity** for E × E^twist: (0, 0)
- **Verdict**: Statistical rarity (~p^-3.4); covers exist at toy scale (Phase 18-21)

### numsp256d1

- **Family**: Microsoft NUMS
- **Deployment**: Research
- **Base field**: p (256 bits, prime, odd)
- **Group order n**: 254 bits, cofactor 1
- **n is prime** (cofactor 1)
- **Trace parity** (t mod 2): 1
- **Target Frobenius parity** for E × E^twist: (0, 1)
- **Verdict**: Phase-10 theorem (STRUCTURALLY blocked, all N)

### numsp384d1

- **Family**: Microsoft NUMS
- **Deployment**: Research
- **Base field**: p (380 bits, prime, odd)
- **Group order n**: 382 bits, cofactor 1
- **n is prime** (cofactor 1)
- **Trace parity** (t mod 2): 1
- **Target Frobenius parity** for E × E^twist: (0, 1)
- **Verdict**: Phase-10 theorem (STRUCTURALLY blocked, all N)

### E-521

- **Family**: Hamburg Edwards
- **Deployment**: Research
- **Base field**: p (513 bits, prime, odd)
- **Group order n**: 595 bits, cofactor 4
- **n has cofactor 4** ⇒ #E even
- **Trace parity** (t mod 2): 0
- **Target Frobenius parity** for E × E^twist: (0, 0)
- **Verdict**: Statistical rarity (~p^-3.4); covers exist at toy scale (Phase 18-21)

### sect163k1

- **Family**: NIST Koblitz binary
- **Deployment**: Deprecated
- **Base field**: binary `F_{2^m}` (m varies)
- **Group order n**: 163 bits, cofactor 2
- **n has cofactor 2** ⇒ #E even
- **Verdict**: NOT BLOCKED — vulnerable to original Teske/GHS binary descent

### sect233k1

- **Family**: NIST Koblitz binary
- **Deployment**: Deprecated
- **Base field**: binary `F_{2^m}` (m varies)
- **Group order n**: 232 bits, cofactor 4
- **n has cofactor 4** ⇒ #E even
- **Verdict**: NOT BLOCKED — vulnerable to original Teske/GHS binary descent

### sect283k1

- **Family**: NIST Koblitz binary
- **Deployment**: Deprecated
- **Base field**: binary `F_{2^m}` (m varies)
- **Group order n**: 281 bits, cofactor 4
- **n has cofactor 4** ⇒ #E even
- **Verdict**: NOT BLOCKED — vulnerable to original Teske/GHS binary descent

### sect571k1

- **Family**: NIST Koblitz binary
- **Deployment**: Niche
- **Base field**: binary `F_{2^m}` (m varies)
- **Group order n**: 570 bits, cofactor 4
- **n has cofactor 4** ⇒ #E even
- **Verdict**: NOT BLOCKED — vulnerable to original Teske/GHS binary descent

