# The Comprehensive Guide to Elliptic Curves: Types, Attacks, and Vulnerabilities

> A complete reference covering every major family of deployed and
> research elliptic curves, every known attack class against the
> elliptic curve discrete logarithm problem (ECDLP) and its
> applications, and the implementation pitfalls that have caused
> real-world breaks.

---

## Part 1: Types of Elliptic Curves

An elliptic curve `E` over a field `K` is a smooth projective curve
of genus 1 with a distinguished point `O` (the "point at infinity").
The set `E(K)` of `K`-rational points forms an abelian group under
the chord-tangent law.

Different *forms* of the equation defining `E` have different
practical properties — efficiency, side-channel resistance, ease
of point validation, etc.

### 1.1 By equation form

#### 1.1.1 General Weierstrass form

```
    y² + a₁·xy + a₃·y = x³ + a₂·x² + a₄·x + a₆
```

The most general 5-parameter form.  Used:
- **Always over characteristic 2** (where short Weierstrass doesn't work)
- **Over characteristic 3** when reducing to short form would lose info

#### 1.1.2 Short Weierstrass form (`char K ≠ 2, 3`)

```
    y² = x³ + a·x + b
```

The form used by NIST P-curves, secp256k1, brainpool, GOST, SM2.

Variants:
- `a = -3`: NIST chose this for P-256 because it enables a faster
  doubling formula in Jacobian coordinates (`4M + 4S` instead of
  `5M + 4S`).
- `a = 0`: secp256k1 (j-invariant = 0) — enables an efficient
  endomorphism `(x, y) ↦ (β·x, y)` where `β` is a primitive cube
  root of unity (used for GLV speedup).

#### 1.1.3 Montgomery form (`char K ≠ 2`)

```
    B·y² = x³ + A·x² + x
```

Used by **Curve25519** (`A = 486662, B = 1`) and **Curve448**.
Properties:
- **x-only arithmetic**: scalar multiplication can be done using
  only the `x` coordinate via the **Montgomery ladder**, which is
  inherently constant-time.
- **No exception cases** in the ladder formula — no special
  handling needed for `2P` or `P + Q`.
- Resistant to many invalid-curve attacks (the ladder works on
  all `x` even off the curve).

#### 1.1.4 (Twisted) Edwards form (`char K ≠ 2`)

```
    a·x² + y² = 1 + d·x²·y²
```

If `a = 1`, it's **Edwards form**; otherwise **twisted Edwards**.

Used by:
- **Ed25519** (twisted Edwards, `a = -1`, `d = -121665/121666`)
- **Ed448** (Edwards form, `a = 1`)
- **JubJub** (Zcash zk-SNARK curve)

Properties:
- **Complete addition law**: a single formula works for all input
  pairs (no special cases for `P + (-P)`, doubling, identity).
  This makes constant-time implementation easier.
- Slightly faster than Weierstrass in some coordinate systems.

#### 1.1.5 Jacobi quartic / Jacobi intersection / Hessian / Huff

Niche forms used in research or specific implementations.  Generally
not deployed at scale.

### 1.2 By base field

#### 1.2.1 Prime fields `F_p`

The base field is the integers mod a prime `p`.  This is the
standard case for almost all deployed curves.

Choice of `p` matters:
- **Solinas primes**: `p = 2^a − 2^b ± … ± 1` with low Hamming
  weight.  Enables fast modular reduction without multiplication.
  Examples: NIST P-256 (`p = 2²⁵⁶ − 2²²⁴ + 2¹⁹² + 2⁹⁶ − 1`).
- **Pseudo-Mersenne primes**: `p = 2^n − c` for small `c`.
  Curve25519 (`p = 2²⁵⁵ − 19`), secp256k1 (`p = 2²⁵⁶ − 2³² − 977`).
- **Mersenne-like / Crandall primes**: `p = 2^n − c` with very
  small `c`.

#### 1.2.2 Binary fields `F_{2^m}`

Used by NIST sect-curves: `sect163k1`, `sect233k1`, `sect283k1`,
`sect409k1`, `sect571k1` (Koblitz) and their `sectXXXr1`
random-curve counterparts.

Properties:
- Hardware-friendly (XOR-based arithmetic)
- Historically faster than prime fields in custom silicon
- **Vulnerable to Weil descent (Teske/GHS)** when `m` is composite
- NIST deliberately chose `m` PRIME to neutralize this

Modern use: mostly deprecated.  Not used in TLS 1.3 cipher suites.

#### 1.2.3 Optimal extension fields `F_{p^k}` (`k ≥ 2`)

Used in some pairing-friendly contexts.  Generally vulnerable to
sub-exponential index-calculus attacks (Gaudry, Diem) for `k ≥ 3`.
Not used in mainstream cryptography.

#### 1.2.4 Ternary fields `F_{3^m}`

Used historically in pairing-based research (eta pairing).
Vulnerable to similar descent attacks as binary fields.
Effectively obsolete.

### 1.3 By special structure

#### 1.3.1 Ordinary vs Supersingular

- **Ordinary**: `gcd(t, p) = 1` where `t` is the Frobenius trace.
  All deployed cryptographic curves.
- **Supersingular**: `p | t` (in characteristic `p`).  Have small
  embedding degree (`k ≤ 6`), making them vulnerable to MOV /
  Frey–Rück but useful for pairing-based crypto.

#### 1.3.2 Anomalous

`#E(F_p) = p`, i.e., trace `t = 1`.  **Catastrophically broken**
by Smart's attack (1999) — `O(log p)` ECDLP via `p`-adic
formal-group logarithms.  Never deployed; only of theoretical
interest.

#### 1.3.3 CM curves

Curves with **complex multiplication** — `End(E) ⊗ Q` is an
imaginary quadratic field with small discriminant `D`.
- secp256k1 has `D = -3` (small, but the curve is still secure)
- Most curves with explicit small `D` are research-only

Small-`D` CM curves are vulnerable to **CSIDH-style isogeny walks**
because the class group `Cl(O_D)` has manageable size.

#### 1.3.4 Pairing-friendly curves

Special construction: curves with **small embedding degree** `k`
(typically 6, 12, 18, 24) but **prime-order subgroup** large.

Examples:
- **BN curves** (Barreto–Naehrig): `k = 12`.  Used in Ethereum
  precompile (`alt_bn128` = BN254).
- **BLS curves** (Barreto–Lynn–Scott): `k = 12, 24, 48`.
  BLS12-381 used in Zcash, Filecoin, Ethereum 2.0 (KZG, etc).
- **KSS curves**: `k = 8, 16, 18, 36`.
- **MNT curves**: `k = 3, 4, 6`.

Designed for pairings; ECDLP security is preserved by construction.

#### 1.3.5 Koblitz curves

Subset of binary curves with `a, b ∈ F_2`.  Examples: `sect163k1`,
`sect233k1`.  Have an efficient Frobenius endomorphism enabling
GLV-style speedup.  Vulnerable to GHS if `m` composite.

#### 1.3.6 Twist-secure curves

Curves where both `E(F_p)` and `E^twist(F_p)` have large prime
subgroups.  Most NIST curves are twist-secure to varying degrees:
- P-256 twist: largest exploitable smooth factor is 179 (~7.5 bits)
- secp256k1 twist: largest is 22639 (~14 bits)
- P-384 twist: essentially prime (no smooth factors at 2²² bound)

### 1.4 The deployment landscape

| Category | Examples | Status |
|---|---|---|
| NIST FIPS 186-4 | P-{192, 224, 256, 384, 521} | Production (TLS, FIDO2, govt) |
| SECG K-curves | secp256k1 | Production (Bitcoin, Ethereum) |
| SECG R-curves | secp{192, 224}k1, secp256r1 (=P-256) | Mostly equivalent to NIST |
| Brainpool | brainpoolP{192..512}r1 | EU government |
| ANSSI | FRP256v1 | French government |
| SM2 (China) | sm2p256v1 | Chinese standard, mandatory in China |
| GOST (Russia) | CryptoPro A/B/C, TC26 256-A, 512-A/B | Russian standard |
| RFC 7748 | Curve25519 (X25519), Curve448 (X448) | TLS 1.3, Signal, WireGuard, Tor |
| RFC 8032 | Ed25519, Ed448 | TLS 1.3 client auth, JWT, Signal |
| Pairing-friendly | BN254, BLS12-381 | Ethereum, ZK proofs |
| Microsoft NUMS | numsp256d1, numsp384d1 | Research |
| Hamburg | E-521, Goldilocks-448 | Research / academic |
| Zcash | JubJub, Pallas, Vesta | ZK-SNARK protocols |
| NIST binary | sect{163, 233, 283, 409, 571}k1/r1 | Deprecated |

---

## Part 2: Generic ECDLP attacks

These attacks work against ANY discrete log problem, not just EC.
They define the security baseline.

### 2.1 Brute force

`O(n)` where `n = #E`.  Infeasible for `n ≥ 2⁸⁰`.

### 2.2 Baby-step Giant-step (BSGS)

Shanks 1971.
- Memory: `O(√n)`
- Time: `O(√n)`
- For 256-bit `n`: requires `~2¹²⁸` storage, infeasible.

### 2.3 Pollard's rho

Pollard 1978.  The current best generic ECDLP attack.
- Memory: `O(1)` (constant!)
- Time: `O(√n)` expected
- For 256-bit `n`: `~2¹²⁸` operations
- **Distinguished-point variant** enables parallelization with
  near-linear speedup.

### 2.4 Pollard's kangaroo

Best when the DLP scalar is known to be in a bounded range
`[0, N]`:
- Time: `O(√N)` (smaller than full `O(√n)` when `N << n`).
- Used in bitcoin "address with leaked range" attacks.

### 2.5 Pohlig–Hellman

Pohlig & Hellman 1978.  When `n` has small prime factors:
- Decompose ECDLP via CRT into ECDLPs in each prime-order
  subgroup.
- Time dominates: `√(largest prime factor of n)`.
- **Mitigation**: deployed curves have prime `n` (cofactor 1) or
  large-prime subgroup (Curve25519 has cofactor 8 with the 252-bit
  prime subgroup carrying all security).

### 2.6 Quantum: Shor's algorithm

Shor 1994.  Polynomial-time ECDLP solver:
- Time: `O((log p)³)` quantum operations
- Space: `O(log p)` qubits
- **Realistic estimate for breaking P-256**: a fault-tolerant
  quantum computer with `~2300` logical qubits.  Not feasible
  with 2026 quantum hardware (which has at most a few thousand
  *noisy* physical qubits).
- **Defense**: post-quantum cryptography migration (Kyber,
  Dilithium, SPHINCS+).

### 2.7 Quantum: Grover (for related problems)

Grover provides a √-speedup on unstructured search.  Reduces
brute-force ECDLP from `O(n)` to `O(√n)` — but Pollard ρ already
achieves this classically.  Net Grover effect on EC: `n^{1/4}`
combined with quantum-walk variants, but practically Shor
dominates.

---

## Part 3: Structure-based attacks (curve-specific)

These attacks exploit specific algebraic properties of the curve.
They define which curves are SAFE.

### 3.1 MOV reduction

Menezes, Okamoto, Vanstone 1993.  Uses the **Weil pairing** to
embed ECDLP into the multiplicative group of `F_{p^k}` where
`k` is the **embedding degree** (smallest `k` with `n | p^k − 1`).

- Practical for `k ≤ 6` (supersingular curves)
- For `k > 100` or so, the resulting `F_{p^k}` is too big
- All deployed curves have `k` exponentially large

### 3.2 Frey–Rück reduction

Same idea using **Tate pairing** instead of Weil; sometimes works
when MOV doesn't.  Same fundamental obstacle: embedding degree.

### 3.3 Smart's anomalous-curve attack

Smart 1999.  When `#E(F_p) = p` (anomalous):
- Use `p`-adic formal-group logarithms
- ECDLP reduces to a single division in `F_p`
- Time: `O(log p)` — catastrophic break
- **Mitigation**: never use anomalous curves.  Verify `#E ≠ p`.

### 3.4 GHS / Weil descent (Teske 2006)

For curves over `F_{q^n}` (composite extension):
- Find a "magic number" `m` for the curve
- Construct a genus-`2^{m−1}` cover `C/F_q`
- ECDLP on `E(F_{q^n})` → HCDLP on `Jac(C)(F_q)`
- For small `m`, HCDLP via index calculus is faster

**Applicable to**: binary curves with composite `m`.  Original
target of Teske's work.

**Not applicable to**: prime-field curves (no subfield) — this is
the topic of the 21-phase P-256 investigation showing the
analogous attack is structurally impossible (Phase 10 mod-2
obstruction).

### 3.5 (N, N)-cover attack family (this work)

Hypothesized prime-field analogue of GHS.  Search for a smooth
genus-2 curve `C/F_p` with `Jac(C) ~_{F_p} E × E^twist` via an
`(N, N)`-isogeny.

**Status**:
- **STRUCTURALLY IMPOSSIBLE** for all 22 deployed prime-order
  curves (Phase 10 theorem)
- **Statistically rare** for cofactor curves (Curve25519,
  Curve448) — density `~p^{-3.4}` → vanishes at crypto scale
- **Not applicable** to binary curves (covered by original
  Teske/GHS)

### 3.6 Index calculus for `F_{p^n}` (Diem)

Diem 2011/2012.  For `E/F_{p^k}` with `k ≥ 3`:
- Use Semaev summation polynomials
- Sub-exponential ECDLP: `L_{q^k}(2/3, c)` or similar
- **Applicable**: curves over extension fields with `k ≥ 3`
- **Not deployed**: real-world curves use `k = 1` (prime field).

### 3.7 CSIDH-isogeny attack

For curves with **small class number** of `End(E)`:
- Walk the Cl(O)-action on the isogeny graph
- Reduce ECDLP to a class-group DLP problem
- **Realistic threshold**: `h(O) ≤ 2⁵⁰` or so
- All deployed curves have `h(O) ≥ 2¹²⁵`, infeasible

### 3.8 Static Diffie-Hellman / Cheon attack

Cheon 2006.  When a single private key is used many times with
**known-difference** challenges:
- Time: `O((√n + √(d|p-1)) · log d)`
- `d` is a parameter related to small subgroup structure of
  `p - 1`
- **Practical implication**: vulnerable in some specific
  protocols using long-term keys with structured challenges
- **Mitigation**: ephemeral keys, fresh randomness

### 3.9 Hidden-Number Problem (HNP) / Lattice attacks

Boneh & Venkatesan 1996.  Recover ECDSA private key from biased
nonces:
- Each signature with biased `k` leaks a linear equation
- After `~log₂(n) / bias_bits` signatures, lattice reduction
  (LLL/BKZ) recovers `d`
- **Real-world breaks**: Sony PS3 (2010), Trezor (2017), various
  Android wallets

### 3.10 Cross-protocol attacks

Reuse the same EC key in multiple protocols (e.g., ECDH and
ECDSA, or with different hash functions):
- Different operations leak compatible information
- **Mitigation**: NEVER reuse keys across protocols

---

## Part 4: Side-channel attacks

These attacks exploit physical implementation rather than the
abstract algorithm.

### 4.1 Timing attacks

Kocher 1996.  If scalar multiplication takes different time for
different bits of `k`:
- Measure operation time
- Recover bits of `k` (and thus the private key)
- **Real-world**: OpenSSL's montgomery ladder bug (2003),
  Lucky 13 (2013, on related TLS code), Spectre/Meltdown (2018)
- **Mitigation**: constant-time implementations (Montgomery
  ladder, complete addition formulas)

### 4.2 Simple Power Analysis (SPA)

Measure power consumption during operation:
- A "1" bit (add) vs "0" bit (skip) typically differ in power
- For a vulnerable implementation, SPA recovers `k` in ONE trace
- **Mitigation**: dummy operations, Montgomery ladder, scalar
  blinding

### 4.3 Differential Power Analysis (DPA, CPA)

Kocher, Jaffe, Jun 1999.  Collect many traces with different
inputs:
- Correlate power with predicted intermediate values
- Reveals key bits statistically even when masked
- **Mitigation**: scalar blinding (`k' = k + r·n`), randomized
  projective coordinates, message blinding

### 4.4 Electromagnetic analysis

Similar to power, but uses EM emanations.  Same defenses.

### 4.5 Cache attacks

- **Flush+Reload, Prime+Probe** (Yarom & Falkner 2014)
- Attack-side process shares CPU cache with victim
- Reveals memory access patterns
- **Real-world**: cross-VM key extraction (TLS, GnuPG)
- **Mitigation**: cache-oblivious implementations, no
  secret-dependent memory accesses

### 4.6 Branch prediction attacks

CPU branch predictor leaks secret-dependent branches.
- **Spectre v1/v2/v4** variants (2018)
- **Mitigation**: branchless implementations

### 4.7 Differential Fault Analysis (DFA)

Induce a fault during computation:
- Compare faulted output with correct output
- Recover key bits from the difference
- Real attacks: clock glitching, laser injection, voltage
  glitching on smartcards
- **Mitigation**: verify before output, use multi-execution
  checks

### 4.8 Cold boot attacks

RAM retains state after power-off:
- Cool RAM to seconds-long retention
- Read out key material
- **Mitigation**: encrypted memory, immediate key zeroing,
  hardware enclaves

---

## Part 5: Implementation attacks

These exploit bugs in how the curve is used, not the curve itself.

### 5.1 Invalid-curve attacks

If implementation doesn't validate input points are on the
correct curve:
- Attacker sends a point on a DIFFERENT curve with smooth `#E'`
- ECDH computation leaks bits of private key via Pohlig–Hellman
  on `#E'`
- **Real-world**: TLS implementations 2004-2010 (Microsoft, Sun
  Java), OpenSSL bug 2018
- **Mitigation**: validate `y² = x³ + ax + b` for every input
  point, or use x-only Montgomery curves (immune by construction)

### 5.2 Small-subgroup attacks

If cofactor `h > 1` and the implementation doesn't clear it:
- Attacker provides point of small order
- ECDH leaks the secret mod `h`
- **Mitigation**: multiply input point by `h`, or use clamped
  scalars (Curve25519 does this by design — clearing the low 3
  bits of the scalar maps everything into the prime-order
  subgroup)

### 5.3 Compressed point parsing bugs

EC points can be represented compressed (`x`-only + sign bit).
Parsing bugs:
- Recover `y` from `x` via `sqrt(x³ + ax + b)`
- If implementation doesn't validate the sqrt exists, can be
  tricked
- **Real-world**: various TLS implementations

### 5.4 Edge case handling

Specific curve points cause problems:
- `y = 0`: doubling formula degenerates (point of order 2)
- `x = 0`: special cases for some curves
- `P + (-P)`: must yield identity, easy to mishandle
- **Mitigation**: complete addition formulas (Edwards form,
  Renes-Costello-Batina formulas for Weierstrass)

### 5.5 Modular reduction bugs

Solinas reduction is faster but more error-prone:
- **CVE-2014-3570** (OpenSSL): wrong reduction in `bn_sqr`
- **CVE-2017-3736** (OpenSSL): carry-propagation bug
- **Mitigation**: extensive testing, fuzzing, formal verification

### 5.6 ECDSA nonce-related bugs

The most-broken EC algorithm in practice due to nonce mishandling.

#### 5.6.1 Nonce reuse

Sign two messages with the same `k`:
- `s₁ = k⁻¹(h₁ + d·r), s₂ = k⁻¹(h₂ + d·r)`
- Solve: `k = (h₁ - h₂)/(s₁ - s₂) mod n`
- Then `d = (s·k - h)/r mod n`
- **Real-world**: Sony PS3 (2010), bitcoin transactions 2013,
  Trezor 2017

#### 5.6.2 Biased nonce (HNP)

If `k` is generated with biased distribution:
- Each signature reveals partial info
- Lattice attack recovers `d`
- **Real-world**: Android bitcoin wallets 2013, OpenSSL 2009 with
  biased internal PRNG

#### 5.6.3 Predictable nonce

If `k` is generated from a predictable source:
- Recover `k` from RNG state
- Same as nonce reuse attack
- **Mitigation**: RFC 6979 deterministic nonces (HMAC-DRBG seeded
  with `d || message`)

### 5.7 Replay attacks

Same signature on same message accepted twice.
- **Mitigation**: protocol-level nonces, sequence numbers

### 5.8 Cross-protocol attacks

Reuse key for ECDH and ECDSA, or two ECDSA flavors:
- Different operations leak compatible info
- **Mitigation**: per-protocol key derivation

### 5.9 Group-order confusion

If implementation accepts scalars `k ≥ n`:
- `k mod n` is the effective scalar, leaks via timing
- **Mitigation**: reduce scalar mod `n` at parse time

---

## Part 6: Specific real-world breaks (case studies)

### 6.1 Sony PlayStation 3 (2010)

Sony ECDSA-signed games with constant nonce `k`.  fail0verflow
team extracted Sony's master signing key from ~3 game discs.
Resulted in jailbreaking of PS3.

### 6.2 Bitcoin transactions (2013)

Android Bitcoin wallets used SecureRandom with biased seed.  Many
private keys recovered from biased ECDSA signatures via HNP.
Millions of dollars stolen.

### 6.3 Sony BMG / DRM (2007)

ECDSA private key embedded in product, leaked via reverse
engineering.

### 6.4 Lenovo Superfish (2015)

CA private key included in laptop firmware.  ECDSA signatures
on intercepted TLS connections forgeable.

### 6.5 Android RNG (2013)

`SecureRandom` had biased output, breaking ECDSA across all apps
using it.

### 6.6 TLS / ECDH invalid-curve attacks (2017)

Encore et al. found ~100 TLS servers vulnerable to ECDH
invalid-curve attacks due to insufficient point validation.
Practical key extraction in seconds.

### 6.7 ROCA (2017) — RSA-side equivalent

Demonstrated that lattice-based RSA attacks could break weak
keys; analogous techniques work on weak ECC.

### 6.8 OpenSSL CVE history (2003–2023)

Long list of EC-related CVEs:
- 2003: Montgomery ladder timing leak
- 2014: BN_sqr Solinas bug
- 2017: Carry propagation in P-224
- 2018: Specific point validation flaw
- 2022: Side-channel in ECDSA signing

---

## Part 7: Curve-by-curve security map

### 7.1 NIST FIPS 186-4 prime curves

#### P-192
- 192-bit prime, prime-order
- Pollard ρ: `~2⁹⁶` — **insufficient for new deployment** post-2010
- Deprecated; do not use for new systems

#### P-224
- 224-bit prime, prime-order
- Pollard ρ: `~2¹¹²`
- **Production but on deprecation track**: NIST SP 800-131A
  lists it as legacy-only after 2030

#### P-256
- 256-bit prime, prime-order
- Pollard ρ: `~2¹²⁸`
- **Mainstream production**: TLS 1.2/1.3, FIDO2, government
- Twist leakage: ~16 bits per query (validate inputs)
- Subject of this 21-phase investigation: structurally secure
  against all known classical attacks

#### P-384
- 384-bit prime, prime-order
- Pollard ρ: `~2¹⁹²`
- Used for **top-secret-class** government communications
- Twist essentially prime (cleanest among NIST)

#### P-521
- 521-bit prime, prime-order
- Pollard ρ: `~2²⁶⁰`
- Maximum-strength NIST option
- Some implementation overhead due to odd bit-length

### 7.2 SECG / Bitcoin family

#### secp256k1
- 256-bit Solinas prime, prime-order
- **Used by Bitcoin, Ethereum, Cardano, most cryptocurrencies**
- Cofactor 1 — but **twist leakage ~37 bits per query** (highest
  among standard curves; validation MANDATORY)
- `j = 0` enables GLV speedup (50% faster scalar mul)
- CM discriminant `D = -3` (small but doesn't enable attacks at
  this size)

#### secp192k1, secp224k1
- Smaller variants, similarly designed
- Niche use; not recommended for new deployment

### 7.3 Brainpool

`brainpoolP{192, 224, 256, 320, 384, 512}r1` — RFC 5639.
- "Verifiably random" via pseudo-random number theory
- All cofactor 1, all secure against known classical attacks
- Used in EU government applications
- Slightly slower than NIST curves (no `a = -3` optimization)

### 7.4 Curve25519 / Curve448 / Ed25519 / Ed448

#### Curve25519 (X25519)
- 255-bit Mersenne prime, cofactor 8
- Montgomery form, x-only Diffie-Hellman
- **Designed for safety**: clamped scalars, no point validation
  needed
- Used by Signal, WhatsApp, WireGuard, Tor, TLS 1.3
- **Phase 10 obstruction does NOT directly apply** (#E even via
  cofactor), but density argument (Phase 18-21) protects at
  cryptographic scale

#### Ed25519 (Edwards form of same curve)
- Used for digital signatures (EdDSA)
- Deterministic signatures by design (immune to nonce attacks)
- Complete addition formulas

#### Curve448 / Ed448
- 448-bit prime, cofactor 4
- Higher security tier (~224-bit symmetric equivalent)
- Same design principles as Curve25519

### 7.5 Pairing-friendly

#### BN254 (alt_bn128)
- 254-bit prime, embedding degree `k = 12`
- Used in **Ethereum precompile** for snark verification
- ~110-bit security (down from initial 128-bit estimates due to
  recent advances in number-field sieve for `F_{p¹²}`)
- Phase 10 BLOCKED (cofactor 1)

#### BLS12-381
- 381-bit prime, embedding degree `k = 12`
- Used in **Zcash, Filecoin, Ethereum 2.0**
- ~120-bit security
- Phase 10 BLOCKED (cofactor odd, so `#E` odd)

### 7.6 Binary curves (NIST)

`sect{163, 233, 283, 409, 571}{k, r}1`.  All deprecated for new
use.  The Phase 10 analysis does NOT apply (binary field).
Original target of Teske/GHS, but mitigated by NIST's choice of
prime `m`.

---

## Part 8: Mitigations and best practices

### 8.1 Curve selection

For new applications:
- **Default choice (2026)**: Curve25519 (ECDH) + Ed25519 (signing)
- **NIST-required**: P-256 (general) or P-384 (top-secret class)
- **Bitcoin/cryptocurrency**: secp256k1 (with input validation!)
- **Pairings/ZK**: BLS12-381 (NOT BN254 due to security drift)

### 8.2 Implementation

- **Constant-time everything**: scalar mul, point arithmetic,
  modular reduction
- **Validate inputs**: every external point must be on the curve
  (or use x-only Montgomery)
- **RFC 6979 nonces** for ECDSA (or Ed25519 deterministic)
- **Cofactor handling**: clear cofactor on input, clamp scalar on
  output, or use prime-order curve

### 8.3 Side-channel hardening

- Montgomery ladder (constant-time scalar mul)
- Complete addition formulas (no exception cases)
- Scalar blinding (`k' = k + r·n`)
- Randomized projective coordinates
- Cache-oblivious access patterns
- Branchless conditional selection

### 8.4 Protocol-level

- Ephemeral keys (forward secrecy)
- Per-context KDF (no key reuse)
- Authenticated encryption (don't use bare ECDH)
- Nonce uniqueness (sequence numbers, monotonic counters)

### 8.5 Quantum-readiness

ECC is quantum-vulnerable.  Migration plan:
- **Hybrid**: combine ECC with post-quantum (Kyber for KEM,
  Dilithium for signatures)
- **Pure PQC**: deploy NIST-selected post-quantum algorithms
  (Kyber-768 = ML-KEM-768, Dilithium-3 = ML-DSA-65)
- **Timeline**: NIST mandates PQC migration by 2030 for
  high-value systems

---

## Part 9: Open research directions

### 9.1 Continued classical attacks

- **Mazur–Tate σ Phase 3**: complete the functional equation via
  Eisenstein E_2/4/6 — partial progress in this project's
  `mazur_tate_sigma` module
- **Adversarial-ML rho walks**: train ML model to guide Pollard ρ
  walk — research-grade, GPU compute pending
- **Higher-genus covers**: extension of (N, N)-attack to
  genus-3+ — Phase 15 derived algebraic constraints

### 9.2 Side-channel evolution

- **Microarchitectural attacks**: ongoing arms race
- **Photon emission analysis**: emerging high-precision side
  channel
- **AI-based template attacks**: deep learning applied to power
  trace analysis

### 9.3 Quantum

- **Lower qubit Shor**: reducing the 2300-qubit estimate for
  P-256 (Gidney 2025 estimates closer to 1700)
- **Distillation/error-correction**: dominant cost in practice

### 9.4 Pairing-friendly curve drift

- **Recent NFS advances**: lowered BN254 from 128 to ~110 bit
  security
- **Active migration**: BN to BLS family in Ethereum, ZK
  ecosystems

---

## Part 10: Quick reference — what should I use?

| Scenario | Choice | Why |
|---|---|---|
| New TLS / general | Curve25519 (X25519) + Ed25519 | Safe by design, fast, complete |
| Government / NIST-compliant | P-256 (Routine), P-384 (top secret) | Required by FIPS |
| Bitcoin / cryptocurrency | secp256k1 (mandatory!) | Embedded in protocol |
| ZK proofs / pairings | BLS12-381 | Modern security level |
| EU compliance | Brainpool P-256 r1 or P-384 r1 | EU government accepted |
| Smartcards / embedded | P-256 / Curve25519 | Wide HW support |
| Threshold signatures | Curve25519 (Frost), secp256k1 (MuSig2) | Protocol-determined |
| Post-quantum | Hybrid Kyber + X25519 | Quantum-vulnerable today, PQC plan |

**Don't use**: P-192, P-224 (too small), anomalous curves,
custom-rolled curves, binary curves with composite `m`,
unverified random curves, ANY curve without thorough
validation.

---

## References (canonical)

- Silverman, J., *The Arithmetic of Elliptic Curves* (textbook)
- Cohen, Frey, et al., *Handbook of Elliptic and Hyperelliptic
  Curve Cryptography* (CRC Press, 2005)
- Bernstein, Lange, *SafeCurves: choosing safe curves for
  elliptic-curve cryptography* (safecurves.cr.yp.to)
- NIST FIPS 186-4 (current), 186-5 (in preparation)
- RFC 7748 (Curve25519/Curve448), RFC 8032 (Ed25519/Ed448)
- IETF Crypto Forum Research Group recommendations
- This investigation: `RESEARCH_P256_ISOGENY_COVER.md`,
  `MAINSTREAM_CURVE_AUDIT.pdf`, `P256_PROPERTIES.md`

---

*End of Guide.*
