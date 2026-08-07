# The Elliptic Curve Discrete Logarithm Problem: State of the Art (2025–2026)

*Compiled June 2026. A survey of the best known classical algorithms, special-case
attacks, quantum resource estimates, record computations, and the implications for
standardized curves and the post-quantum migration.*

---

## TL;DR

- For a well-chosen curve of group order `n`, the **best known classical attack is still
  generic square-root** — Pollard rho / parallel collision search / kangaroo at
  `~sqrt(n)` ≈ `2^(n_bits/2)` group operations. This has not changed.
- **No subexponential classical algorithm exists for prime-field ECDLP.** Index calculus
  (summation polynomials, point decomposition, Gröbner bases) yields subexponential or
  improved attacks **only over extension fields** `F_{q^n}`, and the most optimistic
  binary-field claims rest on the **first-fall-degree assumption, now widely doubted**.
- **Special-case breaks** (anomalous curves, low-embedding-degree / MOV, Weil descent /
  GHS) and **implementation breaks** (invalid-curve, small-subgroup, twist, biased-nonce
  lattice/Bleichenbacher attacks) do **not** apply to correctly-implemented P-256,
  Curve25519, or secp256k1.
- The **practical record frontier** is ~**112-bit** (prime field) / ~**117-bit** (binary
  field) for genuine ECDLP, and ~**130-bit interval** DLP for the Bitcoin-puzzle kangaroo
  solves — roughly `2^58`–`2^65` operations, about `10^19` times short of breaking a
  256-bit curve.
- **Shor's algorithm breaks ECDLP in polynomial time** on a fault-tolerant quantum
  computer, and **ECC-256 is an easier quantum target than RSA-2048**. Latest logical-qubit
  estimates for ECDLP-256 are ~**1,200 logical qubits** with tens of millions of Toffoli
  gates — but **no cryptographically relevant quantum computer (CRQC) exists**, and the
  hardware gap remains 1–2 orders of magnitude in logical qubits.
- **NIST IR 8547** (deprecate ECC by 2030, disallow by 2035) and **NSA CNSA 2.0** drive the
  migration; **hybrid X25519MLKEM768** is already the dominant PQC key exchange on the web
  (~38% of human HTTPS traffic on Cloudflare by early 2025).

---

## 1. Classical generic algorithms (the standing security bound)

The security of every standardized curve rests on the fact that, with no exploitable
structure, the only attack is generic.

- **`O(sqrt(n))` is the best known and essentially optimal.** Pollard rho, baby-step
  giant-step, and van Oorschot–Wiener parallel collision search all cost
  `~sqrt(pi*n/2)` group operations. Shoup's 1997 generic-group lower bound — any
  `m`-query generic algorithm succeeds with probability `< (m+2)^2/(2p) + 1/p` — matches
  the BSGS upper bound, so generic algorithms are optimal *in the generic model*. Applying
  that bound to a concrete curve is a modeling assumption, but it is the universally-used
  effective security bound.
  *(Galbraith–Gaudry survey, DCC 2015, doi 10.1007/s10623-015-0146-7; Shoup, EUROCRYPT 1997.)*

- **Parallelization is linear.** Van Oorschot–Wiener replaces one long rho walk with many
  short walks to *distinguished points*, collisions detected via a shared list; expected
  work divides cleanly across processors. This is what all record computations use.
  *(van Oorschot & Wiener, J. Cryptology 1999.)*

- **The negation map gives a `sqrt(2)` (~1.41×) speedup** by walking on classes `{P, -P}`,
  but induces *fruitless cycles* that trap the walk; realized speedups are
  implementation-dependent (reported ~1.29× up to near-2× with careful SIMD cycle
  handling), generally short of the theoretical `sqrt(2)`.
  *(Bernstein–Lange–Schwabe 2011; Bos–Kleinjung–Lenstra, Information Sciences 2012.)*

- **Interval ECDLP** (logarithm known to lie in an interval of width `w`) is solved by
  **Pollard's kangaroo (lambda)** in `~O(sqrt(w))`, parallelizable linearly. Interleaved
  BSGS and Gaudry–Schost variants can edge out plain kangaroo. This is the regime of the
  Bitcoin-puzzle solves (§4).

**Bottom line:** as of 2025–2026 there is **no accepted classical algorithm beating
`O(sqrt(n))` for standard prime-field curves**; recurring "sub-sqrt" preprints have not
survived scrutiny.

### 1a. Precomputation and amortization — beating `sqrt(n)` *online*

Rho's `sqrt(n)` is the bound for a *single* discrete log in an *unknown* group. Two
distinct techniques do better when the group is **fixed** (e.g. a standardized curve) and
either precomputation is amortized or only online cost is counted.

- **Hellman-table / cube-root preprocessing (the `S·T² = N` tradeoff).** Spend a one-time
  precomputation on the fixed group, store an `S`-bit advice string, then answer each query
  in online time `T`, with `S·T² = Θ̃(N)`. At the balanced point `S = T = N^(1/3)`:
  **online time `N^(1/3)`, advice `N^(1/3)`, precomputation `≈ N^(2/3)`.** Algorithm:
  Mihalcik (2010), Bernstein–Lange, *Computing small discrete logarithms faster*
  (LATINCRYPT 2012). **Optimality:** Corrigan-Gibbs & Kogan, *The Discrete-Logarithm Problem
  with Preprocessing* (EUROCRYPT 2018, eprint 2017/1113) proved `S·T² = Ω̃(N)` in the generic
  group model — so `N^(1/3)` online is the best possible with `N^(1/3)` advice. A 2022 ITC
  paper made the table construction fully explicit/optimal. This is the substance of
  Bernstein–Lange, *Non-uniform cracks in the concrete* (ASIACRYPT 2013, eprint 2012/318):
  in the **non-uniform** model P-256's advertised `2^128` is only `~2^85` online. *Caveat:*
  the precomputation for P-256 is `~N^(2/3) ≈ 2^170` — far more than `2^128` — so this is a
  statement about the security *model* and about small-`N`/many-target regimes, **not a
  practical break**. For a one-off log in a random group, rho's `sqrt(n)` still wins.
- **Batch / amortized rho (`sqrt(LN)` for `L` targets).** Kuhn–Struik, *Random Walks
  Revisited* (SAC 2001), share a distinguished-point database across `L` targets in the same
  group: total cost `~sqrt(2LN)`, i.e. `~sqrt(N/L)` amortized per log — a **`sqrt(L)`
  speedup** over independent rho. No expensive separate precomputation; the win is real and
  practical when many keys live on one curve, bounded by the number of targets.

Both interact with **interval ECDLP**: the cube-root preprocessing takes interval-DLP of
width `w` from `sqrt(w)` to `w^(1/3)` online (Bernstein–Lange's "small logarithms" target),
so a precomputed table for a fixed curve speeds the Bitcoin-puzzle-style kangaroo solves
too — again only if the `w^(2/3)` precompute is affordable. None of this dents 256-bit
curves in practice.

---

## 2. Index calculus on elliptic curves (works only over extension fields)

The major research thread that *could* have broken ECDLP — and the reason the prime-field
case is reassuring precisely because the thread failed there.

- **Semaev summation polynomials (2004).** `S_{m+1}` vanishes exactly when `m+1` curve
  points with the given x-coordinates sum to zero. This is the algebraic engine of every
  elliptic index-calculus attempt. *(Semaev, ePrint 2004/031.)*

- **Gaudry / Diem point decomposition.** Over `F_{q^n}`, Weil descent + Gröbner-basis
  solving of the summation system yields, for fixed `n`, heuristic time roughly
  `O(q^{2 - 2/n})`, beating generic `q^{n/2}` once `n` is large enough.
  *(Gaudry, ePrint 2010/157; Diem.)*

- **Diem's subexponential result.** For suitable families over `F_{q^n}` as `n` grows
  (including some binary fields heuristically), genuinely **subexponential** complexity is
  achievable — the strongest theoretical index-calculus success against elliptic curves.

- **Petit–Quisquater (ASIACRYPT 2012)** conjectured subexponential **binary-field** ECDLP
  *under the first-fall-degree assumption*.

- **The refutation.** The **first-fall-degree assumption is now widely doubted**:
  computational evidence and constructed counterexamples (the CRYPTO 2015 "last fall
  degree" work, ePrint 2015/573; ePrint 2015/358 on generalized FFD assumptions) show it
  fails in general and is unlikely to hold for the relevant Weil-descent systems. This is
  the main reason the 2012-era subexponential optimism for binary ECDLP has receded.

- **Prime fields are untouched.** No subexponential algorithm is known for ECDLP over
  `F_p` (or binary fields of *prime* extension degree). Lifting / "xedni calculus"
  approaches fail to produce small enough relations (Miller's argument, strengthened by
  Silverman–Suzuki). Even where index calculus formally applies, generic `O(sqrt(n))` still
  wins at practically-used parameter sizes — the advantage is purely asymptotic and
  assumption-dependent. *(Galbraith–Gaudry survey; ePrint 2017/609.)*

---

## 3. Special-case and implementation attacks (avoidable by construction)

None of these are generic ECDLP breaks; they target curves or implementations with
exploitable structure.

### Structural (curve-class) attacks
| Attack | Target class | Effect | Hits standard curves? |
|---|---|---|---|
| **MOV / Frey–Rück** | low embedding degree `k` | transfers ECDLP to DLP in `F_{p^k}^*`, then subexponential index calculus | **No** — P-256/secp256k1/25519 have huge `k` |
| **Anomalous (SSSA: Smart / Semaev / Satoh–Araki, 1997–99)** | trace 1, `#E = p` | **polynomial-time** via p-adic elliptic log into `(F_p, +)` | **No** — none are anomalous |
| **Weil descent / GHS (2001–02)** | binary `F_{2^n}`, composite `n` | maps to higher-genus curve where index calculus is easier | **No** — prime-field / well-chosen binary curves |

- **exTNFS fallout (Kim–Barbulescu, CRYPTO 2016)** improved finite-field DLP in
  `F_{p^n}` with composite `n`, which **re-scored pairing-friendly curves**: BN256, once
  believed 128-bit, was re-estimated at **~100-bit** security. The community moved to
  larger curves (**BLS12-381**, **BN462**) for ~128-bit pairings. *(This affects
  pairing-based crypto, not the base ECDLP of standard signing/KEX curves.)*

### Implementation attacks (recover keys *without* solving ECDLP)
- **Invalid-curve attacks**: unvalidated input points on a weaker curve leak the private
  key mod small subgroup orders (CRT); demonstrated against TLS-ECDH.
- **Small-subgroup attacks**: exploit cofactor > 1 via low-order points; harmless for
  prime-order (cofactor 1) curves.
- **Twist attacks**: x-only/Montgomery-ladder ECDH with an x-coordinate on the quadratic
  twist leaks bits if the twist has small factors. **secp256k1 and P-256 are not strongly
  twist-secure** and require point validation for x-only ECDH; **Curve25519 was designed
  twist-secure.** A 2026 advisory (GHSA-r6ph-v2qm-q3c2) shows missing SECT-curve subgroup
  validation in pyca/cryptography — these flaws remain live.
- **Biased-nonce ECDSA** (the most practically important): partial nonce leakage reduces
  to the Hidden Number Problem, solved by lattice reduction (LLL/BKZ) or Bleichenbacher
  FFT — **no ECDLP solve required**. **Minerva (2020)** recovered 256-bit keys from
  ~500–2,100 signatures across libgcrypt/wolfSSL/MatrixSSL/SunEC/Crypto++; **LadderLeak
  (CCS 2020)** needed <1 bit of leakage; **2024 work (Osaki–Kunihiro, SAC 2024)** tolerates
  high nonce-error rates via 4-list sum algorithms.

**Takeaway:** secp256k1, P-256, and Curve25519 have prime/near-prime order, huge embedding
degree, and no anomalous/descent structure → best generic attack is `2^128` rho. Their real
residual risk is **implementation** (point validation, nonce bias), not the math.

### 3a. Can isogenies map a prime-field curve to a weaker one? — No (a structural argument)

A natural attack idea: an isogeny `φ: E → E'` is a group homomorphism, so a small-degree
isogeny transports a DLP instance from `E` to `E'` and back — if `E'` were weaker, you'd
win. It fails for well-chosen prime-field curves, and the reason is a **provable invariance**:

- **Tate's isogeny theorem:** two curves over `F_p` are isogenous **iff** `#E(F_p) =
  #E'(F_p)`. So the group order `n`, and hence the trace `t = p+1−n`, is **constant across
  the entire isogeny class.**
- Therefore every structural-attack precondition is a **class invariant**:
  - the large prime subgroup order `r | n` (governs rho `sqrt(r)`) — fixed;
  - **anomalous** status (`t = 1`) — fixed: not anomalous ⇒ nothing isogenous is;
  - **embedding degree** `k` (= order of `p` mod `r`) — fixed, since `p` and `r` are.
- The **only** thing that varies within a class is the **endomorphism ring** (Kohel's
  isogeny volcano). But **no known ECDLP attack exploits the endomorphism ring**; the closest
  thing — efficiently computable endomorphisms (GLV/GLS, e.g. secp256k1's `λ`) — gives only
  a **constant-factor** rho speedup, never an asymptotic one.

So the set of curves reachable from a good prime-field curve by isogeny is exactly the set
sharing its (good) group order — none weak to any known attack.

The isogeny-walk-to-a-weak-curve idea *does* work in **one** place: **GHS Weil-descent over
composite-degree binary fields** `F_{2^n}`, where Hess / Maurer–Menezes–Teske showed you can
hop via a small isogeny to a GHS-vulnerable curve (e.g. the `GF(2^155)` results). That needs
a subfield to descend to and is **absent over prime fields**. Related confusions: the
**SIDH/SIKE break (Castryck–Decru 2022)** recovers an isogeny from torsion-point images (the
opposite direction) and says nothing about ECDLP; **quantum isogeny-path** algorithms
(Childs–Jao–Soukharev) solve isogeny-finding, not ECDLP, where Shor is already polynomial.

**Assessment:** a prime-field isogeny break would require a *new* curve invariant not
determined by the group order, varying within a class (the endomorphism ring is the only
candidate) and admitting a faster discrete log — something the entire GLV/GLS/CM literature
has not produced. This route looks markedly less promising for prime fields than the
(already stalled) index-calculus and (extension-field-only) Weil-descent routes.

---

## 4. Record computations (and what they prove)

> **Crucial distinction:** *full random-instance ECDLP records* (~112–117 bit) are not
> comparable to *interval/kangaroo "puzzle" solves* (up to ~130-bit interval), which only
> work because the key is range-confined **and** the public key was exposed.

### Full ECDLP records
- **Certicom challenges:** ECCp-109 solved Nov 2002 (~10,000 machines, Monico); ECC2-109
  solved Apr 2004 (~2,600 machines). **ECC2K-130 (131-bit class) is the smallest unsolved
  Certicom challenge** — a 12+ institution CPU/GPU/Cell/FPGA effort built the attack but
  never completed it; all 131-bit-and-larger challenges remain open.
- **112-bit prime field (secp112r1):** solved July 2009 by Bos–Kaihara–Kleinjung–
  Lenstra–Montgomery on 200+ PlayStation 3 consoles over ~6 months — **largest classic
  prime-field secp record.**
- **114-bit prime field (Barreto–Naehrig curve w/ automorphisms):** ~2017/18, ~2,000 CPU
  cores over ~6 months — largest *prime-field* ECDLP, on a special curve.
- **117.35-bit binary field (`F_{2^127}`):** Bernstein–Engels–Lange–Niederhagen–Paar–
  Schwabe–Zimmermann, FPGA, 2016 (ePrint 2016/382) — **largest completed binary-field
  ECDL.**
- **No full random-curve ECDLP above ~118 bits has ever been completed.**

### Interval ECDLP — Bitcoin "puzzle" kangaroo solves (secp256k1)
These are `~2^(N/2)`-op interval DLPs, enabled by exposed public keys:
- **#115** (114-bit interval): 16 Jun 2020, Zieniewicz & Pons, JeanLucPons Kangaroo,
  256× Tesla V100, ~13 days (~`2^58` ops).
- **#120, #125**: solved 2023–2024 (kangaroo).
- **#130** (~`2^130` interval, ~`2^65` ops): solved **Sep 2024** (RetiredCoder /
  RCKangaroo — "SOTA" symmetry method, `K≈1.15` vs `2.1` classic, ~8 GKeys/s on RTX 4090).
- **#135** (~`2^67.5` ops): **unsolved as of late 2025**, targeted by distributed pools.
- Lower puzzles *without* exposed pubkeys are brute-forced (BitCrack): #66 (Sep 2024),
  #67 (Feb 2025), #68 (Apr 2025), #69 (Apr 2025) — several "stolen" by RBF/mempool bots
  re-deriving the key from the briefly exposed pubkey at spend time.

### Security-margin implication
Largest completed work ≈ `2^58`–`2^65` operations. A 256-bit curve needs `~2^128`
operations by rho/kangaroo — roughly `2^63` (≈ `10^19`) times harder. **Classically
infeasible.** The puzzle solves do **not** threaten secp256k1 in normal use: a standard
unused address (a *hash* of the pubkey, full 256-bit range) offers no shortcut.

---

## 5. Quantum attacks (the actual threat) and resource estimates

- **Shor solves ECDLP in polynomial time** via period-finding on the group; for generic
  prime-order groups no quantum algorithm asymptotically beats it. **Grover is irrelevant**
  to ECDLP (only `O(sqrt(N))` unstructured search; it matters for symmetric/hash sizes).

- **ECC-256 is an easier quantum target than RSA-2048.** Because the best *classical* ECC
  attack (rho, `sqrt`) is far costlier per bit than GNFS on RSA, ECC keys are much smaller,
  and Shor's cost scales with operand size — so P-256 needs **~2.6× fewer logical qubits
  and ~100–150× fewer gates** than RSA-2048/3072.

**Trajectory of 256-bit ECC Toffoli-count estimates (illustrative, one curve):**

| Year | Work | Logical qubits | Toffoli gates | Architecture/notes |
|---|---|---|---|---|
| 2017 | Roetteler–Naehrig–Svore–Lauter | ~2,330 | ~1.26 × 10¹¹ | seminal concrete estimate; `9n + O(log n)` qubits |
| 2020 | Häner–Jaques–Naehrig–Roetteler–Soeken | improved trade-offs | lower depth/T-count | windowed arithmetic, Q# implementation |
| 2023 | Litinski (PsiQuantum) | — | **~44–50 × 10⁶** | "active-volume" photonic; ~1 key / 10 min on ~6,000 modules |
| 2023 | Gouzien et al. (Alice & Bob) | — | — | **126,133 cat qubits, ~9 hours** |
| 2026* | Chevignard–Fouque–Schrottenloher (EUROCRYPT 2026) | **~1,193** (≈ 3.12n) | ~2³⁸ (quartic) | Legendre-symbol compression, avoids modular inversion; space↓ gate↑ |
| 2026* | Google QAI / Ethereum Fdn / Stanford (secp256k1) | **~1,200–1,450** | **~70–90 × 10⁶** | <500k physical qubits, runtime in minutes (superconducting) |

\* 2026-dated items rest partly on secondary reporting + arXiv preprints (some details of
the Google ECDLP paper reportedly disclosed only via a ZK proof) — **treat exact figures as
provisional** until proceedings are confirmed.

- **Hardware gap is still large.** Google Willow (Dec 2024) showed below-threshold error
  correction with 105 physical qubits; public logical-qubit counts are in the *tens*
  (Quantinuum/Microsoft ~12, Atom/Microsoft ~24). IBM's roadmap targets ~200 logical qubits
  (~10,000 physical) around 2028–29. Breaking ECC-256 needs **~1,200+ logical qubits /
  hundreds of thousands of physical qubits** — a 1–2 order-of-magnitude gap.
- **Expert timeline (GRI / Mosca Quantum Threat Timeline 2025):** a CRQC is "quite possible"
  (28–49%) within 10 years and "likely" (51–70%) within 15; ~92% of surveyed experts put
  ≥50% probability at 20 years. Timeline viewed as accelerating.

---

## 6. Standardized curves & the post-quantum migration

- **Classical status:** P-256/P-384, Curve25519/Ed25519, secp256k1 have **no practical
  classical break**; ~128-bit (P-384: ~192-bit) security stands. secp256k1's GLV
  endomorphism is only a ~1-bit constant-factor speedup. The **NIST P-curve seeds remain
  unexplained** (fails SafeCurves "rigidity," a lingering post-Dual_EC trust debate), but
  **no hidden weakness has ever been shown**; Curve25519 satisfies SafeCurves fully.
- **Consensus:** ECDLP is **not** expected to fall to classical math advances — the
  migration is driven **exclusively by Shor**. Koblitz–Menezes counsel humility given the
  history of special-instance surprises, but no cryptographer is predicting a classical
  break.
- **The SIDH/SIKE break (Castryck–Decru 2022)** exploited torsion-point info via Kani's
  criterion — it killed an *isogeny* KEM and has **no bearing on ECDLP** or
  P-256/25519/secp256k1.
- **NIST IR 8547** (initial public draft Nov 2024): RSA/ECDSA/EdDSA/ECDH at 112-bit
  **deprecated after 2030, disallowed after 2035**; ≥128-bit ECC also disallowed after 2035
  for federal use. *(As of mid-2026 it still appears to be at IPD status — secondary reports
  conflict; verify the final at csrc.nist.gov/pubs/ir/8547.)*
- **NSA CNSA 2.0** (2022, updated 2025): direct migration to PQC (no hybrid required for
  NSS); networking-gear exclusivity ~2030, OS/apps/cloud ~2033, all NSS quantum-resistant by
  2035.
- **Harvest-now-decrypt-later** is the stated reason to migrate **ECDH key exchange now**:
  recorded traffic is retroactively decryptable once a CRQC exists.
- **PQC replacements:** FIPS 203 (ML-KEM), 204 (ML-DSA), 205 (SLH-DSA) finalized Aug 2024;
  **HQC** selected as fifth algorithm (code-based KEM backup) Mar 2025. **Hybrid
  X25519MLKEM768** is the dominant deployed TLS 1.3 KEX (Chrome default since v124, Apr
  2024; OpenSSL/Go/Apple OSes). By early 2025 Cloudflare measured **~38% of human HTTPS
  traffic** using hybrid PQC key agreement, though origin-side support lagged (~4–10%).

---

## 7. Active research directions (2025–2026)

1. **Quantum resource minimization for ECDLP-256** — the live frontier: driving logical
   qubits toward ~1,200 and Toffoli counts toward tens of millions (Litinski active-volume;
   Chevignard et al. Legendre-symbol compression; the Google/Ethereum secp256k1
   spacetime-volume optimization). Benchmarking proposals like *"Brace for impact: ECDLP
   challenges for quantum cryptanalysis"* (arXiv 2508.14011, 2025) define challenge ladders
   to track the quantum threat to Bitcoin.
2. **Extension-field index calculus & the last/first-fall-degree question** — refining
   complexity bounds, probing whether any salvageable subexponential binary-field attack
   survives the collapse of the first-fall-degree assumption. **See §8 for the detailed
   frontier** (symmetrized Semaev, Gröbner bottleneck, cover/gluing attacks, EDS equivalence,
   and the prime-field "Holy Grail").
3. **Better generic-attack engineering** — kangaroo symmetry methods (RCKangaroo SOTA),
   negation-map cycle handling, GPU/FPGA throughput; relevant to records and interval DLP.
4. **Nonce-leakage / HNP attacks** — increasingly powerful Bleichenbacher-FFT and
   lattice methods tolerating tiny or noisy bias (the practically dangerous class).
5. **PQC migration engineering** — hybrid KEX deployment, agility, key-transparency, and
   the 2030/2035/2033 compliance timelines.

---

## 8. The algebraic-attack frontier in detail (mapped to this repo's modules)

This section goes deeper on the index-calculus / descent program — the directions this
repository actually pursues — and states, per direction, the current literature consensus and
the concrete obstruction. **Net: every one of these attacks is confined to extension fields
or stalls on a Gröbner/degree bottleneck; none beats Pollard rho at any cryptographic size,
and the prime-field case (`secp256k1`, P-256) has no working index calculus at all.**

### 8.1 Symmetrized summation-polynomial index calculus over `F_{q^n}`
*(repo: `symmetrized_semaev.rs`, `semaev_higher.rs`, `groebner_f4.rs`)*

- **The symmetrization gain is real but bounded.** Faugère–Gaudry–Huot–Renault (J. Cryptology
  2014) get an **exponential `2^{ω(n−1)}` factor** on the point-decomposition step for twisted
  Edwards / Jacobi-intersection models, via the `S_n` action and elementary-symmetric /
  power-sum variable changes (plus small-torsion tricks for more compact polynomials). This is
  exactly the density reduction the repo's `decompose_symmetric` targets.
- **But the Gröbner step is the wall, and it's a brutal constant.** The dominant cost is F4/F5
  solving the Semaev system, *not* the linear algebra. Galbraith–Gebregiyorgis (INDOCRYPT
  2014) measured `<1000` decomposition systems/sec versus `25,000–150,000` rho point-additions/sec,
  and concluded **"Pollard rho is still much faster than index calculus for ECDLP over
  `F_{2^n}` of reasonable size."** This is the empirical reality the repo's `groebner_f4`
  performance envelope (≈6 vars before the Buchberger brake) runs into; a research-grade
  F4/F5 + FGLM + sparse linear algebra is what would even make the comparison interesting.
- **Asymptotics.** Gaudry/Diem give point-decomposition `Õ(q^{2−2/n})` for fixed `n` (with the
  Gaudry–Thomé–Thériault–Diem double-large-prime variation), beating generic `q^{n/2}` only
  once `n` is large — and even then the "n treated as constant" framing hides the Gröbner cost
  that prevents a concrete win at small fixed `n = 2,3,4,5`.
- **Incremental 2012–2018 progress, no class change:** Joux–Vitse cover-and-decomposition
  (composite `n`, see §8.3); Vitse–Wallet *Improved Sieving on Algebraic Curves* (LATINCRYPT
  2015, ~3× constant speedup); Amadori–Pintore–Sala one-relation variant (§8.5). **Consensus:
  symmetrized index calculus threatens no standardized curve.**

### 8.2 The first-fall-degree question — is binary ECDLP subexponential?
*(repo: `RESEARCH_FFD_*`, `semaev_sat.rs`)*

- Petit–Quisquater (ASIACRYPT 2012) claimed `O(2^{c·n^{2/3} log n})` for `F_{2^n}` **conditional
  on the first-fall-degree assumption** (degree of regularity ≈ first fall degree). Even granting
  it, the **crossover with rho is non-cryptographic** — their own estimate beats rho only for
  `n ≳ 2000`, vs deployed `n ≈ 160–600`.
- **The assumption is now widely distrusted.** Kosters–Yeo (*Notes on summation polynomials*,
  arXiv 1503.08001, 2015) give experimental evidence the degree of regularity **grows with `n`**
  (a rise observed around `n ≈ 45` for `m=2`), contradicting `D_reg ≈ D_FirstFall`.
  Huang–Kosters–Yeo (*Last fall degree, HFE, and Weil descent attacks on ECDLP*, CRYPTO 2015,
  eprint 2015/573) introduce the monomial-order-independent **last fall degree**, prove an
  *unconditional* poly-time solve for HFE, and **construct Weil-descent systems where the
  first-fall-degree assumption is unlikely to hold.** Huang–Petit–Shinohara–Takagi (eprint
  2015/358) try to *salvage* a generalized assumption — the live tension is "likely false"
  (Kosters–Yeo) vs "refine it" (Petit et al.).
- **Concrete bounds that *do* hold:** Kousidis–Wiemers (J. Math. Cryptol. 2019) prove the first
  fall degree of the Semaev Weil-descent system is `≤ m² − m + 1` (improving `m² + 1`); Huang
  (arXiv 2103.07282, 2021) bounds the last fall degree of *linearized* descent systems — but
  the **Semaev/ECDLP case remains open.** The repo's FFD harness is measuring exactly the
  invariant whose growth is the crux.
- **A provable negative:** the *naive* Semaev index calculus **cannot beat Pollard rho** over
  an arbitrary finite field (J. Math. Cryptol. jmc-2019-0029, 2020) — and experimentally not
  even brute force.
- **A closed-off side route:** quasi-subfield polynomials (Huang–Kosters–Yeo 2020) were killed
  by Euler–Petit's non-existence theorem (FFA 2021, arXiv 1909.11326): the needed parameter
  regime can't be met, so **no speedup.**

### 8.3 Cover / descent / genus-2 gluing — and why prime fields are immune
*(repo: `RESEARCH_MESTRE_HOWE.md`, Richelot / Howe-gluing autolab work, `secp256k1_cm`)*

- **All cover attacks need a subfield.** GHS and its generalizations (Hess; Galbraith–Hess–Smart
  isogeny extension; Diem's odd-characteristic GHS) transfer ECDLP over `F_{q^n}` to a
  higher-genus Jacobian DLP over a subfield. Diem showed the odd-char GHS **fails for prime
  extension degree ≥ 11** (genus blows up); only `d = 2,3,5,7` need analysis. Joux–Vitse
  cover-and-decomposition (EUROCRYPT 2012) broke a real **156-bit curve over `F_{p^6}`** — but
  it is a *composite-degree* attack.
- **Newer genus-2/3 covers for prime-*order* curves still need an extension field.** Song Tian,
  *Cover Attacks for Elliptic Curves over Cubic Extension Fields* (J. Cryptology 2023, arXiv
  2012.07173): an `F_q`-rational `(ℓ,ℓ,ℓ)`-isogeny from the Weil restriction to a **genus-3**
  Jacobian solves DLP in `Õ(q)` for *some* prime-order curves over `F_{q^3}`. Fan (2022)
  constructs genus-2 covers over *quadratic* extensions. **Neither touches prime fields.**
- **The gluing machinery the repo is exploring is constructive, not an attack vector for
  `F_p`.** Howe–Leprévost–Poonen `(2,2)`-gluing of two elliptic curves into a genus-2 Jacobian,
  and Richelot `(2,2)`-isogenies between genus-2 Jacobians, are used to *build* curves and in
  *isogeny-based* PQC; the Castryck–Decru SIDH break uses them **with exchanged torsion-point
  images** — side information plain ECDLP doesn't expose. Scholten's construction
  (`E/F_{p²} → ` genus-2 curve over `F_p` via Weil restriction) and "weak Weil restriction"
  curves give an obstruction **for `E/F_{p²}`, not for prime-field curves.** Concretely, for
  the repo's secp256k1 sextic-twist gluing: `secp256k1` lives over a *prime* `F_p`, so there is
  **no Weil restriction to descend through**, and any cover correspondence preserves the group
  order (Tate) while the index-calculus cost in a genus-`g` Jacobian is `Õ(q^{2−2/g})` — needing
  small `g`, whereas prime-order covers force large (non-hyperelliptic) genus. That is the
  **dual obstruction**: (i) no subfield, (ii) genus growth makes any cover costlier than rho.

### 8.4 Elliptic divisibility sequences / elliptic nets
*(repo: `RESEARCH_EDS_RESIDUE.md`)*

- Shipsey's and Swart's EDS attacks solve ECDLP **only when the point order equals `q−1`** —
  i.e. the same multiplicative structure MOV/Frey–Rück already exploit.
- Lauter–Stange (SAC 2008, eprint 2008/099) define **EDS Association, EDS Residue, and EDS
  Discrete Log** and prove each is subexponential **iff ECDLP is** — a rigorous *equivalence*,
  not an advantage. The repo's `EDS_RESIDUE` work is therefore probing a problem provably no
  easier than the target; index calculus via elliptic nets is in turn **equivalent to Semaev's
  summation polynomials**. **Consensus: no advantage over generic `sqrt(n)` for prime fields.**

### 8.5 Direct prime-field index calculus — the "Holy Grail"
*(repo: `diem_descent.rs` and its `Amadori–Pintore–Sala` reference, `pkm_criterion.rs`)*

- The structural problem is stark: over `F_p` there is **no subfield** to define a small factor
  base, so Semaev's polynomials have nowhere to land. Petit–Kosters–Messeng (PKC 2016) were the
  **first to make summation-polynomial IC run at all over prime fields**, via factor bases from
  compositions of small-degree rational maps — but **only for tiny parameters**, with asymptotics
  bounded by the heuristic Gröbner behavior. Amadori–Pintore–Sala (FFA 2018, eprint 2017/609)
  cut the relations needed to **one** (a Las Vegas algorithm, success ≈ 0.6) and outperform
  *other IC variants* in the prime-field case — **but not Pollard rho.** Kudo–Yokota (CANS 2018,
  title literally *"…and Its Limitation"*) accelerate it and document that it still doesn't beat
  generic methods.
- **Refuted/closed historically:** Silverman's **xedni calculus** was decisively refuted
  (Jacobson–Koblitz–Silverman–Stein–Teske, DCC 2000) — an absolute bound on lifted-relation
  coefficient size makes it asymptotically almost-certain to fail. Semaev's 2015 subexponential
  claim (eprint 2015/310) was **binary-only and rests on the now-distrusted first-fall-degree
  assumption** (§8.2), not a prime-field result.
- **Status of the Holy Grail:** *open but unmoved.* No subexponential prime-field algorithm
  exists and none is refuted as impossible; decades of attempts have produced only
  small-parameter demonstrators. This is the honest framing for the repo's prime-field
  ambitions — the existence question is genuinely open, but there is no published toehold.

**What would actually constitute progress** (and what the repo's bets ride on): a Gröbner/
solving-degree breakthrough that makes Semaev decomposition sub-rho at small fixed `n`; a proof
(or refutation) of a first-/last-fall-degree bound for the *Semaev* system; a prime-field
factor base that doesn't need a subfield; or a genuinely new curve invariant (cf. §3a) tying
ECDLP to an easier problem. Absent one of these, `secp256k1`/P-256 retain full `2^128` classical
security and the operative threat remains Shor.

---

## Sources (selected, by section)

**Classical/generic:** Galbraith–Gaudry survey (DCC 2015, doi 10.1007/s10623-015-0146-7);
Shoup EUROCRYPT 1997; van Oorschot–Wiener (J. Cryptology 1999); Bernstein–Lange–Schwabe
negation map (cr.yp.to/elliptic/negation-20110102.pdf); Bos–Kleinjung–Lenstra (Inf. Sci.
2012).
**Precomputation/amortization:** Corrigan-Gibbs–Kogan EUROCRYPT 2018 (eprint 2017/1113);
Bernstein–Lange *Non-uniform cracks in the concrete* ASIACRYPT 2013 (eprint 2012/318);
Bernstein–Lange *Computing small discrete logarithms faster* LATINCRYPT 2012; Mihalcik 2010;
fully-constructive optimal preprocessing (ITC 2022); Kuhn–Struik *Random Walks Revisited*
SAC 2001.
**Isogenies & ECDLP:** Tate isogeny theorem; Kohel isogeny volcanoes (thesis 1996);
Galbraith–Hess–Smart / Maurer–Menezes–Teske isogeny-extension of GHS; Castryck–Decru
(eprint 2022/975, distinct problem); Childs–Jao–Soukharev quantum isogeny path (2014).
**Index calculus:** Semaev ePrint 2004/031; Gaudry ePrint 2010/157; Petit–Quisquater
EUROCRYPT/ASIACRYPT 2012; last-fall-degree ePrint 2015/573; ePrint 2015/358; prime-field
ePrint 2017/609.
**Algebraic-attack frontier (§8):** Faugère–Gaudry–Huot–Renault (J. Cryptology 2014, HAL
hal-00700555); symmetrized/small-torsion (doi 10.1007/978-3-642-55220-5_3);
Galbraith–Gebregiyorgis INDOCRYPT 2014 (eprint 2014/806); Gaudry–Thomé–Thériault–Diem double
large prime; Kosters–Yeo arXiv 1503.08001; Huang–Kosters–Yeo CRYPTO 2015 (eprint 2015/573);
Kousidis–Wiemers (J. Math. Cryptol. 2019, arXiv 1906.05594); Huang arXiv 2103.07282; naive-IC
lower bound (jmc-2019-0029); Euler–Petit quasi-subfield (FFA 2021, arXiv 1909.11326);
Joux–Vitse EUROCRYPT 2012 (eprint 2011/020); Vitse–Wallet LATINCRYPT 2015; Tian (J. Cryptology
2023, arXiv 2012.07173); Fan (Math. Probl. Eng. 2022); Howe–Leprévost–Poonen gluing; Scholten
construction; Lauter–Stange SAC 2008 (eprint 2008/099); Petit–Kosters–Messeng PKC 2016;
Amadori–Pintore–Sala FFA 2018; Kudo–Yokota CANS 2018; xedni refutation (Jacobson–Koblitz–
Silverman–Stein–Teske, DCC 2000).
**Special-case/impl:** SafeCurves (safecurves.cr.yp.to transfer/twist); Smart/SSSA;
Gaudry–Hess–Smart ePrint 2001/084; Kim–Barbulescu exTNFS (CRYPTO 2016); Minerva
(minerva.crocs.fi.muni.cz); LadderLeak ePrint 2020/615; Osaki–Kunihiro SAC 2024;
GHSA-r6ph-v2qm-q3c2.
**Records:** Bos et al. 112-bit (joppebos.com/presentations/112bitECDLP.pdf); 114-bit BN
(Springer 10.1007/978-3-319-78556-1_13); binary-field ePrint 2016/382; ECC2K-130 ePrint
2009/466 & ecc-challenge.info; JeanLucPons/Kangaroo; RetiredC/RCKangaroo;
privatekeys.pw/puzzles.
**Quantum:** Roetteler et al. arXiv 1706.06752; Häner et al. ePrint 2020/077; Litinski
arXiv 2306.08585; Gouzien et al. arXiv 2302.06639; Gidney 2025 (RSA) arXiv 2505.15917;
Chevignard–Fouque–Schrottenloher (EUROCRYPT 2026); arXiv 2508.14011; GRI Quantum Threat
Timeline 2025.
**Standards/PQC:** NIST IR 8547 (csrc.nist.gov/pubs/ir/8547/ipd); NSA CNSA 2.0; FIPS
203/204/205; NIST HQC selection (Mar 2025); Cloudflare "State of the post-quantum Internet
in 2025"; Castryck–Decru ePrint 2022/975; Koblitz–Menezes ePrint 2015/1018.

> **Methodology / caveats.** Synthesized from nine parallel literature searches (five
> broad-survey passes for §§1–7, four targeted deep-dives for §8: symmetrized index calculus,
> cover/gluing attacks, first-/last-fall-degree, and EDS + prime-field IC). Many primary PDFs
> (eprint, arXiv, NIST, Springer, Cloudflare, Wikipedia) returned HTTP 403 to automated fetch,
> so some figures rest on authoritative search snippets + corroborating secondary sources;
> these are flagged inline. The least-settled points are: (a) exact post-exTNFS bit levels for
> BLS12-381/BN462; (b) the **final** status of NIST IR 8547 (still IPD as far as could be
> confirmed); (c) all **2026-dated quantum estimates**, pending published proceedings; (d) a
> few §8 snippet-only figures (Joux–Vitse exact bit-record, Amadori–Pintore–Sala largest
> instance, the "8th symmetrized polynomial" detail). The core qualitative claims — generic
> `sqrt(n)` classical hardness, no prime-field subexponential attack, the first-fall-degree
> assumption being distrusted, cover/gluing attacks confined to extension fields, Shor as the
> operative threat, ~112–117-bit record frontier — are high-confidence and corroborated across
> independent sources.
