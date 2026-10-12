---
title: "sect113r1 weak-curve audit: exact pre-admission diagnostic"
author: "crypto research harness"
date: "2026-10-06"
lang: en-US
toc: true
toc-depth: 3
geometry: margin=0.72in
fontsize: 10pt
mainfont: "DejaVu Sans"
monofont: "DejaVu Sans Mono"
colorlinks: true
linkcolor: blue
urlcolor: blue
---

# Status

**PRE-ADMISSION DIAGNOSTIC COMPLETE — NO ADMITTED SCIENTIFIC RUN OR
INDEPENDENT IMPLEMENTATION VALIDATION EXISTS.** A native diagnostic emitted
deterministic exact certificates and replayed them with the same implementation.
Its machine-readable artifact is
[`diagnostics/pre-admission-certificate.json`](diagnostics/pre-admission-certificate.json),
with embedded semantic certificate SHA-256
`1cf4a3733e98c6acdc2219c24f29c15771fc67578bdb33dd426ab643e12572d8`.
The artifact index separately records the raw JSON byte SHA-256. The artifact
explicitly records `admitted_scientific_run=false` and
`evidence_status=PRE_ADMISSION_DIAGNOSTIC_EXACT_CERTIFICATES`.
Its `implemented_diagnostic_gates_passed=true` field means only that the
bounded gates implemented by this diagnostic replayed; it is not an admission,
coverage, exploitability, or independent-validation verdict.

This is a completed **diagnostic**, not a completed experiment. The protocol
and admission requirements remain in [`README.md`](README.md). Exact endpoint
records and ordered routes are recorded in
[`degree5_curve_records.json`](degree5_curve_records.json) and
[`isogeny-routes.json`](isogeny-routes.json). Reproduction commands, artifact
hashes, toolchain details, and validation receipts are indexed in
[`ARTIFACT_INDEX.md`](ARTIFACT_INDEX.md). The editable diagram is
[`evidence-flow.dot`](evidence-flow.dot), with rendered
[`SVG`](evidence-flow.svg) and [`PDF`](evidence-flow.pdf) views.

![Evidence lanes for the exact source-class result, valid degree-5 transfer,
and implementation-input findings. Solid green edges are exact certificate
obligations; dashed amber edges remain modeled or conditional.](evidence-flow.svg){width=95%}

The diagnostic supports one unsurprising class-wide result, two exact transfer
controls, and one conditional implementation result:

- At a frozen 80-bit classical threshold, the exactly certified prime subgroup
  has 56.325748 bits of modeled generic collision work (55.825379 with the
  stated negation optimization). The diagnostic verdict is `CLASS_WEAK`. This is a
  legacy parameter-size classification, independently corroborated by the
  published 2015 `sect113r1` discrete-log computation. It is not a new
  isogeny-induced break.
- Exactly two rational degree-5 kernels, codomains, full-coordinate map
  implementations, order-`n` images,
  forward planted-target homomorphisms, duals, and dual composition on both
  `G` and the planted target replayed. Both routes are
  `TRANSFER_ONLY_SPEEDUP_NOT_ESTABLISHED`; no codomain DLP was solved and no
  representative-specific attack or measured performance advantage was
  supplied.
- The repository's low-level binary point formulas do not use `b` during
  addition/doubling. One actual unchecked multiplier call on the singular
  same-`a`, `b=0` cubic was followed by exact BSGS/CRT recovery of the bounded
  planted scalar. This supports `CONDITIONAL_IMPLEMENTATION_WEAK` only if an
  attacker can
  reach that low-level API with a chosen point and distinguish its full-point
  output. No production protocol exposure was established. The singular cubic
  is not an elliptic curve and is not isogenous to `sect113r1`.

# Evidence vocabulary

- **EXACT DIAGNOSTIC:** deterministic algebraic equality, certificate, or replay
  emitted and checked locally before admission.
- **MODELED:** analytic work estimate; not an executed attack or timing sample.
- **PRIOR ART:** a claim made by a cited external source, not measured here.
- **CONDITIONAL:** depends on an interface exposure not established here.
- **ADMITTED / INDEPENDENTLY VALIDATED:** absent. The verifier replays the
  artifact through the same native implementation and is not an independent
  implementation.

The diagnostic can establish exact equalities within its implementation, but
it cannot promote them to admitted or independently reproduced scientific
evidence. No wall-clock performance result, randomized rho run, prevalence
estimate, or independent verdict is reported.

# Exact source object

The source curve is SECG `sect113r1`:

$$
E/\mathbb F_{2^{113}}:\quad y^2+xy=x^3+ax^2+b,
$$

with polynomial basis modulus `z^113+z^9+1` and SEC 2 coefficients:

```text
a  = 003088250CA6E7C7FE649CE85820F7
b  = 00E8BEE4D3E2260744188BE0E9C723
Gx = 009D73616F35F4AB1407D73562C10F
Gy = 00A52830277958EE84D1315ED31886
n  = 0100000000000000D9CCEC8A39E56F
h  = 2
```

The source is [SEC 2 version
1.0](https://www.secg.org/SEC2-Ver-1.0.pdf), §3.2.1. The frozen PDF SHA-256
is:

```text
d1b16728ad83888fd656d16b99dc71bcd5541d42d848ffd0de7c62c19010d8c3
```

Resolved source representation identity:

```text
icv1-f2m113-tm122610772499221213-97df4ac6
ICV1 (join wrapped lines without spaces):
ICV1:f2m-113-99967757:
-122610772499221213:
10384593717069655379671765157661406:
0x6942e38fc45c62366c09aa8204cd:
unk:unk:r:97df4ac684cb
EC1N113Csect113r1hf529f17bd191
curve UID (join wrapped lines without spaces):
urn:ec-record:1:sha256:
f529f17bd1913792333a661e3557ad6b8e0ca2d4d02b939bc069d17d9fd94d97
```

The diagnostic resolved these source strings. Each valid codomain has a
different representation identity in
[`degree5_curve_records.json`](degree5_curve_records.json), and the ordered
maps are bound in [`isogeny-routes.json`](isogeny-routes.json). The twist has a
different curve identity. The singular control receives no curve identity.

# Exact diagnostic source/class certificates

## Order and generic security

The diagnostic exactly recomputed:

```text
q = 2^113
  = 10384593717069655257060992658440192
n = 5192296858534827689835882578830703
N = #E(F_q) = 2*n
  = 10384593717069655379671765157661406
t = q + 1 - N = -122610772499221213
```

The native diagnostic verified a complete `n-1` primality certificate for `n`,
the generator's exact order, and the unique Hasse-interval multiple supporting
`#E(F_q)=2n`. It also verified irreducibility of `z^113+z^9+1`. These checks
were replayed by the same implementation, not an independent point counter or
validator. The generic references below are formulas, not benchmarks:

| reference | modeled log2 group operations |
|---|---:|
| `sqrt(pi*n/2)` expected collision work | 56.32574806473616 |
| `sqrt(2*n*ln(20))` 95% birthday quantile | 57.29145434846796 |
| frozen policy threshold | 80 |

Thus the diagnostic classification is `CLASS_WEAK`. Because a
rational isogeny preserves the trace and point count, the size classification
applies to every valid representative in this class. A representative can be
worse for a special attack, but it cannot restore more than the class's generic
baseline.

**PRIOR ART.** Wenger and Wolfger state that they computed a full discrete
logarithm on `sect113r1` with ten Kintex-7 FPGAs, with their design processing
900 million iterations per second; see [ePrint
2015/143](https://eprint.iacr.org/2015/143). This report will not rebrand those
published observations as measurements from the present run.

## Frobenius and endomorphism order

Exact diagnostic arithmetic:

```text
D_pi = t^2 - 4*q
     = -26504973335422840129609279124569399
     = -7 * 47 * 411643769 * 195708607903277035369399
```

The diagnostic's exact factor/primality certificates make `D_pi` squarefree
and fundamental. Thus `Z[pi]` is maximal and the diagnostic found no conductor
levels to hunt in the class. The minimum exact norm of a non-integer
endomorphism is

```text
(abs(D_pi)+1)/4 = 6626243333855710032402319781142350
                 = approximately 2^112.35181831539929.
```

This excludes a cheap non-scalar CM/GLV endomorphism under the checked model.
It does not exclude unmodeled algorithms.

## Structural attack ledger

| Attack family | Exact diagnostic result | Evidence boundary |
|---|---|---|
| Pohlig–Hellman | `n` prime with complete `n-1` certificate | exact in diagnostic; no independent replay |
| MOV/Frey–Rück | `ord_n(q)=(n-1)/2` exactly | exact in diagnostic; no independent replay |
| anomalous | `N != q` | exact arithmetic |
| supersingular | ordinary; `t` is odd | exact arithmetic |
| subfield/GLS | 113 is prime; no degree-113 trace from an `F_2` curve equals `t` | exact recurrence in diagnostic |
| GHS/Hess | `ord_113(2)=28` | exact orbit parameter only; magic number, type, genus, and applicability were not computed, so no attack conclusion |
| CM/GLV | fundamental `D_pi`; no low-norm non-scalar element | exact in diagnostic; scoped to checked model |
| direct Semaev/index calculus | no end-to-end result | unavailable; no negative claim |
| cover/Jacobian/extension transfer | generated catalog contains same-field degree-3 `H -> E` cover certificates for the source and both endpoints | cover existence only; subgroup transfer, inverse/recovery, end-to-end cost, and DLP advantage are untested (`dlp_advantage=null`) |

No row may be converted from “not demonstrated” to “impossible” without the
corresponding proof. In particular, equation counts or asymptotic discussion
do not establish an executable index-calculus attack.

# Valid-isogeny diagnostic

Modulo 5, the exact Frobenius polynomial factors as

$$
X^2-tX+q=(X-3)(X-4)\pmod 5.
$$

The diagnostic derived exactly two rational cyclic order-5 kernels from the
fifth division polynomial rather than accepting expected coefficients:

```text
Xq  = x^q mod psi_5
h1  = gcd(psi_5, Xq + x)
Xq2 = Xq^q mod psi_5
h2  = gcd(psi_5, Xq2 + x) / h1
```

Both degree-2 factors passed squarefreeness, divisibility, Frobenius, and
point-order checks. Each route produced a nonsingular Vélu codomain, full-point
map, order-`n` image of `G`, and separately constructed dual. The diagnostic
checked `dual(phi(G))=[5]G`, `dual(phi([d]G))=[5d]G`, and recovered the original
planted target point as `[5^{-1} mod n]dual(phi([d]G))`. This is an exact
homomorphism and point-pullback control; it is not a discrete-log solve on the
codomain. The native certificate recomputes each full ICV1/EC1/UID from the
codomain model and generator. The endpoint and route JSON files are derived
views cross-referenced to that certificate; their byte hashes and the source
revision are bound by the artifact index, not by a separate native bundle
verifier.

The diagnostic evaluated dual composition on `G` and the planted target, not
on the broader protocol set of kernel points, infinity, and sign partners.
Kernel/formula certificates and those two subgroup checks are strong transfer
evidence, but they do not discharge that broader full-morphism validation
obligation for a later admitted run.

Any future speedup claim must still charge the route cost

$$
C^*=C_{\rm construct}+C_{\rm map\ in}+C_{\rm attack}
   +C_{\rm map\ out}+C_{\rm recover}+C_{\rm verify}.
$$

No component may be silently treated as free. The diagnostic did not execute a
matched source/codomain attack or collect timing samples, so it makes no
performance claim.

Exact diagnostic route table:

Endpoint A has `b=0109267245489e254e8f14002629a1`; endpoint B has
`b=0162b1a595685a1387c82647bf44cf`. Their complete identities and kernel
digests are in the linked record and route JSON files.

| route | valid curve/map certificate | subgroup replay | result |
|---|---|---|---|
| source | source order and generator | exact | incumbent; no timing claim |
| endpoint A | forward/dual kernels and formulas | order `n`; target composition and pullback; broader point-set check pending | transfer only; speedup not established |
| endpoint B | forward/dual kernels and formulas | order `n`; target composition and pullback; broader point-set check pending | transfer only; speedup not established |

These are successful transfer controls, not new weak curves. The diagnostic
found no representative-specific advantage in the tested routes, but it did
not exhaust the isogeny class and cannot assert that no such representative
exists.

# Twist safety experiment

The exact diagnostic twist order is

```text
10384593717069655134450220159218980
= 2^2 * 5 * 11 * 17 * 449 * 883 * 1493 * 4690844705882931102817.
```

The diagnostic certified the factorization and primality of the largest factor
`4690844705882931102817`. Its generic collision work of about `2^36.32092` is
a model, not an executed solve. The twist is a different trace class and is
classified `OFF_CLASS_CONTROL_NOT_AN_ISOGENOUS_REPRESENTATIVE`; no interface
exposure was established.

# Singular-input experiment

## Algebraic object and cost

The deliberate control

$$
S:\ y^2+xy=x^3+ax^2
$$

has `b=0` and discriminant zero. It is a singular cubic, not an elliptic curve.
The exact diagnostic found `Tr(a)=1`, making the smooth locus a nonsplit torus
with order

```text
q+1 = 10384593717069655257060992658440193
    = 3 * 227 * 48817 * 636190001 * 491003369344660409.
```

The point `P=(a,0)` passed `[q+1]P=O` and `[(q+1)/r]P != O` for every certified
prime factor `r`. The largest factor has about 58.769 bits; its generic
reference is about `2^29.710` expected work, 26.616 bits below the nominal
subgroup's modeled generic reference.

That comparison is an algebraic work model, not a wall-clock speedup. Any
future end-to-end performance claim must count the projection, each
discrete-log solve, CRT, transcript operations, and final replay.

## Fixed public fixtures

The small-factor product is

```text
M = 3 * 227 * 48817 * 636190001 = 21149740236874377.
```

The bounded proof fixture is:

```text
domain = "crypto/sect113r1-singular/demo/v1\0"
d_demo = 1 + OS2IP_BE(SHA256(domain)) mod (floor(M/2)-1)
d_demo = 6599291615786184 = 0x177205508734c8
```

Native BSGS on all four factors plus CRT recovered `d_demo` exactly from one
actual call to the unchecked multiplier. The recovered residues were:

| factor | residue | BSGS giant steps | exact replay |
|---:|---:|---:|---|
| 3 | 0 | 1 | yes |
| 227 | 58 | 4 | yes |
| 48817 | 26624 | 121 | yes |
| 636190001 | 487863039 | 19342 | yes |

CRT returned `6599291615786184`, and the final nominal public-key replay
matched. The expected scalar was not used by the solver path. The same native
implementation replayed each equation; an independent implementation has not
yet done so.

The full-width fixture is:

```text
seed ASCII = "SECT113R1-SINGULAR-V1"
domain = "crypto/sect113r1-singular-recovery/seed/v1\0"
SHA256(domain || seed) = c4cfd5aa9d48d2079788cbdd3c6c7f579cf7d41cd002106eb2bab540ad639f08
d_full = 1 + OS2IP_BE(SHA256(domain || seed)) mod (n-1)
d_full = 1108386698526968361534713354743651
       = 0x36a5ce943bb0bf6c27642c5b6f63
Qx = 0x00AF667D2D64DF714FCB7229D4509F
Qy = 0x0071C688F4F136AAF006193D1DC24F
Ux = 0x0078AED9A682BE670A024A24D76B1D
Uy = 0x00D567FA101A19AF65CE7F58C0D4A4
```

The public key `Q` and unchecked singular result `U` are post-search replay
controls. They must not be used to prune BSGS or future rho candidates.

For the full-width fixture, one unchecked call produced the frozen output and
the four executed components recovered
`d_full mod M = 13735027495063591`. The remaining prime is
`491003369344660409`; its solve was explicitly **not executed**. The reported
29.710003-bit rho work is modeled. This is residue recovery, not full secret
recovery.

## Implementation boundary and controls

The checked code-level issue is that low-level binary addition, doubling, and
scalar multiplication depend on `a` but do not validate `b` or subgroup
membership. The diagnostic supports `CONDITIONAL_IMPLEMENTATION_WEAK` at this low-level
surface because:

- the actual low-level multiplier accepted `P`;
- its result agreed with a separate torus-law implementation, including an
  exhaustive `GF(2^8)` control;
- the bounded scalar and the declared full-width residue were exactly replayed;
- the nominal equation and checked public-point predicate rejected `P`;
- the predicate also rejected infinity and points outside the order-`n`
  subgroup; and
- caller reachability and a distinguishing production-protocol output remain
  explicit exposure obligations.

No result here licenses “weak isogenous curve.” The strongest result without a
reachable wrapper is “conditional implementation weakness at the tested
low-level surface.”

# Nominal wrong-subgroup control

The on-curve point

```text
T = (0, 0x01DF93AA20579D377F97203F385B11)
```

has exact order 2 in the diagnostic: `[2]T=O`, `[n]T=T`, and the subgroup
predicate rejected it. Raw multiplication distinguished scalar parity. This is
a standard unchecked-surface/wrong-subgroup control on the nominal curve, not
evidence that the prime-order subgroup is easier; production exposure was not
established.

# How confidence will be earned

Exact algebra dominates this study. The diagnostic gains implementation-level
confidence from machine-checkable certificates, mutation tests, a separate
torus law, and final point replay. It does **not** have source-independent
implementation validation: producer and replay verifier share the native code
base. Repeating deterministic arithmetic does not create statistical
confidence if both paths share a bug.

For any performance comparison, the admitted design must pin hardware and
software, use isolated cores, randomize paired execution order, retain warmups
and all failures, include an A/A check, and report every repetition plus a
bootstrap interval on the paired effect. A speedup must use the same operation
boundary and include construction and recovery.

For randomized rho, seeds and censoring budgets are frozen in advance. Success
counts receive exact binomial intervals; time-to-hit uses survival analysis.
The analytic birthday quantiles are not empirical confidence intervals.

For an isogeny walk, adjacent samples are dependent. Without a mixing theorem
or defensible effective sample size, the report will state only the exact
visited set, route coverage, and sample yield. A zero-hit walk is not evidence
that the whole class contains no special representative.

# Diagnostic verdict matrix

| Object / claim | Diagnostic result | Boundary |
|---|---|---|
| prime subgroup at frozen 80-bit threshold | class weak; exact `n`, modeled 56.325748-bit rho work | legacy size; prior full DLP is external |
| two explicit valid degree-5 neighbors | exact kernels, maps, order-`n` images, duals, forward target homomorphism, and point pullback | same-implementation diagnostic; no codomain DLP solve |
| cheaper attack through a valid neighbor | not demonstrated; transfer only | no matched attack or timing study |
| singular same-`a` cubic accepted by unchecked multiplier | bounded scalar exactly recovered; full-width residue recovered | conditional implementation weakness at the tested low-level surface; attacker reachability/output remain obligations |
| full-width singular recovery | not completed | 59-bit component unexecuted; rho cost modeled |
| public-point validation predicate rejects invalid/wrong-subgroup points | exact diagnostic control holds | no checked multiplier or production wrapper integration asserted |
| entire class has no other special attack | no claim | finite degree-5 controls are not exhaustive |

The canonical scoreboard and progress timeline were checked but not updated:
this diagnostic is neither an admitted IC/performance result nor a benchmark
measurement. The registry-driven IC leaderboard and lab-browser data were
regenerated because the two exact endpoint identities expand their curve
roster. The study-local evidence-flow diagram was updated to show exact,
modeled, external-prior-art, and conditional lanes.

# References

1. Standards for Efficient Cryptography Group, [SEC 2 version
   1.0](https://www.secg.org/SEC2-Ver-1.0.pdf), §3.2.1.
2. E. Wenger and P. Wolfger, [“Harder, Better, Faster, Stronger — Elliptic
   Curve Discrete Logarithm Computations on FPGAs,” ePrint
   2015/143](https://eprint.iacr.org/2015/143).
3. P. Gaudry, F. Hess, and N. Smart, “Constructive and destructive facets of
   Weil descent on elliptic curves,” *Journal of Cryptology* 15 (2002).
4. F. Hess, “Generalising the GHS attack on the elliptic curve discrete
   logarithm problem,” *LMS Journal of Computation and Mathematics* 7 (2004).
5. A. Menezes and E. Teske, “Cryptographic implications of Hess' generalized
   GHS attack,” *Applicable Algebra in Engineering, Communication and
   Computing* 16 (2006).
6. J. Vélu, “Isogénies entre courbes elliptiques,” *Comptes rendus de
   l'Académie des sciences* 273 (1971).
7. J. Tate, “Endomorphisms of abelian varieties over finite fields,”
   *Inventiones Mathematicae* 2 (1966).

External sources establish parameters, prior art, and theorem context. They do
not replace the exact run certificates, explicit routes, exposure evidence, or
matched end-to-end costs required above.
