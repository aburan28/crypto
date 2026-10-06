---
title: "sect113r1 weak-curve audit: prospective report"
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

**DRAFT — NO ADMITTED SCIENTIFIC RUN, PRODUCER ARTIFACT, OR INDEPENDENT
VERDICT EXISTS.** This report records the question, exact candidate arithmetic,
and evidence gates before execution. It must not be cited as a completed
experiment. The authoritative prospective contract is
[`README.md`](README.md), and the editable object/claim diagram is
[`evidence-flow.dot`](evidence-flow.dot).

Pre-admission calculations predict one unsurprising class-wide result and one
potential implementation result:

- At a frozen 80-bit classical threshold, the prime subgroup with `log2(n)` just
  above 112 has only about 56 bits of generic collision security. If the order
  certificate passes,
  the entire rational isogeny class is `CLASS_WEAK_LEGACY_SIZE`. This is a
  legacy parameter-size classification, independently corroborated by the
  published 2015 `sect113r1` discrete-log computation. It is not a new
  isogeny-induced break.
- The repository's low-level binary point formulas do not use `b` during
  addition/doubling. An unchecked input on the singular same-`a`, `b=0` cubic
  is therefore a candidate invalid-curve/singular-input oracle. A successful
  planted recovery would be `IMPLEMENTATION_WEAK` at that interface, conditional
  on reachability and observable output. The singular cubic is not an elliptic
  curve and is not isogenous to `sect113r1`.

Two explicit degree-5 neighbors are expected. Their purpose is to prove actual
valid isogeny transport and test whether any representative-specific advantage
survives full accounting. Merely finding a codomain or sharing an order does
not establish such an advantage.

# Evidence vocabulary

- **FROZEN:** fixed in the admitted protocol before producer execution.
- **CANDIDATE DERIVATION:** pre-admission arithmetic to be independently
  recomputed; not yet certified evidence.
- **PRODUCER:** emitted by the admitted native run and bound to its content
  snapshot.
- **VALIDATED:** independently rederived against the snapshot.
- **PRIOR ART:** a claim made by a cited external source, not measured here.
- **OPEN:** required evidence or implementation is unavailable.
- **PENDING:** no admitted result exists.

The words “proved,” “certified,” and “measured” are reserved for artifacts that
pass the admitted producer and independent review. Candidate numbers in this
draft are expectations and falsifiable test vectors.

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
1.0](https://www.secg.org/SEC2-Ver-1.0.pdf), §3.2.1; the frozen PDF SHA-256 is
`d1b16728ad83888fd656d16b99dc71bcd5541d42d848ffd0de7c62c19010d8c3`.

Proposed representation identity:

```text
icv1-f2m113-tm122610772499221213-97df4ac6
ICV1:f2m-113-99967757:-122610772499221213:10384593717069655379671765157661406:0x6942e38fc45c62366c09aa8204cd:unk:unk:r:97df4ac684cb
EC1N113Csect113r1hf529f17bd191
urn:ec-record:1:sha256:f529f17bd1913792333a661e3557ad6b8e0ca2d4d02b939bc069d17d9fd94d97
```

The registry must resolve these strings before admission. Each valid codomain
gets a different representation identity. The twist gets a different curve
identity. The singular control receives no curve identity.

# Pre-admission analytic baseline

## Order and generic security

Candidate exact arithmetic:

```text
q = 2^113
  = 10384593717069655257060992658440192
n = 5192296858534827689835882578830703
N = #E(F_q) = 2*n
  = 10384593717069655379671765157661406
t = q + 1 - N = -122610772499221213
```

`n` is expected prime and must receive a complete primality certificate. The
identity `#E(F_q)=2n` needs an independently replayable point-count certificate,
not only the standard's stated cofactor. If certified, squarefree `2n` and the
elliptic-group invariant-factor constraint imply that `E(F_q)` is cyclic. The
generic references, computed without benchmarking, are:

| reference | candidate log2 group operations |
|---|---:|
| `sqrt(pi*n/2)` expected collision work | 56.32574806473616 |
| `sqrt(2*n*ln(20))` 95% birthday quantile | 57.29145434846796 |
| frozen policy threshold | 80 |

Thus the candidate classification is `CLASS_WEAK_LEGACY_SIZE`. Because a
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

Candidate arithmetic:

```text
D_pi = t^2 - 4*q
     = -26504973335422840129609279124569399
     = -7 * 47 * 411643769 * 195708607903277035369399
```

If all four factors are certified prime, `D_pi` is squarefree and fundamental.
Then `Z[pi]` is maximal, so there are no conductor levels to hunt anywhere in
the class. The minimum candidate norm of a non-integer endomorphism is

```text
(abs(D_pi)+1)/4 = 6626243333855710032402319781142350
                 = approximately 2^112.35181831539929.
```

That would exclude a cheap non-scalar CM/GLV endomorphism. It does not exclude
unmodeled algorithms.

## Structural attack ledger

| Attack family | Candidate exact result | Pending evidence boundary |
|---|---|---|
| Pohlig–Hellman | `n` prime | complete primality/order certificates |
| MOV/Frey–Rück | `ord_n(q)=(n-1)/2=2596148429267413844917941289415351` | exact modular-order proof |
| anomalous | `N != q` | recompute from frozen parameters |
| supersingular | ordinary; `t` is odd | explicit ordinary criterion |
| subfield/GLS | 113 is prime; no degree-113 trace from an `F_2` curve equals `t` | independent recurrence over base traces `-2..2` |
| GHS/Hess | `ord_113(2)=28`; candidate minimum genus `2^28-1` under the stated model | theorem assumptions and independent derivation |
| CM/GLV | no low-norm non-scalar element | fundamental-discriminant and norm certificates |
| direct Semaev/index calculus | no end-to-end result | explicit solver, recovery, and charged reference required |
| cover/Jacobian/extension transfer | no supplied beneficial map | explicit map, subgroup preservation, inverse, and cost required |

No row may be converted from “not demonstrated” to “impossible” without the
corresponding proof. In particular, equation counts or asymptotic discussion
do not establish an executable index-calculus attack.

# Valid-isogeny experiment

Modulo 5, the candidate Frobenius polynomial factors as

$$
X^2-tX+q=(X-3)(X-4)\pmod 5.
$$

The two eigenlines predict two rational cyclic order-5 kernels. The run must
derive them from the fifth division polynomial rather than from expected
coefficients:

```text
Xq  = x^q mod psi_5
h1  = gcd(psi_5, Xq + x)
Xq2 = Xq^q mod psi_5
h2  = gcd(psi_5, Xq2 + x) / h1
```

Both expected degree-2 factors require squarefreeness, divisibility,
Frobenius, and point-order certificates. Each route then requires a
nonsingular Vélu codomain, full-point map, order-`n` image of `G`, separately
constructed dual, and composition `[dual] o [map]=[5]` on a deterministic
test set. The planted logarithm must transport and replay end to end.

The route cost is

$$
C^*=C_{\rm construct}+C_{\rm map\ in}+C_{\rm attack}
   +C_{\rm map\ out}+C_{\rm recover}+C_{\rm verify}.
$$

No component may be silently treated as free. The comparison is to a matched
source-curve attack in the same operation unit. The identity route remains the
incumbent unless a valid routed attack beats it by the preregistered effect
size.

Pending result table:

| route | valid curve/map certificate | subgroup replay | charged advantage |
|---|---|---|---|
| source identity | `PENDING` | `PENDING` | `PENDING` |
| eigenvalue-4 degree-5 neighbor | `PENDING` | `PENDING` | `PENDING` |
| eigenvalue-3 degree-5 neighbor | `PENDING` | `PENDING` | `PENDING` |

An explicit valid neighbor with no cheaper attack is a successful transfer
control, not a failed experiment and not a new weak curve.

# Twist safety experiment

The candidate twist order is

```text
10384593717069655134450220159218980
= 2^2 * 5 * 11 * 17 * 449 * 883 * 1493 * 4690844705882931102817.
```

The largest factor is about `2^71.99034`; its expected generic collision work
is about `2^36.32092`. The twist therefore deserves an exact invalid-input
check. It is a different trace class and cannot be described as a weak
isogenous representative. A security finding additionally requires a reachable
interface that fails to validate the curve or subgroup.

# Singular-input experiment

## Algebraic object and cost

The deliberate control

$$
S:\ y^2+xy=x^3+ax^2
$$

has `b=0` and discriminant zero. It is a singular cubic, not an elliptic curve.
Candidate `Tr(a)=1` makes the smooth locus a nonsplit torus with order

```text
q+1 = 10384593717069655257060992658440193
    = 3 * 227 * 48817 * 636190001 * 491003369344660409.
```

The candidate full-order point is `P=(a,0)`. Its exact-order proof must check
`[q+1]P=O` and `[(q+1)/r]P != O` for every prime factor `r`. The largest factor
has about 58.769 bits; its generic reference is about `2^29.710` expected work,
26.616 bits below the nominal subgroup's expected generic reference.

That comparison is an algebraic work model, not a wall-clock speedup. The run
must count the projection, each discrete-log solve, CRT, transcript operations,
and final replay.

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

Native BSGS on all four factors plus CRT must recover `d_demo` exactly. The
producer may not read the expected scalar while solving. A validator must
independently replay every group equation.

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

Solving only the four small factors supports recovery modulo `M`, not a full
secret recovery. A projected 30-bit solve for the largest factor is a model
until that solve executes under a frozen budget.

## Implementation boundary and controls

The candidate code-level issue is that low-level binary addition, doubling,
and scalar multiplication depend on `a` but do not validate `b` or subgroup
membership. A supported `IMPLEMENTATION_WEAK` finding requires all of:

- the actual low-level multiplier accepts `P`;
- its result agrees with independent singular smooth-locus arithmetic;
- the fixed scalar or declared residue is recovered and exactly replayed;
- the nominal equation rejects `P` before checked multiplication;
- the public-point predicate rejects infinity, off-curve inputs, and points
  outside the order-`n` subgroup; and
- a caller-controlled input and distinguishing output are either demonstrated
  or explicitly left as an exposure obligation.

No result here licenses “weak isogenous curve.” The strongest result without a
reachable wrapper is “conditional implementation weakness at the tested
low-level surface.”

# Nominal wrong-subgroup control

The candidate on-curve point

```text
T = (0, 0x01DF93AA20579D377F97203F385B11)
```

has expected order 2. The run must prove `[2]T=O`, `[n]T=T`, and rejection by
the subgroup predicate. If raw multiplication output is distinguishable, this
point reveals scalar parity. This is a standard subgroup-validation issue on
the nominal curve, not evidence that the prime-order subgroup is easier.

# How confidence will be earned

Exact algebra dominates this study. Confidence comes from machine-checkable
certificates, mutation tests, source-independent recomputation, and final point
replay. Repeating deterministic arithmetic does not create statistical
confidence if both runs share the same bug.

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

# Prospective verdict matrix

| Object / claim | Maximum verdict after its own evidence passes | Current state |
|---|---|---|
| prime subgroup at frozen 80-bit threshold | `CLASS_WEAK_LEGACY_SIZE` | `PENDING` |
| explicit valid degree-5 neighbor | `VALID_ISOGENOUS_REPRESENTATIVE` | `PENDING` |
| cheaper attack through a valid neighbor | `REPRESENTATIVE_WEAK` or `TRANSFER_SPEEDUP` | `PENDING` |
| singular same-`a` cubic accepted by unchecked multiplier | `IMPLEMENTATION_WEAK` | `PENDING` |
| checked entry point rejects invalid/wrong-subgroup points | validation control `HOLDS` | `PENDING` |
| finite valid-neighbor search with no advantage | `NO_WEAKNESS_FOUND_WITHIN_SCOPE` | `PENDING` |
| entire class has no other special attack | no finite-walk verdict allowed | `OPEN` |

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
