# sect113r1 weak-curve audit protocol (prospective draft)

Status: **DRAFT — NO SCIENTIFIC RUN HAS BEEN ADMITTED OR EXECUTED**

Date drafted: 2026-10-06

This directory specifies a reproducible audit of `sect113r1`, its rational
degree-5 isogeny neighbors, its quadratic twist, and two deliberately invalid
or wrong-subgroup inputs.  Numbers marked **candidate derivation** below came
from pre-admission arithmetic used to design the test.  They are not producer
results.  The native producer and an independent validator must rederive them
from the frozen parameters before the report can promote them to certified
facts.

The audit asks three different questions and never substitutes one for
another:

1. **Class-wide legacy strength.**  Is the prime-order subgroup below a frozen
   classical security threshold?  Every valid curve in the same
   `F_(2^113)`-isogeny class has the same point count, so this classification is
   class-wide.
2. **Valid-representative weakness.**  Does an explicitly constructed,
   nonsingular isogenous representative enable an attack cheaper than the
   source curve after charging construction, maps, subgroup transport, and
   recovery?
3. **Implementation-input weakness.**  Does the low-level multiplication
   surface accept a point from a singular same-`a` formula-compatible cubic or
   the nominal curve's order-2 subgroup?  Such a finding concerns validation
   at that interface.  A singular cubic is not an elliptic curve and cannot be
   called an isogenous weak curve.

The editable evidence-flow diagram is
[`evidence-flow.dot`](evidence-flow.dot).  [`REPORT.md`](REPORT.md) is a
prospective report shell; its result cells must remain pending until admitted
artifacts exist.

## 1. Admission gate

The experiment is not admitted until all of the following are true:

- this protocol, the producer, tests, exact expected schema, independent
  validator, and review assignment are committed and pushed;
- the commit SHA, dirty-state check, toolchain versions, CPU identity, resource
  limits, and SHA-256 of every frozen input are recorded;
- the standard parameter source and its digest are pinned;
- a repository experiment ID, decision ID, run ID, lane, and archive target
  have been allocated by the repository workflow rather than invented in this
  document;
- the producer emits machine-readable evidence and never uses the candidate
  values below as an acceptance oracle;
- the validator receives only the frozen inputs and producer artifact, then
  independently recomputes every scientific joint; and
- failures, timeouts, OOMs, missing receipts, and mismatches remain visible and
  resolve to `INDETERMINATE` or `BREAKS`, never to a strength claim.

Until that gate is satisfied, local executions are diagnostics only.  Their
outputs must not be pasted into the result table as scientific evidence.

## 2. Frozen object identities

### 2.1 Valid source curve

The standard source is SECG `sect113r1` over the polynomial-basis field

```text
F = GF(2)[z] / (z^113 + z^9 + 1)
a = 0x003088250CA6E7C7FE649CE85820F7
b = 0x00E8BEE4D3E2260744188BE0E9C723
Gx = 0x009D73616F35F4AB1407D73562C10F
Gy = 0x00A52830277958EE84D1315ED31886
n = 0x0100000000000000D9CCEC8A39E56F
h = 2
```

The model is

```text
E: y^2 + x*y = x^3 + a*x^2 + b.
```

The parameter source is SECG, [SEC 2 version
1.0](https://www.secg.org/SEC2-Ver-1.0.pdf), printed pages 24–25 (physical PDF
pages 30–31).  The source PDF must be pinned as SHA-256
`d1b16728ad83888fd656d16b99dc71bcd5541d42d848ffd0de7c62c19010d8c3`
by the admitted archive.

The proposed repository identity is:

```text
slug: icv1-f2m113-tm122610772499221213-97df4ac6
ICV1: ICV1:f2m-113-99967757:-122610772499221213:10384593717069655379671765157661406:0x6942e38fc45c62366c09aa8204cd:unk:unk:r:97df4ac684cb
EC1 alias: EC1N113Csect113r1hf529f17bd191
curve UID: urn:ec-record:1:sha256:f529f17bd1913792333a661e3557ad6b8e0ca2d4d02b939bc069d17d9fd94d97
field SHA-256: da55718b2ae51e38fc5d836fcf62bbca5907b877c61b29d5e3c30a5e94b60ee3
```

The exact ICV1 model-hash preimage is sorted-key, compact UTF-8 JSON:

```json
{"a":"0x3088250ca6e7c7fe649ce85820f7","b":"0xe8bee4d3e2260744188be0e9c723","field":"f2m-113-99967757","form":"y^2+xy=x^3+a*x^2+b","modulus":"0x20000000000000000000000000201","v":"1"}
```

Its candidate SHA-256 is
`97df4ac684cbe5e6ab7686805d46ec4aa74945314512fd4ed79598337d243117`.
The admitted run must resolve this identity through the committed registry.
The identity is representation-specific; it is not proof of any mathematical
property.

### 2.2 Objects that must keep separate identities

- Each valid degree-5 codomain receives its own ICV1, EC1 alias, and curve UID
  after its coefficients and a route certificate are emitted.
- The quadratic twist has trace `-t`, a different order, and a distinct curve
  identity.  It is not in the source curve's rational isogeny class.
- The `b=0` same-`a` cubic is singular.  It receives a typed
  `singular-model/v1` artifact identity only—never ICV1, EC1, or a curve UID.
- The order-2 point is on the nominal curve, but it is not in the subgroup
  generated by `G`; its probe identity must bind both the point and intended
  subgroup.

## 3. Candidate exact derivations to re-prove

These values freeze what the producer and validator must check; they are not a
record that the checks have passed.

### 3.1 Source class

Let `q=2^113`.  Candidate arithmetic gives

```text
q = 10384593717069655257060992658440192
#E(F_q) = N = 2*n
          = 10384593717069655379671765157661406
t = q + 1 - N
  = -122610772499221213
D_pi = t^2 - 4*q
     = -26504973335422840129609279124569399
     = -7 * 47 * 411643769 * 195708607903277035369399
n - 1 = 2 * 7 * 53 * 547 * 2848799 * 8757107 * 512797404440011
```

Required certificates:

1. Rabin irreducibility certificate for `z^113+z^9+1` over `GF(2)`.
2. An independently replayable point-count certificate for `#E(F_q)=2n`, not
   only reliance on the standard's stated cofactor.
3. Exact curve-equation and `[n]G=O` checks, plus `[n/r]G != O` for every
   distinct certified prime factor `r` if the producer claims exact order.
4. A complete, replayable primality certificate for `n`; a probable-prime flag
   is insufficient.
5. Complete primality certificates for every claimed prime factor of `D_pi`.
6. Recompute `t`, `D_pi`, and all factorizations from the frozen parameters.

Because the candidate `N=2n` is squarefree, the standard invariant-factor
constraint `d_1^2 | N` predicts that `E(F_q)` is cyclic.  That conclusion also
depends on the point-count and primality certificates.

If all factors of `D_pi` are certified prime, then `D_pi` is squarefree and
`D_pi = 1 (mod 4)`.  The candidate conclusion is that it is fundamental,
`Z[pi]` is already maximal, and every ordinary curve in the rational isogeny
class has the same maximal endomorphism order.  The smallest candidate norm of
a non-integer element is

```text
(abs(D_pi) + 1) / 4 = 6626243333855710032402319781142350
log2 = 112.35181831539929...
```

This would rule out a low-degree non-scalar GLV/CM endomorphism anywhere in the
class; it does not prove that every conceivable attack is expensive.

### 3.2 Generic attack and frozen strength threshold

Freeze `80` classical bits as the policy threshold before execution.  The
nominal prime subgroup has `log2(n)` just above 112 (and a 113-bit integer
encoding).  The analytic generic references are

```text
sqrt(pi*n/2): log2 work = 56.32574806473616...
sqrt(2*n*ln(20)): log2 work = 57.29145434846796...  (95% collision quantile)
```

The second expression is a birthday-model quantile, not a measured solver
confidence interval.  An implementation using negation has a smaller constant;
the audit records, but does not disguise, that distinction.

At the frozen 80-bit threshold, a certified `n` is sufficient for the analytic
verdict `CLASS_WEAK_LEGACY_SIZE`.  This is not a newly discovered special
structure: Wenger and Wolfger report the historical full `sect113r1` discrete
log computation using ten Kintex-7 FPGAs in [ePrint
2015/143](https://eprint.iacr.org/2015/143).  The audit must report that prior
result as external corroboration, not as its own measurement.

### 3.3 Other structural routes

The producer and validator must independently evaluate the following exact
gates:

| Route | Candidate obligation | Claim allowed after certificate |
|---|---|---|
| Pohlig–Hellman on `G` | certify `n` prime | no subgroup-order decomposition gain |
| MOV/Frey–Rück | prove `ord_n(q)=(n-1)/2=2596148429267413844917941289415351` | no low embedding degree |
| anomalous | check `N != q` | anomalous attack not applicable |
| supersingular | check `t` odd and ordinary criteria | supersingular transfer not applicable |
| subfield/GLS | prove the only proper subfield is `F_2`; enumerate the five possible degree-113 traces from base traces `-2..2` and show none equals `t` | no `F_2`-defined representative or GLS orbit gain |
| GHS/Hess | prove `ord_113(2)=28`; derive the prime-degree model's minimum candidate genus `2^28-1=268435455` with its exact theorem/model assumptions | standard GHS route is not competitive |
| low-norm CM/GLV | certify fundamental `D_pi` and the norm bound above | no low-degree non-scalar endomorphism in the class |
| Semaev/direct index calculus | provide an end-to-end charged solver or label unavailable | no conclusion from equation generation alone |
| field extension, cover, or Jacobian transfer | explicit homomorphisms, subgroup preservation, inverse/recovery, and whole-path cost | no speedup from representation change alone |

The five candidate traces for an `F_2`-defined curve after degree-113 base
extension are

```text
-2^57, -1267584991505179, 0, 1267584991505179, 2^57.
```

The admitted report must say that the GHS value is a derivation under the cited
Menezes–Teske model, not a theorem quoted for this exact parameter outside the
paper's displayed range.

### 3.4 Twist

Candidate arithmetic gives

```text
#E_twist(F_q) = q + 1 + t
              = 10384593717069655134450220159218980
              = 2^2 * 5 * 11 * 17 * 449 * 883 * 1493
                * 4690844705882931102817
```

The largest candidate factor has `log2=71.99033773270085...`; its generic
reference is `log2(sqrt(pi*r/2))=36.320916931086586...`.  Complete factor and
primality certificates are required.  A twist result is a twist-safety or input
validation finding, never a same-isogeny-class result.

## 4. Valid degree-5 isogeny experiment

The Frobenius characteristic polynomial modulo 5 is expected to be

```text
X^2 - t*X + q = (X - 3)(X - 4) mod 5.
```

It has two distinct eigenlines, predicting exactly two rational cyclic
degree-5 kernels.  This prediction is not an edge certificate.  The producer
must construct the fifth division polynomial `psi_5` (degree 12) and derive
kernels without using precomputed expected coefficients:

```text
Xq  = x^q mod psi_5                         (113 modular squarings)
h1  = gcd(psi_5, Xq + x)                    (expected degree 2)
Xq2 = Xq^q mod psi_5                        (113 more squarings)
g2  = gcd(psi_5, Xq2 + x)                   (expected degree 4)
h2  = g2 / h1                               (expected exact degree 2)
```

For each `h_i`, require all of:

- monic, degree 2, squarefree, exact divisor of `psi_5`, and invariant under
  the claimed Frobenius action;
- kernel nonzero points have exact order 5;
- the characteristic-2 Vélu codomain is nonsingular and satisfies
  `b' = b + s + s^2` for the independently derived kernel trace `s`;
- explicit full-point source-to-codomain map, not only an x-map;
- image of `G` lies on the codomain and has order `n`;
- a separately constructed dual route composes to `[5]` on a deterministic
  point set, including `G`, a planted subgroup point, kernel points, infinity,
  and sign partners;
- route artifacts bind ordered source/destination curve UIDs, degree,
  direction, kernel polynomial, map formula/version, and SHA-256;
- a planted discrete-log instance transported to the codomain and back keeps
  the same scalar; and
- all construction, mapping, recovery, and verification costs are charged in
  one declared unit before any speedup claim.

The expected outcome is a transfer control: valid neighbors preserve the
already-small prime subgroup and do not by themselves make the logarithm
easier.  If explicit routes pass, the allowed statement is “two valid
degree-5 neighbors were certified and the subgroup was transported.”  It is
not “an isogeny broke sect113r1.”

## 5. Invalid-input and wrong-subgroup experiment

### 5.1 Singular same-`a` companion

Freeze the formula-compatible cubic

```text
S: y^2 + x*y = x^3 + a*x^2       (b=0).
```

For this binary model, `b=0` makes the discriminant zero.  `S` is singular and
not an elliptic curve.  Candidate arithmetic gives `Tr(a)=1`, so its smooth
locus is a nonsplit torus of order

```text
q + 1 = 10384593717069655257060992658440193
      = 3 * 227 * 48817 * 636190001 * 491003369344660409.
```

Use `P=(a,0)`.  Required order certificates are `[q+1]P=O` and
`[(q+1)/r]P != O` for every distinct factor `r`, with complete factor
primality certificates.  The largest factor's generic reference is

```text
log2(sqrt(pi*r/2)) = 29.710003333569247...
95% birthday quantile = 30.675709617301045...
```

The nominal subgroup reference is about 26.615744731166913 bits more work.
That is an operation-model ratio, not a timing measurement.

The fixed small-factor product is

```text
M = 3 * 227 * 48817 * 636190001
  = 21149740236874377.
```

For a bounded deterministic proof-of-capability, use the public fixture

```text
domain = "crypto/sect113r1-singular/demo/v1\0"
SHA256(domain) = fa4be280c9e9f1e2112ae090729ed3b4faf78e0c1133db4a71836e266896b463
d_demo = 1 + OS2IP_BE(SHA256(domain)) mod (floor(M/2) - 1)
d_demo = 6599291615786184
hex(d_demo) = 0x177205508734c8
```

The producer must derive the scalar from that frozen domain and mapping,
query only the actual low-level multiplier, solve each small factor by native
BSGS, combine with CRT, and recover `d_demo` exactly because `d_demo < M/2`.
The validator recomputes the transcript and every BSGS equation independently.

A second full-width public fixture may establish only a residue:

```text
seed ASCII = "SECT113R1-SINGULAR-V1"
domain = "crypto/sect113r1-singular-recovery/seed/v1\0"
SHA256(domain || seed) = c4cfd5aa9d48d2079788cbdd3c6c7f579cf7d41cd002106eb2bab540ad639f08
d_full = 1 + OS2IP_BE(SHA256(domain || seed)) mod (n - 1)
d_full = 1108386698526968361534713354743651
hex(d_full) = 0x36a5ce943bb0bf6c27642c5b6f63
Qx = 0x00AF667D2D64DF714FCB7229D4509F
Qy = 0x0071C688F4F136AAF006193D1DC24F
Ux = 0x0078AED9A682BE670A024A24D76B1D
Uy = 0x00D567FA101A19AF65CE7F58C0D4A4
```

`Q=[d_full]G` is the legitimate public-key replay target and
`U=[d_full]P` is the single unchecked-call transcript. These frozen values are
post-search controls only; neither BSGS nor a future rho walk may consult them
to select or prune candidates.

Unless the large factor is actually solved within a preregistered budget, the
allowed claim is recovery modulo `M`, not full scalar recovery.  A modeled rho
cost is not a completed attack.

The implementation claim requires evidence that the nominal low-level
multiplier produces the same result as the singular formula on `P`, and that
the checked public-point entry point rejects `P` before multiplication.  It is
conditional on an attacker being able to supply a point to that low-level
surface and observe a distinguishing result.  No protocol exposure is assumed.

### 5.2 Nominal order-2 control

The candidate nominal-curve point

```text
T = (0, 0x01DF93AA20579D377F97203F385B11)
```

must satisfy `[2]T=O` and `[n]T=T`, and the checked subgroup predicate must
reject it.  If a raw multiplication result is distinguishable, `[d]T` leaks
the parity of `d`; that is a wrong-subgroup validation result, not a weakness
in the prime subgroup.

## 6. Controls and falsification

Required controls:

- SEC 2 source bytes and parameter round trip;
- field irreducibility and independent field arithmetic implementation;
- nominal generator acceptance and exact subgroup order;
- singular probe rejected as off-curve by the nominal equation;
- order-2 point accepted on-curve but rejected by `[n]P=O` subgroup checking;
- every valid isogeny image accepted by the codomain and subgroup predicate;
- kernel points map to infinity; non-kernel points do not;
- dual composition, sign, infinity, and mutated-kernel negative controls;
- mutation of one bit in every transcript, certificate, and route artifact
  makes validation fail;
- no expected scalar, residue, kernel coefficient, or codomain coefficient is
  consulted by the producer's search path; and
- re-execution by the validator from the frozen public inputs.

Falsification rules:

- A failed exact certificate falsifies the dependent claim.
- A valid neighbor without an end-to-end cheaper attack does not establish a
  representative weakness.
- A singular-input recovery with no demonstrated reachable low-level surface
  supports only a conditional implementation finding.
- Failure to find a weakness in a finite walk supports only
  `NO_WEAKNESS_FOUND_WITHIN_SCOPE`; it cannot certify the entire class strong.
- A timeout, missing artifact, or validator disagreement is `INDETERMINATE` or
  `BREAKS`, not a negative result.

## 7. Measurement and statistical confidence

Most decisive statements in this audit are exact.  They need certificates and
independent reproduction, not p-values:

- field irreducibility;
- point, curve, order, and primality checks;
- trace, discriminant, factorization, embedding degree, and endomorphism norm;
- kernel/codomain/map/dual-composition equations; and
- planted-scalar or residue replay.

Keep four measurement types separate:

1. **Algebraic cost.** Count field multiplications, squarings, inversions,
   additions, point operations, table entries, and bytes.  State conversion
   rules; do not add unlike units.
2. **Wall time and memory.** If compared, use the repository isolated runner,
   pinned cores, frequency/governor metadata, warmups, randomized interleaving,
   A/A control, and a preregistered timeout.  Report every repetition, median,
   dispersion, and bootstrap interval.  These measurements describe that
   implementation and host only.
3. **Randomized rho success.** Freeze independent seeds and budgets before
   running.  Report censored failures.  Estimate success probability with an
   exact binomial interval and time-to-hit with survival methods; never discard
   failed seeds.  The birthday formulas above remain analytic references.
4. **Graph prevalence.** An exhaustive certified finite set may use an exact
   count.  A walk sample is dependent: without a mixing bound or defensible
   effective sample size, report sample yield and path coverage only—no
   binomial prevalence interval.  Zero hits in a bounded walk is not a
   class-wide absence certificate.

The primary “weak” assertion uses the frozen 80-bit threshold and certified
subgroup order; it needs no timing hypothesis test.  Any claim that one valid
representative is faster to attack than another requires matched end-to-end
costs, a preregistered effect size, and uncertainty on the paired difference.

## 8. Result schema (initial state)

Every admitted producer artifact must populate these rows without deleting
failures:

| Claim | Required evidence | Initial state |
|---|---|---|
| legacy class below 80-bit threshold | certified `n`; frozen generic formula | `PENDING_ADMITTED_RUN` |
| nominal order and CM facts | exact certificates and independent replay | `PENDING_ADMITTED_RUN` |
| two rational 5-isogenies | kernels, codomains, full maps, duals, subgroup replay | `PENDING_ADMITTED_RUN` |
| special valid representative weakness | charged end-to-end attack below matched source reference | `PENDING_ADMITTED_RUN` |
| singular small-factor scalar recovery | actual unchecked surface, exact BSGS/CRT replay, controls | `PENDING_ADMITTED_RUN` |
| full-width singular recovery | actual solve and final replay, not projected cost | `PENDING_ADMITTED_RUN` |
| checked-path rejection | nominal curve and subgroup predicate controls | `PENDING_ADMITTED_RUN` |
| twist safety | exact twist order/factors and explicit interface relevance | `PENDING_ADMITTED_RUN` |

## 9. Primary references

- Standards for Efficient Cryptography Group, [SEC 2 version
  1.0](https://www.secg.org/SEC2-Ver-1.0.pdf), §3.2.1 (`sect113r1`).
- E. Wenger and P. Wolfger, [“Harder, Better, Faster, Stronger — Elliptic
  Curve Discrete Logarithm Computations on FPGAs,” ePrint
  2015/143](https://eprint.iacr.org/2015/143).
- P. Gaudry, F. Hess, and N. Smart, “Constructive and destructive facets of
  Weil descent on elliptic curves,” *Journal of Cryptology* 15 (2002).
- F. Hess, “Generalising the GHS attack on the elliptic curve discrete
  logarithm problem,” *LMS Journal of Computation and Mathematics* 7 (2004).
- A. Menezes and E. Teske, “Cryptographic implications of Hess' generalized
  GHS attack,” *Applicable Algebra in Engineering, Communication and
  Computing* 16 (2006).
- J. Vélu, “Isogénies entre courbes elliptiques,” *Comptes rendus de
  l'Académie des sciences* 273 (1971).
- J. Tate, “Endomorphisms of abelian varieties over finite fields,”
  *Inventiones Mathematicae* 2 (1966), for the finite-field isogeny criterion.

These sources establish parameters, prior art, and theorem context.  They do
not substitute for the run's exact certificates or for a demonstrated
end-to-end speedup.
