# SafeCurves / SECG exotic-trapdoor slots (registered 2026-10-17)

**Status.** Protocol with two executed rounds (a sextic-twist-orbit census of
the SECG `j = 0` curves, and a toy-scale grind demonstration) and registered
follow-ups. No weakness in any deployed curve is claimed beyond what the
census rows state; no key recovery is performed. Continues
`research/nist_trapdoor_slot_program_20261011/PROTOCOL.md` (slots S1-S8 for
the NIST prime curves) by asking the mirror question: **which trapdoor
families could *pass* the SafeCurves criteria (or the weaker SECG criteria),
and where would they have to live?**

**Companion files.** `census.md` / `census.json` (round A),
`grind_demo.md` / `grind_demo.json` (round B), `twist_orbit_census.py`,
`sextic_grind_demo.py` (both SageMath; every value recomputed, nothing
pinned by hand), `slots.svg`, `PROTOCOL.pdf`.

---

## 0. Question

SafeCurves certifies, per curve: ECDLP rho cost, no pairings transfer
(embedding degree), CM discriminant `|D| > 2^100`, twist security (rho on the
quadratic twist `> 2^100`), ladder support, completeness, and rigidity of
generation. SECG (SEC 2) requires far less: essentially prime-ish order and
documented parameters, with no twist-security or rigidity requirement at all.
If a designer wanted a curve that *passes* these checks while retaining a
secret sub-exponential or intermediate-cost advantage, which criterion is the
one that cannot see it?

## 1. Criteria-versus-trapdoor map (derived; slots T1-T6)

| slot | trapdoor family | killed by | can it pass SafeCurves? |
|---|---|---|---|
| T1 | group representation `(p, a, b)` up to isogeny: hidden homomorphism / model-dependent solver advantage (Teske class) | nothing in SafeCurves; discriminant check kills only class-group maps, transfer kills pairings | **yes** — this is the only slot with no public check |
| T2 | special prime `p` (hidden-SNFS analogue in `F_p^*`, reached only via a transfer) | transfer criterion (embedding degree) unless the transfer is not a pairing | yes, if the transfer is not a pairing-based one |
| T3 | seed grinding onto a secret coefficient-recognisable weak class | rigidity criterion *in spirit*, but rigidity is satisfied by any *documented* deterministic procedure with unexplained constants (BADA55 pattern) | yes — "explained" generation does not mean "ungrindable" generation |
| T4 | quadratic-twist weakness (invalid-curve leak) | twist criterion | no (SafeCurves); **yes for SECG**, which does not require twist security |
| T5 | **sextic/quartic twist weakness on automorphism-rich models (`j = 0`, `j = 1728`)** | nothing: the twist criterion certifies only the quadratic twist, one of six classes at `j = 0` | **yes** — this round's finding |
| T6 | base-point protocol trapdoor (Dual_EC class) | not an ECDLP criterion anywhere | yes, but not curve-mathematical |

T5 is the round's new slot. It is not "the standard forgot to say
something about twists": the twist criterion *as specified* certifies one
member of a six-member orbit, and round B shows the designer's cost to
exploit that is negligible.

## 2. Round A, executed: sextic-twist-orbit census of the SECG j = 0 curves

**Command.** `sage twist_orbit_census.py` (SEA point counting, BPSW
probable primes; `census.md` / `census.json` are its output).

For `p = 1 mod 6` the curves `y^2 = x^3 + b*d` fall into six `F_p`-isomorphism
classes with **six different group orders** `{p+1 +/- t_i}`, the six traces
satisfying the classical `j = 0` CM relation `t_1 = t_2 + t_3` (verified on
secp256k1 from the census data; the whole orbit is therefore derivable from
the public `(p, n)` alone, which is what makes a public check cheap).

Observed (rho cost = half the bit length of the largest prime factor, per
class):

| curve | published | quadratic twist (the criterion's target) | weakest *other* sextic class |
|---|---:|---:|---:|
| secp192k1 | 2^96 | 2^68.5 | **2^28.0** (witness d = 3) |
| secp224k1 | 2^112.5 | 2^111.5 | **2^55.5** (witness d = 8, the criterion's own class, and d = 14 at 2^56.5 is another) |
| secp256k1 | 2^128 | 2^110 | **2^65.5** (witness d = 3) |

Reading notes, kept honest:

- The quadratic twist of secp256k1 is *strong* (2^110); the well-known weak
  secp256k1 twist figure belongs to a **cubic** twist (class d = 3), which
  the SafeCurves-style check never looks at.
- secp224k1's quadratic twist is the weak one *and* two more classes are at
  2^55.5-2^56.5, so for that curve the criterion happens to catch a weak class
  -- by certifying only one of three.
- Reachability: x-only ladders reach only the quadratic twist (a point with
  `x^3+b` a non-residue lies on exactly that class), so Montgomery-style
  implementations are covered. Implementations that accept a supplied `y`
  without validating `y^2 = x^3+b` can be fed points on *any* class (and
  indeed on any curve sharing `a`), which is the classic invalid-curve
  setting; the marginal content of T5 is that for `j = 0` curves the weak
  class is a *published-parameter-determined* one, at a measured cost far
  below the certified twist. This is an implementation-level observation
  about what "twist security" certifies, reported as such.

## 3. Round B, executed: the grind is nearly free (toy scale)

**Command.** `sage sextic_grind_demo.py`.

Model: the designer picks `p = 1 mod 6` (free under SECG-style generation)
and the sextic class of `b`. Scaled SafeCurves-style acceptance: published
class and its quadratic twist both have largest prime factor `>= 0.79 *
bitlen(p)`; the withheld trapdoor is some other class with largest prime
factor `<= 0.40 * bitlen(p)`.

| size | primes scanned | published class largest factor | twist class | weakest sibling (rho) |
|---:|---:|---:|---:|---:|
| 48 bits | 2 | 48 bits | 40 bits | 13 bits (2^6.5) |
| 96 bits | 3 | 95 bits | 79 bits | 28 bits (2^14) |

**Prediction (registered before scaling further):** at 256 bits the expected
number of `p = 1 mod 6` primes to scan before a curve passing the scaled
pair-check with a sibling below `2^102` (rho `2^51`) is below `2^12`; each
candidate costs six SEA counts. **Falsified if** no such pair exists in
`2^13` scanned primes. Round C (pending) runs this natively in Rust per the
repository's no-Python rule, alongside the census cross-check via the CM trace
relation instead of SEA.

## 4. What survives for SafeCurves-conforming models (the T5 boundary)

Curves that pass the *ladder* criterion are Montgomery/Edwards models with
`j != 0, 1728`; those have exactly one nontrivial twist (the quadratic one),
which the twist criterion fully covers. So T5 lives only on
short-Weierstrass `j = 0` / `j = 1728` models — exactly the SECG `secp*k1`
family — and the exotic-trapdoor search for *SafeCurves-grade* criteria
reduces to T1 (representation) and T3 (rigidity that is not actually
manipulation-proof), with T5 the SECG-grade instance. This is consistent
with the NIST program's S1 conclusion and sharpens it: the only slot with no
public check at either criteria level is the representation itself.

## 5. Registered follow-ups

- **R1 (round C, native).** Rewrite rounds A and B in Rust (SEA or the CM
  trace derivation `4p = L^2 + 27M^2`, repository factoring tools), run the
  256-bit grind prediction of §3, and register the full twist-orbit check
  (`largest prime factor of all six classes >= 2^{100/2}`-scaled) as a
  public rigidity-adjacent criterion candidate.
- **R2.** Quartic analogue at `j = 1728` (`y^2 = x^3 + b`, `p = 1 mod 4`):
  same census and grind on the four-class orbit; no deployed `j = 1728`
  curve is known to us, so this is criteria-completion only.
- **R3.** Reachability experiment: an instrumented x-only ladder and a
  y-accepting implementation on the round-B 96-bit curve, confirming which
  classes each can be pushed onto, so §2's reachability note is measured
  rather than argued.

## 6. Claim hygiene

Per `AGENTS.md`: synthetic fixtures only; no deployed-curve-security claim
beyond the census rows themselves; no key recovery; the grind demonstration
is a construction at toy scale, not evidence that any SECG curve was so
generated. The twist-orbit facts are measured; the T1/T3 conclusions are
derived reasoning over published criteria, and the 256-bit grind prediction is
a registered extrapolation until round C lands.

## 7. References

- D. J. Bernstein, T. Lange, *SafeCurves: choosing safe curves for
  elliptic-curve cryptography*, safecurves.cr.yp.to (criteria tables).
- SEC 2 v2.0, *Recommended Elliptic Curve Domain Parameters*.
- E. Teske, *An elliptic curve trapdoor system*, J. Cryptology 19(1), 2006.
- D. J. Bernstein et al., *How to manipulate curve standards*, SSR 2015
  (BADA55: passing "explained-generation" bars with a ground seed).
- Repository: `research/nist_trapdoor_slot_program_20261011/PROTOCOL.md`
  (slots S1-S8), `RESEARCH_P256_BACKDOOR_SPECULATION.md`.
- Classical `j = 0` twist theory: the six sextic twists and their traces
  `{p+1 -/+ t_i}`; Silverman, *The Arithmetic of Elliptic Curves*, Ex. 5.23
  region (twists and `F_p`-isomorphism classes).
