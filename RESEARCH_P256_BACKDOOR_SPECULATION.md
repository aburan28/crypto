# Research Note: What the "NIST backdoor in P-256" Is, Is Not, and Would Have To Be

**Status.** Speculation with falsifiers. One fact in this note is newly
*executed* here (the X9.62 seed → `b` derivation, §2); everything else is
reasoning over the public record. No weakness in P-256 is claimed or found.
**Code.** `src/cryptanalysis/p256_backdoor_map.rs` encodes every record below
(`cargo run --release --example p256_backdoor_map` regenerates the Markdown;
`--json` gives the machine-readable form). Round-2 groundable-locus probes:
`src/cryptanalysis/p256_trapdoor_round2.rs`
(`cargo run --release --example p256_trapdoor_round2`). Companion:
`src/cryptanalysis/p256_speculation.rs` (the public concern vectors and their
probes; that file is not currently declared in `cryptanalysis/mod.rs`, so the
new module carries its own copies of its `PublicResearchStatus` /
`ProbeVerdict` enums).

---

## 0. The one-paragraph answer

There is no demonstrated backdoor in the P-256 *curve*. The documented NSA
backdoor is **Dual_EC_DRBG**, a random-number generator built on two P-256
points, and it says nothing about ECDLP on the curve. The P-256 *concern* is
that the curve's `b` was derived from a seed nobody can explain, so the
"verifiably random" generation is verifiable only from the seed onward. If a
curve-level trapdoor existed it would have to be one of four things (§4), and
each is either excluded by public probes (several of them implemented in this
repository) or unfalsifiable in principle.

## 1. Confirmed: Dual_EC_DRBG

| | |
|---|---|
| Standard | NIST SP 800-90A (2006), instantiated on P-256 with fixed points `P`, `Q` |
| Mechanism | Output = x-coordinates of multiples of `P`, `Q`. If `Q = [d]P` with `d` secret, ~30 output bytes reveal the state; all later output is predictable |
| Shown by | Shumow–Ferguson, CRYPTO 2007 rump session |
| Attribution | 2013 Snowden material (Bullrun); reported RSA BSAFE default-on arrangement |
| Withdrawn | NIST, 2014 |
| Field incident | Juniper ScreenOS, 2015: an unknown party swapped `Q` for its own and could decrypt VPN traffic |
| Implicates curve math | **No** |

## 2. Unresolved: the seed

P-256 was generated at NSA (Jerry Solinas, c. 1999) by the ANSI X9.62 A.3.3.1
procedure with `a = −3` and the Solinas prime `p = 2²⁵⁶ − 2²²⁴ + 2¹⁹² + 2⁹⁶ − 1`
fixed in advance:

```
seed = c49d360886e704936a6678e1139d26b7819f7e90
s = ⌊(256−1)/160⌋ = 1,  v = 256 − 160·s = 96
H  = SHA-1(seed);  c0 = rightmost 96 bits of H;  W0 = c0 with its top bit cleared
W1 = SHA-1((seed + 1) mod 2^160)
r  = W0 ‖ W1 = 0x7efba1662985be9403cb055c75d4f7e0ce8d84a9c5114abcaf3177680104fa0d
check: r · b² ≡ a³ (mod p)      → holds for the published b
```

`verify_p256_seed_derivation()` re-executes this and the test pins `r`. So:

- **What verifies.** `b` really is the X9.62 output of SHA-1 on the published
  seed. The curve was not hand-picked *after* fixing the seed.
- **What does not.** Nothing constrains the seed. A generator who knew a secret
  weak class of curves of density `2^−k` could try `~2^k` seeds until SHA-1
  landed in it (Bernstein–Lange, BADA55; Bernstein et al., *How to manipulate
  curve standards*). The public cannot distinguish that from an honest seed.
- **Counter-argument** (Koblitz–Menezes, *A riddle wrapped in an enigma*): the
  weak class would have to stay unknown to the open community for 25 years,
  and NSA put P-256 and P-384 into Suite B for its own classified traffic.
- **Mundane account.** A former NSA IAD technical director has said the seeds
  were SHA-1 of English phrases that Solinas later lost. A public bounty for
  the preimages (2023) was unclaimed at the time of writing.

## 3. Grind budget — what "could have ground the seed" means in numbers

Assumptions (parameters of `GrindBudget::default_1999`, change them and
re-read the table): 1 000 CPU-years total; `1 µs` per candidate when
membership in the weak class is testable from `(p, a, b)` alone; `600 s` per
candidate when it needs the group order (one 256-bit SEA point count).

| grind model | log₂ candidates | thinnest reachable weak class |
|---|---:|---:|
| coefficient test (SHA-1 + algebraic check) | ≈ 54.8 | 2^−55 |
| point count per candidate | ≈ 25.6 | 2^−26 |

The gap is the whole story: a trapdoor that needs the Frobenius trace to
recognise could only have been ground onto a class of density above about
`2^−26` among prime-order curves, which is 2^100 times denser than any public
weak class; a trapdoor recognisable from the coefficients could be ~2^−55.
Hence the ranking below.

## 4. Curve-level hypotheses, ranked

Rank 1 is the most compatible with the public record and 1999 compute.

### 1. A weak class recognisable from `(p, a, b)` without point counting — `CurveCoefficient`

- **Mechanism.** Grind seeds until `b` lands in a class whose membership is an
  algebraic property of the coefficients alone.
- **Requires.** A secret, cheaply testable invariant of `(p, a, b)` implying an
  ECDLP weakness, unknown to the public for 25 years.
- **Evidence against.** Every public coefficient-level weak class (anomalous,
  small embedding degree, small CM discriminant, singular) is checked and
  absent. Residues `b mod q` over 10⁵ primes pass a Kolmogorov–Smirnov
  uniformity test (`b_seed_profile`); Solinas-reduction micro-bit correlations
  are null (`solinas_correlations`). No detectable algebraic selection on `b`.
- **Verdict.** `NoPublicAnomaly`.

### 2. The prime, not the curve — `FieldPrime`

- **Mechanism.** The Solinas shape was chosen openly for fast reduction; a
  weakness tied to special-form primes needs no grinding, only the knowledge
  that such primes are weak.
- **Requires.** An ECDLP analogue of the finite-field special-prime trapdoor
  (hidden-SNFS primes, Fried–Gaudry–Heninger–Thomé 2016). None is known.
- **Evidence against.** The strongest structure in this prime family is
  cyclotomic smoothness of `p ± 1` (P-224's `p − 1` is 46-bit smooth,
  `RESEARCH_NIST_SOLINAS_STRUCTURE.md`). That accelerates DLP in `F_p^*`,
  which P-256's huge embedding degree never reaches; summation-polynomial index
  calculus measured no gain from it (`RESEARCH_NIST_SOLINAS_EXPERIMENTS.md`,
  experiment 4).
- **Verdict.** `StructurePresentNoKnownAttack`.

### 3. Hidden endomorphism or isogeny structure — `EndomorphismOrIsogeny`

- **Mechanism.** Grind for a curve whose endomorphism ring or isogeny class
  admits a transfer (class-group action, descent to a weak model, small-degree
  isogeny to a special curve).
- **Requires.** A point count per candidate (the invariant is the trace), and
  a weak class of density far above the public ones — capped near `2^−26` by §3.
- **Evidence against.** Huge embedding degree, not anomalous, enormous CM
  discriminant and class number (no CRS/CSIDH-style action), and the
  isogeny-class walk (`p256_isogeny_cover`, `p256_isogeny_walk`) found no weak
  curve at small degree.
- **Verdict.** `NoPublicAnomaly`.

### 4. A non-public ECDLP algorithm for generic prime-field curves — `UnknownAlgorithm`

- **Mechanism.** No trapdoor in the parameters at all.
- **Requires.** A result 25 years of open research has not reproduced; the
  prime-regime index-calculus ladder in this repository
  (`docs/ic/BOUNDARY_TARGETS.md`, regime C) remains above rho at every
  measured size.
- **Evidence against.** Not falsifiable by a finite public check.
- **Verdict.** `NotPubliclyTestable` — the unfalsifiable residual.

## 5. Round 2 (2026-10-12): groundable-locus probes for hypotheses 1–3

Round 1 checked the *public* weak classes. Round 2
(`src/cryptanalysis/p256_trapdoor_round2.rs`, example
`p256_trapdoor_round2`) asks what a *secret but cheaply testable*
class could have been ground into the seed in 1999: every probe is an
invariant of `(p, a, b)` testable in microseconds, exactly the
predicate shape a 2^55-seed grind admits. Recorded outcomes
(`research/p256_trapdoor_round2_20261012/`):

- **Multiplicative orders (H1).** `p − 1` is fully factored:
  `2 · 3 · 5² · 11 · 101 · 281 · 1531 · 65537 · 490719 · 6602113 · Q`
  with `Q` a 160-bit prime. The exact orders are `ord(b) = ord(r) =
  ord(−3)` (253 bits, `Q` included) and `ord(j)` (252 bits, `Q`
  included). No seed-derived value is a small-order or smooth-order
  element mod `p`; a "small-order b" grind locus is dead. The order
  equalities are forced by the shared `Q`-part, not anomalous.
- **Residuosity (H1).** `b`, `r` are quadratic residues and non-cubes;
  `j` is a quadratic and cubic residue. Each is a cheap 1–2-bit
  predicate a grinder *could* have targeted; each is recorded so the
  locus is covered. The `±b` sign choice in the standard was itself a
  free 1-bit grind (`−b` flips the QR bit).
- **Sparse/power forms (H1).** Hamming weights of `b` (z = 0.81),
  `r` (z = −0.31), `j` (z = −0.12) are binomial-normal; none is a
  perfect power; the nearest `c·2^k` is ~2^49 bits away. The one
  near-threshold statistic is `n` (z = 4.88, from the `FF`-run
  structure shared with `p`), but `n` is order-derived — selecting on
  it costs a point count per candidate, capping that grind at the
  §3 budget of ~2^26 candidates.
- **Continued fractions (H1).** `b/p`, `r/p`, `j/p` show max partial
  quotients 1275, 62, 107 over 64 terms — normal-height rationals,
  nowhere near a quadratic irrational.
- **Field prime (H2).** `p − 1` as above (note the Fermat prime 65537
  and 2^96 in `p + 1 = 2^96 · 7 · 273921 · (140-bit composite)`);
  `ord_p(2) = (p−1)/2`; `p ≡ 3 mod 4` with no other notable residue
  profile. No twist subgroup of order ≥ 2^20 has embedding degree
  ≤ 100.
- **Endomorphism/isogeny (H3).** GLS gap `4p − t²` is 258 bits (no GLS
  structure); the trace has Hamming weight 59 (z ≈ −0.9) and no
  sparse form; `x³ − 3x + b` has exactly one rational root, so exactly
  one rational 2-isogeny — the generic case.

**Round-2 verdict.** No new groundable locus found: every
microsecond-testable invariant of the seed-derived values measured
generic. The hypothesis-1 search space is now covered both by round-1
statistics (`b mod q` uniformity, Solinas micro-correlations) and by
round-2 algebraic invariants (orders, residuosity, powers, heights).
What remains for a curve-level trapdoor is either a weak class whose
membership test is *not* one of these cheap shapes, or hypothesis 4.

### Round 3 (2026-10-12): honest-distribution comparison

Round 3 (`src/cryptanalysis/p256_trapdoor_round3.rs`, example
`p256_trapdoor_round3`) tests the *distributional* signature of a
grind: a 1999 grinder who kept only seeds whose `b` satisfied a secret
predicate shifts `b` off the honest X9.62 A.3.3.1 distribution for
every statistic correlated with that predicate. We simulated 2^20
deterministic random seeds, derived `r` and the valid `±√(−27/r)`
candidates (`p ≡ 3 mod 4`, one exponentiation), pooled 1,048,170
roots, and placed the published `b` in each empirical distribution
with Bonferroni correction over the six-probe family
(`research/p256_trapdoor_round2_20261012/round3_records.json`):

- **QR(b)**: honest pool 0.500; published true. Cubic residuosity:
  honest pool 0.332; published false.
- **Hamming weight of b**: honest 128.00 ± 8.00; published 134
  (z ≈ +0.75). Of `r`: honest 127.50 ± 7.99; published 125.
- **CF max quotient of b/p**: honest median bucket ≈ 11 bits, 51% of
  honest samples at or above the published 1275; published 1275 (11
  bits) — dead center.
- All corrected p-values: 1.000.
- **Seed-preimage phrase probe**: SHA-1 over a 1000-phrase curated
  list (Solinas/NIST-era words and suffixes, ≤ 20 bytes) found no
  preimage of `c49d360886e704936a6678e1139d26b7819f7e90`; expected,
  the public bounty remains the reference for this direction.

**Round-3 verdict.** The published `b` is a *median* member of the
honest X9.62 candidate distribution on every probed statistic — there
is no distributional residue of a grind on the round-1/2 invariant
family. Combined with rounds 1–2, a P-256 curve-level trapdoor now
requires a weak class whose membership predicate is neither a public
weak class, nor an algebraic invariant of `(p, a, b)` (orders,
residuosity, powers, heights), nor visible as an outlier in the
honest-candidate distribution of those invariants — or hypothesis 4.

## 6. Overall

The P-256 curve is almost certainly clean; the documented sabotage was
Dual_EC_DRBG plus pressure on implementations and standards. The lasting
damage is that the NIST process cannot *prove* the curve clean, which is why
rigid generation (Curve25519, Brainpool, the CFRG requirements) became the
expectation.

A related observation from this repository's trapdoor work: the NIST binary
curves all use **prime** extension degrees (163, 233, 283, 409, 571), which
excludes the Teske/GHS magic-number trapdoor implemented in `ec_trapdoor.rs`
(it needs a composite degree). That reads as defensive design.

**What would change this assessment.** A seed preimage (closes the question
in the mundane direction); a public ECDLP method that beats rho on generic
prime-field curves (opens hypothesis 4); a non-uniform statistic on `b` or a
coefficient-level invariant with an attack behind it (opens hypothesis 1).

## 7. References

- ANSI X9.62-1998, Annex A.3.3 (pseudo-random curve generation and
  verification); FIPS 186-4, Appendix D.1.2.3; SEC 2 v2.0, §2.4.2.
- D. Shumow, N. Ferguson, *On the possibility of a back door in the NIST
  SP800-90 Dual Ec Prng*, CRYPTO 2007 rump session.
- S. Checkoway et al., *On the practical exploitability of Dual EC in TLS
  implementations*, USENIX Security 2014; *A systematic analysis of the Juniper
  Dual EC incident*, CCS 2016.
- D. J. Bernstein, T. Lange, *SafeCurves* (rigidity criterion), 2014;
  Bernstein, Chou, Chuengsatiansup, Hülsing, Lambooij, Lange, Niederhagen,
  van Vredendaal, *How to manipulate curve standards: a white paper for the
  black hat*, SSR 2015 (the BADA55 curves).
- N. Koblitz, A. Menezes, *A riddle wrapped in an enigma*, IEEE Security &
  Privacy 2016 (eprint 2015/1018).
- J. Fried, P. Gaudry, N. Heninger, E. Thomé, *A kilobit hidden SNFS discrete
  logarithm computation*, EUROCRYPT 2017.
- E. Teske, *An elliptic curve trapdoor system*, J. Cryptology 19(1), 2006.
- Public seed-preimage bounty announced by F. Valsorda, 2023.
