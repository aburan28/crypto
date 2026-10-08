# Research Note: What the "NIST backdoor in P-256" Is, Is Not, and Would Have To Be

**Status.** Speculation with falsifiers. One fact in this note is newly
*executed* here (the X9.62 seed → `b` derivation, §2); everything else is
reasoning over the public record. No weakness in P-256 is claimed or found.
**Code.** `src/cryptanalysis/p256_backdoor_map.rs` encodes every record below
(`cargo run --release --example p256_backdoor_map` regenerates the Markdown;
`--json` gives the machine-readable form). Companion:
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

## 5. Overall

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

## 6. References

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
