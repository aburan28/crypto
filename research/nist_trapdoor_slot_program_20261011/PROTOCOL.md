# NIST prime-curve trapdoor slots: a preregistered experimental portfolio (registered 2026-10-11)

**Status.** Protocol with one executed round (the census, §3) and seven
pending experiments (§4–§11).  No weakness in any NIST curve is claimed
or found.  Everything in this file is one of: *executed* (with a command
and a hash), *derived* (a formula or a published theorem, cited),
*proposed* (a hypothesis with a falsifier), or *not publicly testable*.

**Question.** If a covert ECDLP trapdoor had been inserted into a NIST
prime curve (P-192 … P-521) so that one party could solve ECDLP in
sub-exponential time while the public could neither solve it nor
recover the trapdoor cheaper than solving it, *where in the published
parameters would the trapdoor have to live, and what experiment, run
on public data, would bound or exclude each location?*

**What this adds to `RESEARCH_P256_BACKDOOR_SPECULATION.md`.** That note
ranks four hypotheses for P-256.  This protocol (i) enumerates every
parameter slot of the X9.62 / FIPS 186-4 generation procedure rather
than only the curve equation, (ii) runs the public-parameter checks on
all five prime curves rather than P-256 alone, and (iii) attaches a
registered experiment, with budget and falsifier, to each slot.  It
reuses the repository's existing probes wherever one exists and names
the gap where one does not.

Companion file: `README.md` (reading order and evidence index);
`census.md` / `census.json` (executed round 0); `slots.svg` (diagram).

---

## 0. The organising principle: invariants cannot hide, models can

Abstractly `E(F_p)` of prime order `n` is `Z/n`.  ECDLP is hard only
because of the *representation*.  A trapdoor is therefore, by
definition, a withheld efficiently-computable homomorphism from the
published representation into one with a sub-exponential DLP.  Two
consequences fix the shape of every experiment below.

1. **Isogeny invariants are public and checkable.**  `#E`, the trace
   `t`, the embedding degree `k = ord_n(p)`, the Frobenius discriminant
   `t² − 4p` and its square-free part, and hence the CM field, are
   constant on the whole `F_p`-isogeny class.  Any weakness that is a
   function of these (anomalous, MOV/Frey–Rück, small class group) is
   excluded by a public check and cannot be hidden by an isogeny walk.
2. **Model-dependent properties can be hidden.**  A weakness that
   depends on the specific Weierstrass model or on the field of
   definition (the GHS magic number, a specific cover, a specific
   factor-base structure) is *not* isogeny-invariant.  The designer
   starts at a weak model `E₀`, walks a secret isogeny chain to `E`, and
   publishes `E`.  Recovering the chain from `E` alone, with the target
   class unknown, is the asymmetry the trapdoor needs.  This is
   Teske's construction (J. Cryptology 2006) over `F_{2^n}`, `3 | n`,
   and it is the only published instance.

So the only slot that can hold a sub-exponential ECDLP trapdoor is the
representation itself (S1).  Every other slot is either a public
invariant (checkable), a protocol-level trapdoor that does not touch
ECDLP (S4), an implementation-level weakness (S5), or non-weakening
(S7).  The experiments are arranged so that the checkable slots are
*executed* and the representation slot is *bounded* at toy scale by
registered, falsifiable predictions.

## 1. Slot taxonomy

| slot | parameter | kind of trapdoor possible | isogeny-invariant? | status in this protocol |
|---|---|---|---|---|
| S1 | the group representation `(p, a, b)` up to isogeny | hidden homomorphism to a sub-exponential group | **no** (that is the point) | S1a–S1d registered, pending |
| S2 | the field prime `p` | special-prime analogue of hidden-SNFS | n/a (field, not curve) | census executed; S2 pending |
| S3 | the seed → `b` map, including the root sign | selection freedom for S1 | n/a | census executed (new fact); S3 pending |
| S4 | the base point `G` | protocol trapdoor (Dual_EC class), not ECDLP | n/a | census executed; S4 pending |
| S5 | the quadratic twist | invalid-curve key leak, not ECDLP | twist order is invariant; its *security* is a cofactor fact | census executed; S5 pending |
| S6 | `n`, cofactor, `n ± 1` | none known | yes | census executed |
| S7 | seed as a fingerprint | non-weakening authorship proof | n/a | census executed (public dictionary only) |
| S8 | a better generic algorithm | capability, not trapdoor | n/a | not publicly testable; see boundary ledger |

## 2. Frozen inputs for every round

- Curves: `CurveParams::p192()`, `p224()`, `p256()`, `p384()`, `p521()`
  from `src/ecc/curve_zoo.rs` (`a = −3`, FIPS 186-4 Appendix D.1.2).
- Seeds: FIPS 186-4 / SEC 2 v2.0 §2.4 values, pinned in
  `examples/nist_slot_census.rs::seed_for`.
- Derivation: ANSI X9.62-1998 A.3.3.1 as implemented in
  `src/cryptanalysis/p256_backdoor_map.rs::x962_candidate_r` /
  `verify_seed_derivation` (already pinned by a test for P-256).
- Repository revision for round 0: `b5bd97f2d302370be58d2cf31f75bd1490a96c99`
  plus this PR's files.  Toolchain: rustc 1.97.0 (2d8144b78 2026-07-07),
  2 × Intel Xeon @ 2.10 GHz, Linux 6.18.
- Trial-division bound `10⁶`; Pollard rho budget `2²²` iterations per
  split; embedding-degree scan `k ≤ 1000`.

## 3. Round 0, executed: the public-parameter census

**Command.**

```text
cargo build --release --example nist_slot_census
./target/release/examples/nist_slot_census          > census.md
./target/release/examples/nist_slot_census --json   > census.json
```

Wall time 35 s (the Pollard rho budget on the five twist orders).

**Hashes.**

```text
sha256(examples/nist_slot_census.rs) = c59d84f8071c2bebd512abdeb206ccac6aa6c3165e216e7f8520ef021116746a
sha256(census.json)                   = 3467a6d68eae2b57f76a3db020ca64794bbccd713462ce9f5e060ec710d701db
sha256(census.md)                     = 9f7419e3668ea20529f6950c5053ebe0580779816689cf0b111a8266e09c2a56
```

**Observed (from `census.md`; every row recomputed, nothing pinned).**

Group invariants (S2, S5, S6).  All five curves: prime order, cofactor
1, trace positive, embedding degree `> 1000`.  `|D| = 4p − t²` has no
square factor from primes below `10⁶` except P-224 (`3³` contributes a
4-bit square part); no residue is a perfect square; the P-192 residue
is a 184-bit probable prime.  So no public small-class-group or
small-conductor structure is reachable by trial division on any of the
five.  `n − 1` and `n + 1` carry only small factors below the bound with
large residues (P-224's `n − 1` residue is a 196-bit probable prime).

Seed → `b` and the root sign (S3).  `r(seed) · b² ≡ a³ (mod p)` holds
for all five published seeds.  The procedure fixes `b²`, not `b`; the
generator chose a sign.  For P-192, P-256, P-384 and P-521, `p ≡ 3
(mod 4)`, so `−1` is a non-residue and `b ↦ −b` is the *quadratic
twist*: the sign choice selects between a curve and its twist after the
hash output is known.  For P-224, `p ≡ 1 (mod 4)`, `i ∈ F_p`, and
`b ↦ −b` is an `F_p`-isomorphism `(x, y) ↦ (−x, iy)`: no freedom.  The
published root is the numerically smaller one for P-192, P-256 and
P-521 and the *larger* one for P-224 and P-384.  No published
convention explains the choice; this is one extra, previously
unremarked bit of generator freedom on four of the five curves, beyond
the seed itself.  It is not evidence of grinding: it is one bit.

Generator (S4).  `G` is on every curve.  `G.x` is never the
lexicographically first valid `x` (that is `x = 0` on P-192, P-256,
P-384, P-521 and `x = 3` on P-224), and `G.y` is the larger of the two
roots on P-224 and P-521.  So `G` is not the deterministic
"first point" a rigid procedure would produce; its provenance is
unexplained on all five curves, as the P-256 literature already notes.

Twist (S5).  Twist orders are composite for P-192, P-224, P-256,
P-521 and prime for P-384.  With the rho budget: P-256's twist is
`3·5·13·179·q₂₄₁` (`q₂₄₁` a 241-bit probable prime), P-521's is
`5·7·69_697_531·635_884_237·q₄₆₁`, P-224's is `3²·11·47·3_015_283·40_375_823·c₁₆₆`
with a 166-bit composite the budget did not split, and P-192's is
`23·c₁₈₈`.  The P-224 twist is the known weak one (SafeCurves: rho on
the twist near `2^58`); the unsplit composite is exactly where that
figure lives and is S5's pending item.

Seed relations (S7).  84 public nothing-up-my-sleeve candidates per
curve (SHA-1 of the curve's names, of `p` and `n` in decimal and hex,
with and without newline; the other four seeds; SHA-1 of each seed and
of each seed `+ 1 … 4`; `name ‖ counter`): zero matches on all five.
This excludes only those families; a `SHA-1(secret ‖ counter)` seed is
undetectable by construction.

**What round 0 does and does not establish.**  It establishes that no
public, trial-division-reachable weak class is present on any of the
five curves and that the generation procedure leaves two unexplained
degrees of freedom per curve besides the seed: the root sign (four
curves) and `G` (all five).  It does not bound S1, which is the only
slot that can carry a sub-exponential trapdoor; S1 is the subject of
the pending rounds.

## 4. S1a — prime-field cover route (pending)

**Hypothesis H1a.** There is no cover `C → E` with `C` defined over
`F_p` (or over `F_{p^k}` with `E` base-changed) whose Jacobian index
calculus beats rho on `E(F_p)`; the complexity exponent
`Õ(p^{2 − 2/g})` of double-large-prime index calculus on `Jac(C)` is
above rho's `1/2` for every genus, and base change to `F_{p^k}` only
raises it.  *Derived* already as Prediction A-1 in
`research/cover_route_prime_binary_20261008/PROTOCOL.md`.

**What is not derived.** Whether a *specific* small-degree cover can be
hidden by an isogeny walk and still transfer the subgroup.  Over `F_p`
the answer is "irrelevant" if H1a holds (the destination solver is
never faster), so the experiment is a toy-scale confirmation of the
exponent, not a search.

**Experiment E1a.** For `p ∈ {2^20 + c, 2^24 + c, 2^28 + c}` (first
primes above those powers), 32 random prime-order curves each.  For
each, construct the degree-2 covers `C: w² = f(x)` with `f` the
cubic pulled back through `x ↦ x²` and `x ↦ x² + u` (genus 2) and the
degree-3 cover through the Kummer line (genus 3 or 4), using
`src/cryptanalysis/jv_cover.rs` and
`research/hyperelliptic_cover_infrastructure_20261006`.  Run the
repository's hyperelliptic index calculus on `Jac(C)` with the standard
large-prime variant and rho on `E`, both in ecbench's counted unit.
**Prediction:** measured `log_p(cost(Jac))` slope across the three `p`
is within `0.1` of `2 − 2/g` and above `0.5 + 0.05` for every cover.
**Falsified if** any cover family has slope below `0.55`.  Stop
condition: 3 `p` values × 32 curves × 3 cover families, or 48 CPU-hours.

## 5. S1b — class-group map (pending)

**Hypothesis H1b.** No efficiently computable homomorphism
`E(F_p) → Cl(O)` exists for `O = End(E)` unless `gcd(n, h(O)) > 1`,
because a homomorphism of prime-order `Z/n` into a group of order `h`
is trivial when `gcd(n, h) = 1`.  *Derived* (Lagrange).  For the NIST
curves `h(O)` is not computable (`|D|` is 194–523 bits with no public
small-conductor structure, census §3), so H1b cannot be tested on them;
it can only be tested on the *density* question: how often does a
random prime-order curve have `gcd(n, h(O)) > 1`, and could a 1999
grind reach that class?

**Experiment E1b.** For `p` the first prime above `2^28`, sample 2¹⁶
random prime-order curves (point count by the repository's SEA or BSGS
at that size), compute `D = t² − 4p`, its square-free part and `h(D)`
by the repository's class-number routine (`hilbert_class_poly.rs`
machinery or a baby-step enumeration of reduced forms at 30 bits), and
record `gcd(n, h)`.  **Prediction:** the fraction with `gcd(n, h) > 1`
is below `2^{−12}` and every such case has `h` with a small factor
dividing `n` trivially (since `n` is prime, `gcd > 1` means `n | h`,
which for `h ≈ √|D| ≈ p^{1/2} < n` is impossible).  So the registered
prediction is **zero** such curves; **falsified if** any curve has
`n | h(O)`.  This is a sanity bound, recorded so that S1b has a
measured row rather than an argument.  Cost: one point count and one
class-number computation per curve at 30 bits, ≈ 2 CPU-hours.

## 6. S1c — hidden factor-base structure for summation-polynomial index calculus (pending)

**Hypothesis H1c.** The point-decomposition cost in Semaev index
calculus over `F_p` depends on `p` and the decomposition arity, not on
the curve: isogenous curves have equal decomposition success
probability and equal Gröbner / SAT solve cost to within noise, so no
model in an isogeny class is a "better factor base carrier" than
another.  *Proposed.*

**Experiment E1c.** For `p` the first prime above `2^24`, take 8
random prime-order curves and for each walk 16 steps of the
repository's non-backtracking `ell`-isogeny walk (`isogeny_walk/`,
`p256_isogeny_walk.rs` generalised to small `p`).  On each of the 136
vertices, run the frozen Semaev arity-3 decomposition harness
(`semaev_decomp.rs`, the boundary-ledger protocol for regime C) with
the same factor-base construction (`x`-coordinate below `p^{1/3}`),
1024 random targets each.  Record decomposition success rate and mean
solve time in the counted unit.  **Prediction:** within each isogeny
class the success rate is constant to within binomial `3σ` and the
solve time within `10 %`.  **Falsified if** any vertex differs from
its class mean by more than that on both metrics.  A falsification
would be the first evidence that a model-dependent factor-base
advantage exists, and would be the prime-field analogue of the GHS
magic number.  Cost: 136 × 1024 decompositions at 24 bits, ≈ 6
CPU-hours.

## 7. S1d — isogeny-hiding asymmetry (pending)

**Hypothesis H1d.** (i) Every public weak-class predicate (anomalous,
`k ≤ 20`, square-free part of `D` below `2^40`) is constant along the
walk; (ii) recovering a hidden chain of `L` steps of degree `ell` to a
*known* target costs `O(√L · ell)` by meet-in-the-middle and to an
*unknown* target costs no less than enumerating the class, so the
asymmetry a Teske-style trapdoor needs exists over `F_p` *for the
hiding step*, even though no weak class exists to hide.  (i) is
*derived* (Tate); (ii) is *proposed* at the measured level.

**Experiment E1d.** For `p` the first prime above `2^28`, 8 random
prime-order curves, walks of `L ∈ {8, 16, 32}` steps at the smallest
split `ell`.  Assert (i) on every vertex with the repository's
invariant checks.  For (ii), implement the bidirectional
meet-in-the-middle recovery with the target known, and the
class-enumeration recovery with only the *predicate* known (here the
predicate is "has `j`-invariant in a random 1 % subset", standing in
for a secret class), both counting modular-polynomial root
extractions.  **Prediction:** known-target recovery cost grows as
`ell^{L/2}` within a factor 4; predicate-only recovery grows as the
class size divided by the subset density, i.e. is independent of `L`
and ≈ `100 × √|D|/ell`-ish at this `p`; the ratio between the two exceeds
`2^10` at `L = 32`.  **Falsified if** predicate-only recovery is within
`2^5` of known-target recovery at `L = 32`.  Cost: ≈ 4 CPU-hours.

## 8. S2 — the prime (pending, one new control)

**Hypothesis H2.** The Solinas shape of `p` confers no advantage to
any ECDLP method the repository can run.  Covered for
summation-polynomial index calculus by
`research/notes/cm-isogeny/RESEARCH_NIST_SOLINAS_EXPERIMENTS.md` experiment 4 (measured no gain).
`p ± 1` smoothness only reaches `F_p^*`, which `k > 1000` (census)
never touches.

**Experiment E2.** The one control that file did not run: matched rho.
Run `p256_negation_rho` (`p256_attacks.rs`) on P-256 and on a random
256-bit prime-order curve over a random (non-Solinas) 256-bit prime
with the same ecbench spec, 2²⁰ iterations each, and record the
per-iteration counted cost.  **Prediction:** the Solinas reduction is
faster per field multiplication (that is its purpose) and the
*iteration count to collision* is statistically identical; the
protocol is a stage diagnostic and makes no end-to-end claim.
**Falsified if** the Solinas curve needs fewer iterations at `3σ`.

## 9. S3 — seed grinding density at toy scale (pending)

**Hypothesis H3.** The X9.62 procedure at toy scale reproduces the
grind-budget table of `RESEARCH_P256_BACKDOOR_SPECULATION.md` §3: the
density of prime-order outputs is `≈ 1/ln p`, and the density of each
public weak class among prime-order outputs is as the literature gives
(anomalous `≈ 1/√p`-ish per Hasse interval, `k ≤ 20` `≈ 20/p`,
square-free part of `D` below `B` ≈ `B/√p`).

**Experiment E3.** Implement the full X9.62 A.3.3.1 loop natively for
`p` of 40 bits (first prime above `2^40`), `a = −3`, 2²⁴ sequential
seeds.  For each candidate compute `b` (both roots), the curve order
(the repository's point counter at 40 bits), and the public weak-class
predicates.  Record densities with Wilson intervals.  **Prediction:**
each density within a factor 2 of the formula.  **Falsified if** any
public class is denser than `2^{−16}` among prime-order outputs.  A
second output of E3 is the *root-sign* statistic: whether the smaller
root yields a prime-order curve more or less often than the larger
(prediction: no difference at `3σ`).  Cost: 2²⁴ point counts at 40
bits, ≈ 20 CPU-hours; run on isolab per `ecbench-independent-runner`.

## 10. S4 — base point (pending, two items)

**Experiment E4a (rigid-derivation test).** Pin the ANSI X9.62-1998
Annex A.3.2 base-point derivation and the SEC 1 / FIPS 186-4 text on
`G` (none of which specify a seeded derivation for the NIST `G`s; this
must be confirmed from the standards' text and recorded in the round
report).  Then test each published `G` against every candidate
deterministic rule the round can state exactly: lexicographically
first `x` (done, census: no), `x = SHA-1(seed ‖ i)` reduced for
`i ≤ 2^16`, `x = SHA-1(seed)` with the X9.62 `r` construction, and
`G = [h]·(first point)` with `h = 1`.  **Prediction:** no match.  A
match would close S4 in the mundane direction for that curve.

**Experiment E4b (protocol audit, no computation).** Enumerate deployed
protocols on NIST prime curves that expose scalar multiples of *two*
fixed points or of `G` plus a long-term point in a way that makes a
secret relation `Q = [d]G` exploitable (the Dual_EC pattern):
Dual_EC_DRBG (withdrawn), ECDSA with deterministic nonce derivation
(no), ECIES/ECDH (no second fixed point), the ECQV implicit-certificate
scheme, and the EC-based KEMs in the ISO 18033 family.  Deliverable is
a table with the trapdoor condition per protocol.  This is a literature
round and makes no ECDLP claim.

## 11. S5 — twist completion (pending)

**Experiment E5.** Finish the factorisation of the P-192 and P-224
twist orders (the 188-bit and 166-bit unsplit composites of census §3)
with ECM in the repository's native factoring tools (or a bounded
`cargo` ECM), so the "rho on the twist" cost is a measured number
rather than a literature citation.  **Prediction:** the P-224 twist's
largest prime factor is near `2^117` (SafeCurves `2^58.4` rho cost ×2)
and P-192's is above `2^150`.  Then re-run the invalid-curve harness
(legacy Python in the repository; per `AGENTS.md` it must be replaced
by a native path before being extended, so this round includes that
replacement) against a test ECDH endpoint that skips point validation,
and record the per-query leak in bits, which the census already bounds
from below by the small factors (`15.1` bits for P-256, `12.2` for
P-224).  This is an implementation-level result and is reported as
such.

## 12. S6, S7, S8 — closed or not testable

- **S6.** Executed in census §3.  No pending item: `n` prime, cofactor
  1, no known consequence of `n ± 1` smoothness.
- **S7.** The public dictionary is executed (zero matches).  The
  remaining family, `SHA-1(secret ‖ counter)`, is undetectable; the
  only closing evidence is a published preimage (the 2023 bounty).  No
  experiment can be registered.
- **S8.** A classified generic improvement is not a parameter trapdoor
  and no public parameter check bears on it.  The repository's
  evidence on this is the boundary ledger
  (`docs/ic/BOUNDARY_TARGETS.md`, regime C), which records the measured
  prime-field index-calculus frontier against rho.  Nothing is
  registered here.

## 13. Claim hygiene and what success would mean

Per `AGENTS.md`: every round reports in the counted unit against the
ledger, fail-closed, synthetic fixtures only; no key recovery, no
sub-rho asymptotic claim, no deployed-curve-security claim.  The
portfolio's possible outcomes are:

- all predictions hold: the slots S2–S7 are *measured* clean on all
  five curves and S1 is *bounded* at toy scale; the residual is S8
  and the undetectable S7 family;
- a falsification in E1c or E1d: a model-dependent property that
  survives isogeny hiding over `F_p` exists at toy scale and becomes a
  research lead for the index-calculus ladder, reported as a stage
  diagnostic with the exact curve identities;
- a falsification in E3: the X9.62 procedure at toy scale is biased in
  a way the literature formulas do not predict, which would re-open the
  grind budget of `RESEARCH_P256_BACKDOOR_SPECULATION.md` §3.

None of these outcomes would by itself establish that a NIST curve is
or is not trapdoored.  The portfolio bounds *where* a trapdoor could
be and *what it would cost to have placed one*, which is the most a
public check can do.

## 14. References

- ANSI X9.62-1998, Annex A.3.3; FIPS 186-4, Appendix D.1.2; SEC 2 v2.0 §2.4.
- E. Teske, *An elliptic curve trapdoor system*, J. Cryptology 19(1), 2006.
- D. J. Bernstein et al., *How to manipulate curve standards*, SSR 2015.
- N. Koblitz, A. Menezes, *A riddle wrapped in an enigma*, eprint 2015/1018.
- P. Gaudry, F. Hess, N. P. Smart, *Constructive and destructive facets of Weil descent*, J. Cryptology 2002.
- C. Diem, *On the discrete logarithm problem in elliptic curves*, Compositio 2011.
- S. D. Galbraith, F. Hess, N. P. Smart, *Extending the GHS Weil descent attack*, EUROCRYPT 2002.
- D. Shumow, N. Ferguson, CRYPTO 2007 rump session (Dual_EC).
- D. J. Bernstein, T. Lange, *SafeCurves*, twist-security tables.
- Repository: `RESEARCH_P256_BACKDOOR_SPECULATION.md`,
  `research/cover_route_prime_binary_20261008/PROTOCOL.md`,
  `research/p256_isogeny_walk_20261004/PROTOCOL.md`,
  `src/cryptanalysis/ec_trapdoor.rs`.
