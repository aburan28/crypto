# Plan: what is left for summation-polynomial decomposition on ECC2K-130

Written 2026-10-16. **This is a plan, not a protocol.** Each workstream below becomes
its own experiment directory, with a `PROTOCOL.md` frozen and committed before
anything it cites is run. Nothing here is a measurement, and nothing here claims work
below `2^61`.

**Success criterion** (owner, 2026-10-15; AGENTS.md). Any verified, measured improvement
over the current index-calculus baseline is a success. Each result reports its gap to
the matched rho as the distance remaining, not as pass or fail. The kill criteria below
close a branch only when it cannot improve the baseline.

**Implementation language** (AGENTS.md, "Implementation language: no Python"). Every
workstream is built in Rust, inside the existing library and its Cargo tests.
- The closed experiments listed below were written in Python. They stay as legacy
  evidence, with their Python provenance labelled.
- Any piece a workstream reuses from them gets a native port first: `F_{2^n}` arithmetic,
  the `S3` linearization and the pair oracle. Use the library's existing binary-field and
  Koblitz code where it already covers them.

## Where things stand

The curve is `y² + xy = x³ + 1` over `F_{2^131}`, with `#E = 4r` and a matched rho cost
of `2^60.8090`.

The best decomposition found uses the linear `(A, P)` pair oracle over a geometric factor
base `V = θ⟨1, g, …, g^{l−1}⟩`
([`../geometric_v_linear_20261006`](../geometric_v_linear_20261006/README.md)). It stays
far above rho:

| | gap above rho at `n = 131` |
|---|---|
| bits only | 42.47 bits |
| charged with linear algebra (`ω = 2`) | 55.22 bits |

These branches are closed, each by a frozen experiment:

| branch | result | where |
|---|---|---|
| higher-degree XL on the bilinear core | best charged gap 71.82 bits | #1420 `xl_bilinear_core_20261006` |
| XL on the trilinear system, `t` unknown | 87.16 bits, worse than the bilinear core | #1429 `xl_trilinear_t_unknown_20261006` |
| three-summand linearization of `S4 = G²` | `31l − 24` unknowns, so `l ≤ 5` | #1442 `geometric_v_linear_20261006` |
| graph-cycle IC (no linear algebra) | rho-dominated: needs about `2^{n−l}` tests | #1510 `graph_cycle_ic_20261007` |
| torsion and model invariants | at most 2 bits per summand | #1509 `cheap_closures_20261007` |
| excess additive energy of geometric `V` | none: ratio 0.94–1.01 to a size-matched random set | #1509 `cheap_closures_20261007` |

## Three facts found while scoping

These came from throwaway scoping checks run in the planning session. Those checks were
Python, so they are not committed and not cited. Facts 2 and 3 are re-established in
Rust under frozen protocols: the bound check in W1, the budget table in W3.

Fact 1 is exact arithmetic (the order of 2 modulo `n`), but the W2 protocol still
recomputes it.

1. **No small Frobenius-closed factor base exists at `n = 131`.**
   - The order of 2 modulo 131 is 130, so `x^131 − 1 = (x − 1)·Φ₁₃₁(x)` with `Φ₁₃₁`
     irreducible over `F₂`.
   - Squaring is `F₂`-linear on `F_{2^131}`. Its invariant subspaces are therefore only
     `{0}`, `F₂`, the trace-zero hyperplane, and the whole field.
   - The toy fields `n = 13, 19` behave the same way (orders 12 and 18). `n = 23` does
     not: the order of 2 is 11, so `F_{2^23}` has 11-dimensional Frobenius-invariant
     subspaces. That makes `n = 23` a ready-made positive control and an unrepresentative
     test field.
2. **The geometric factor base is already optimal for the pair test.**
   - Every random `V` sampled at `n = 13, 19, 23`, `l = 3–5` has `dim span(V·V) ≥ 2l − 1`.
     Geometric `V` reaches exactly `2l − 1`.
   - This is the linear analogue of Kneser's theorem (Hou–Leung–Xiang). With no
     intermediate subfield it reads `dim(V·W) ≥ dim V + dim W − 1`.
   - So the pair oracle needs at least `3l − 1` unknowns for every `V`, not only geometric
     ones. Its reach `l ≤ (n+1)/3` is a ceiling rather than an artefact of our choice of `V`.
3. **The cost model sets the target for any four-summand test.**
   - These are the per-attempt costs, in group operations, that still fit the rho budget
     under the [`gcost.py`](../geometric_v_linear_20261006/gcost.py) heuristics:

     | test resolves | `m = 4` | `m = 5` |
     |---|---|---|
     | same-conjugate decompositions only | `≤ 2^11.0` at `c = 29` | `≤ 2^32.5` at `c = 28` |
     | mixed Frobenius conjugates too | `≤ 2^33.8` at `c = 27` | `≤ 2^39.8` at `c = 21` |

   - Mixed conjugates are what the linear oracle cannot do cheaply. For geometric `V` at
     `n = 19, l = 5`, `dim(V·V^{2^k})` is 13–19 against `2l − 1 = 9`.
   - So W2 and W3 meet at one question: **is there a test for mixed-conjugate 4- or
     5-summand decompositions that costs at most about `2^25`–`2^34` per attempt?**

## Workstreams

### W1. Write up the ceiling (analytic; do first)

**Goal.** Turn the closures into one scoped boundary result. The other two workstreams
then cite it for their kill criteria.

**Claims to state, prove and scope:**

- **C1 (reach).** For any `F₂`-subspace `V ⊂ F_{2^n}` with `n` prime, the `(A, P)`
  linearization of `S3` needs at least `3·dim V − 1` unknowns. So a unique-solution pair
  test exists only while `dim V ≤ (n+1)/3`.
  - Proof: the Hou–Leung–Xiang bound with no intermediate subfield.
  - The equality case (linear Vosper, Bachoc–Serra–Zémor) characterizes geometric
    progressions. That is why geometric `V` is the optimum.
- **C2 (test count).** Take any method that finds relations by testing group elements
  with a pair test, so each test succeeds with probability at most about `2^{2l−n}`. It
  needs about `2^{n−l}` tests for `2^l` relations. That beats `√r` only if `l > n/2`, and
  C1 forbids `l > n/2`.
  - So every pair-test method over a subspace factor base is rho-dominated at prime `n`.
  - This covers the linear-algebra route through the `c ≈ 30` cap, and the cycle route
    measured in #1510.
- **C3 (transfers).** Weil descent (GHS) needs an intermediate subfield, and 131 is prime.
  MOV and other pairing transfers need a small embedding degree. Record both as closed,
  with citations.

**Scope stated explicitly.** C1 and C2 say nothing about:
- tests that are not linearizations;
- `m ≥ 3`-summand tests that do not reduce to pair tests;
- factor bases that are not subspaces;
- tests that resolve mixed Frobenius conjugates.

Those are exactly W2 and W3.

**Verification.** The frozen protocol checks C1's bound numerically: random and
geometric `V`, prime `n ∈ {13, 17, 19, 23, 29, 31}`, `l = 3–⌊(n+1)/3⌋`.
- **Prediction:** no violation of `dim(V·V) ≥ 2l − 1`, and equality for geometric `V`.
- **Falsifier:** a single violation means the citation is misapplied.
- Before freezing, check the bound against the papers themselves, not the summary above.

**Deliverables.**
- `research/decomposition_ceiling_<date>/`, containing `NOTE.md` (statements, proofs,
  scope, citations) and the results of a Rust verifier.
- An `identity.certificate/v1` record for C1's bound
  (`docs/identity-certificates/README.md`), cited by its `IDC1h…` id.
- The report, diagram and PDF the `research-visuals` skill requires.
- A `knowledge/` finding through `/curate-knowledge` once the note is validated.

**Cost.** Small: the Rust verifier, a few CPU-minutes of runs, plus the writing.

### W2. Frobenius-orbit factor base (analytic first, then toy measurement)

**Corrected question.** Fact 1 rules out a closed `V`. Instead, use `V` together with its
131 conjugates `τ^j(V)`, with `log(τ^j P) = λ^j log P`. That gives 131 times the factor
base for the same unknowns, but the gain is real only if the test sees mixed pairs
`τ^i A + τ^j B`.
- That test needs `dim span(V·V^{2^k})` small for some `k ≠ 0`.

**Analytic prediction.** Equality in the Kneser-type bound for `V` and `W = V^{2^k}`
needs progressions with the same ratio. That means `g^{2^k} = g^{±1}`, which forces
`g ∈ F_{2^{gcd(k,131)}} = F₂`. So the minimum is never reached.
- The open question is how large the excess `e(l) = min over V of dim(V·V^{2^k}) − (2l − 1)`
  is.
- If `e(l)` stays bounded by a constant, the mixed pair test survives and the orbit row of
  the cost map applies.
- If `e(l)` grows with `l`, the branch is closed.

**Experiment.**
- **Fields:** `n ∈ {13, 19}`, where 2 is a primitive root, as at 131.
- **Positive control:** `n = 23`. A `V` inside its 11-dimensional Frobenius-invariant
  subspace must show `e = 0` for the detector to count.
- **Families of `V`:**
  - geometric, every ratio `g` exhaustively at `n = 13`;
  - normal-basis-aligned;
  - random;
  - a small hill-climb that minimizes `dim(V·V^{2^k})`.
- **Measure:** `e(l)` for `l = 3–6` and every `k`.
- **Kill:** `e(l) ≥ l − 2` for every family at both fields. The mixed test's reach then
  falls below the `l` where the orbit row beats the plain row.
- **Reopen:** `e(l) ≤ 2` with no trend in `l`. Recompute the orbit cost map with the
  measured `e`, then carry the result into W3.

**Cost.** Small: exhaustive at `n = 13`, sampled at 19, all CPU-minutes.

### W3. A four- or five-summand test (bounded search with a hard kill)

**Target, frozen from the cost model before any search.** A test that decides `m = 4` or
`m = 5` decomposability of one target, within the per-attempt budget in fact 3 for the
conjugate model it handles.
- Same-conjugate only: `≤ 2^11` for `m = 4`, at `c ≈ 29`, where the linear algebra is
  already about `2^60`.
- Mixed conjugates: `≤ 2^25`–`2^34`.
- The protocol rederives these figures from `gcost.py` and freezes them. It does not copy
  them from this plan.

**Negative control, known to fail.** Glue two pair oracles by collision:
`R − (P₁ + P₂) ∈ pair set`, tested with the linear pair oracle.
- It costs about `2^{n−2c}` tests per relation, the same as the pair oracle on fresh
  targets. So it must come out at least as costly as the pair oracle.
- If the harness shows it doing better, the harness is wrong.

**Candidate generation.** Run one `/propose-ideas` round, through the idea-generator
subagent, seeded with:
1. The `S4 = G²` square-root structure, with `G = D² + D√σ₄ + σ₁σ₄ + σ₃` of lower degree
   than `S4`. Use it as a test of a single condition rather than a linearization (already
   closed at `l ≤ 5`).
2. λ-coordinate halving trees: `λ² + λ = x(2P)` is `F₂`-affine, so halving is a linear
   solve. Look for decompositions of `R` along a halving chain.
3. Trace and linear-part prefilters on `S5`. These are cheap necessary conditions that
   reject most targets before any solve. They reduce cost per test, not the number of
   tests, so they help only if the filtered test lands under budget.
4. Mixed-conjugate tests that use `λ`, the Frobenius eigenvalue on `⟨P⟩`, to align
   conjugates. Use this only if W2 reopens.

**Per-candidate gate (pre-compute falsification).** Before any code is written:
- an analytic unknown count, or a cost bound at `n = 131`;
- then a toy measurement at `n ∈ {13, 19}` of cost per attempt against `l`.

A candidate is killed as soon as either one exceeds the frozen budget's slope.

**Success short of rho.** A candidate that does not reach the rho budget still counts
under the success criterion if it lowers the measured charged gap below the current
55.22 bits. It is landed with its distance to rho.

**Branch kill.** No candidate passes the gate within the round, capped at 6 candidates,
and none lowers the baseline gap. The decomposition branch is then recorded as closed
under known test families, with W1 as its boundary.

**Cost.** Moderate: one ideation round, up to 6 analytic gates, and toy runs only for
candidates that pass.

### W4. Rho constants (outside decomposition; optional, independent)

Matched rho already uses the negation and Frobenius speedups. What is left is constants
only:
- step cost in `F_{2^131}` arithmetic;
- the distinguished-point fraction;
- batch inversion.

Measure them with `/ecbench-measure` and seal the result with `/ecbench-bounds`. For a
second host, use `/ecbench-independent-runner`.

**Claim tier:** constant-level only. Nothing here changes an exponent.

## Order and gates

1. **W1 first.** C1–C3 are cited by every other kill criterion.
2. **W2, and the budget freeze of W3, alongside W1.** Both are cheap and analytic.
3. **The rest of W3 after W2.** W2 decides which budget row applies, and whether seed 4
   is in scope.
4. **W4 at any time.** It doesn't interact with the others.

Rules for every experiment:
- its own directory;
- `PROTOCOL.md` frozen before any cited run;
- smoke tests disclosed and not cited;
- failures kept, corrections additive;
- for any comparative claim, negative and positive controls per the `controlled-experiment`
  skill;
- one PR per experiment.

## Files here

| file | role |
|---|---|
| `PLAN.md` | this plan |
| `reports/plan-20261016/` | the plan report: workstream diagram (SVG and its source) and PDF |
