# Protocol I-4: does index-calculus relation yield vary between isogenous curves?

Frozen 2026-10-08, before any instrument is built.  **Boundary; an advance
only if variation survives every control, and then only as a candidate
for its own protocol at the next size.**  `S`, end-to-end cost and
speedup are **unset**.

## Derivation (stated before measuring)

Index calculus reads the equation: the factor base is a set of
abscissae, and the decomposition oracle solves `S_{m+1}` on the chosen
model.  Transport across `φ` sends the factor base of `E` to a set of
points of `E′` that is *not* a factor base of the same family on `E′`
(a subspace of abscissae does not map to a subspace of abscissae), so
nothing forces the relation yield `γ` of a fixed factor-base family to
agree between isogenous curves.  Two hypotheses:

- **H4-null.**  For a fixed family (standard subspace of dimension `ℓ_fb`
  over binary fields; the leaderboard's prime-field factor base for
  `F_p`), the yield `γ(E′)` and the whole-pipeline `S_IC(E′)` across a
  class lie within the A/A interval of a single curve.  Mechanism: the
  product law of `RESEARCH_ECC2K130_DECOMPOSITION.md` §5 is a counting
  statement about `m`-sums of a `2^{ℓ_fb}`-point set in a group of order
  `r`, and the model enters only through the oracle's constant.
- **H4-alt.**  `γ` varies with the model beyond the A/A interval.
  Mechanism candidates, each falsifiable: coefficients in a subfield
  align the factor base with a Frobenius-stable subspace (binary; test by
  the I-3 screen); `j` small or in a subfield changes the summation
  polynomial's symmetric structure; a vertical step changes `End(E)` and
  with it the endomorphisms available for factor-base folding.

If H4-alt held on prime fields with none of the named mechanisms, it
would be the first model-dependent prime-field IC effect in this
repository and would be classed as a **reproducible unexplained
anomaly** until a mechanism is named and the effect reproduced on fresh
classes.

## Instrument (Rust, follow-on PR)

- The native `ic` pipeline (`docs/ic/README.md`) with `--solver sat` and
  `--solver enumerate` as the independent cross-check, one target per
  run, the five IC online phases charged, `Q = [d]P` verified on the
  original `E` after transport back.
- Per class: arms on `E`, two small-`ℓ` neighbours, one large-`ℓ`
  neighbour, one vertical neighbour where one exists, and the
  isomorphic-model control of one neighbour; the same factor-base family
  and dimension on every arm; `ecbench` spec with strong rho on `E` as the
  matched reference.
- Yield is measured as relations per trial over a fixed trial count, with
  the ceiling of the product law beside it.

## Frozen inputs

| item | value |
|:--|:--|
| binary | `icv1-f2m13-t181-515ee569` (`ℓ_fb = 5`, `m = 3`), `icv1-f2m19-t797-b6cf2467` (`ℓ_fb = 7`), one composite class over `F_{2^{15}}` from I-3 |
| prime field | the I-1 classes at `p ≈ 2^{20}` and `2^{24}` with the leaderboard's factor base at the size's registered dimension |
| trials per arm | 2,000 at `2^{13}`/`2^{20}`; 500 at `2^{19}`/`2^{24}` |
| targets per arm | 8, hidden scalars |
| seed | 20261014 |
| unit | ecbench group operations; `S_IC` whole pipeline; `γ` relations per trial |

## Predictions (pass/fail)

- **Q1 (null).**  For every class, the paired bootstrap 95% interval of
  `γ(E′)/γ(E)` contains 1 for every arm, and so does that of
  `S_IC(E′)/S_IC(E)`.
- **Q2 (ceiling).**  Every `γ` is at or below the product-law ceiling for
  its `(ℓ_fb, m, r)`; a yield above the ceiling is a bug.
- **Q3 (oracle agreement).**  `sat` and `enumerate` agree on every
  decomposition verdict on every arm.
- **Q4 (positive control).**  On the composite binary class, the arm with
  the smallest `m(b)` from I-3 shows a `γ` outside the A/A interval
  *only if* the factor base is the descent-aligned one; with the standard
  subspace it does not.  (The GHS mechanism acts through the descent, not
  through the subspace factor base.)

## Decision rule (registered)

- Q1–Q4 pass: **boundary**.  Relation yield of a fixed family is
  class-invariant at these sizes; the model enters only through the
  oracle constant, and the weak-curve methodology's `N(E′)` may be taken
  as constant across a class for this family.
- Q1 fails on an arm: run the controls in order — isomorphic model, second
  field representation, `enumerate` oracle, a fresh class of the same
  shape — and classify: an effect that dies under any control is an
  implementation effect or a statistical artifact; one that survives all
  four is a **reproducible unexplained anomaly**, recorded with its
  effect size, BH-adjusted confidence and the mechanism candidates ruled
  out, and it earns a protocol at the next size.  It is not an advance:
  no exponent moved.
- Q2 or Q3 fails: a bug; halt.

## Stop condition and inadmissible moves

Bounded: the classes above, the arm set above, one run plus the control
sequence on any failing arm.

Inadmissible: changing `ℓ_fb`, `m`, the factor-base family or the trial
count between arms; reporting a yield without its ceiling; reading a
solver-stage difference as a yield difference; a claim at `p ≤ 2^{24}` or
`n ≤ 19` about any registered target curve.
