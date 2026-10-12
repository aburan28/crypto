# PS-1: algebra spec, bit-degree and root-test formulations

Part of [`../PLAN.md`](../PLAN.md); frozen when the plan merges.  Repo: crypto.
Needs PS-0.

**Hypotheses.**

- **H1.**  Writing the symmetrised `S₄` as a quartic in the last abscissa,
  every coefficient is a Boolean-vector polynomial in the bits of the other
  two summands of bit-degree at most 4 (and at most 2 for most), so a Gray
  step changes it by a table lookup of fixed depth.
- **H2.**  Some root-test formulation for "the quartic has a root in the
  subspace" costs fewer counted bit operations than `gcd(q, L_V mod q)`
  (candidates: the current gcd; a trace-based reduction to an affine
  quadratic; a resolvent-cubic test followed by a linear-map membership
  check).

**Frozen inputs.**  The repository's symmetrised `S₄` (`binary_semaev_s4`),
which is specialised to `b = 1` (the Koblitz family, `a ∈ {0, 1}`); general
`b` is out of scope and the module's panic on it stays.  Exhaustive checks for
`n ≤ 9` over all triples; for `n = 11, 13` over `10^6` random pairs against
the PS-0 baseline.

**Method.**  Compute the exact bit-degree of each coefficient symbolically;
state each candidate root test as a straight-line program over `F₂`; count
its gates; check equality with direct evaluation of `S₄` and with
`gcd(q, L_V mod q)`.  Every identity relied on gets an
`identity.certificate/v1` record under
[`docs/identity-certificates/`](../../../docs/identity-certificates/README.md).

**Unit.**  Gates per pair (XOR, AND, ternary-logic) for each formulation.

**Success.**  Bit-degrees and at least two root-test formulations certified;
exhaustive agreement at `n ≤ 9`; sampled agreement at `n = 11, 13`; the
cheaper formulation and its gate count named for PS-2.

**Stop and redirect.**  A coefficient of bit-degree 5 or more, or no
formulation under `1.8·10^5` gates per pair: PS-2 and PS-3 are re-scoped to
direct evaluation, and the result is recorded as a finding.

**Inadmissible.**  Claiming an identity without its certificate; testing
only on planted triples.

**Class.**  Accounting.  **Effort.**  3–4 days.
