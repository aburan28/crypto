# Enlarging Smith's 18.57% by isogeny walks: only the rationality conditions can move

- **Status:** `open`
- **Kind:** implementation gap; empirical bottleneck (walk cost unknown)
- **Affected work:** [weak isogeny classes study, 2026-10-13](../research/weak_isogeny_classes_20261013/REPORT.md),
  lead 2
- **Next action owner:** unassigned
- **Last updated:** 2026-10-13

## Problem

We cannot yet say how much of the gap between Smith's 18.57% and the 50.16%
ceiling an isogeny walk inside the hyperelliptic locus recovers. The
per-subgroup rationality conditions (Smith's Proposition 2: a discriminant
is a square; Proposition 4: the leading coefficient of `s(t)` is a square) are
not implemented here. Neither is any explicit hyperelliptic-to-hyperelliptic
isogeny inside the class.

## Evidence and uncertainty

- **Observed:** `smith_kernel_coverage` (exhaustive over every degree-8
  model for `p = 3, 5, 7, 11`; see the study's `data/smith_p*.json`).
  Classes are `L`-polynomials. In every class, either all hyperelliptic
  members have a rational tractable subgroup or none does. The class closure
  equals the direct fraction.
- **Derived (finite check, `tests::reproduces_smith_lemma1_and_constant`):**
  on all 22 cycle types of `S_8`, the characteristic polynomial of Frobenius
  on `J[2]`, which is `chi_pi mod 2` and an isogeny invariant, determines
  whether a rational tractable subgroup exists. So walking inside the
  hyperelliptic locus never creates a tractable subgroup. It can only change
  the outcome of the two square conditions, which Smith models as
  independent fair coins (Hypotheses 2 and 3). Under those hypotheses the
  asymptotic ceiling is `P(#S(H) > 0) = 0.50159` against Smith's `0.18574`.
- **Missing capability:** native Proposition 2 and 4 tests. Explicit
  `(2,2,2)` or `(l,l,l)` isogenies between hyperelliptic genus-3 Jacobians,
  with their cost.
- **What this does not establish:** that the ceiling is reached, or the cost
  of reaching it. It also leaves out walks through non-hyperelliptic members
  of the class.

## What would unblock it

- [ ] Native Proposition 2 and 4 tests, validated by reproducing Smith's
  experiment (18.58% success on random curves; Smith §8).
- [ ] Exhaustive per-class success at `p = 7, 11`: the share of curves whose
  class has a member where some tractable subgroup passes both tests.

## Next discriminating step

- **Action:** implement Smith's matrix `M` (eq. (2)), the discriminant (4)
  and `s(t)` (Lemma 4) in Rust. Add a `full_success` column to `PrimeReport`.
- **Outcomes:** if the class-level success is near 50%, the walk extension is
  worth costing. If it is near 18.6%, the conditions are correlated within a
  class and the extension gains little.
- **Recheck trigger:** that implementation landing.

## Resolution and history

- 2026-10-13: opened with the coverage study. Unresolved.
