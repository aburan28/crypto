# Enlarging Smith's 18.57% by isogeny walks: only the rationality conditions can move

- **Status:** `open` (the 2026-10-13 decision step is done: class-level success
  is near the ceiling, so the walk is worth costing; the walk itself is not
  built or priced)
- **Kind:** implementation gap; empirical bottleneck (walk cost unknown)
- **Affected work:** [weak isogeny classes study, 2026-10-13](../research/weak_isogeny_classes_20261013/REPORT.md),
  lead 2; [Smith Prop. 2/4 native study, 2026-10-14](../research/smith_prop24_native_20261014/REPORT.pdf)
- **Next action owner:** unassigned
- **Last updated:** 2026-10-15

## Problem

We cannot yet say what an isogeny walk inside the hyperelliptic locus costs
to move a curve to a member of its class where Smith's construction succeeds.
Since 2026-10-14 the per-subgroup rationality conditions (Smith's
Proposition 2: a discriminant is a square; Proposition 4: the leading
coefficient of `s(t)` is a square) are implemented natively
(`src/cryptanalysis/smith_prop24.rs`), and the share of curves whose class
holds such a member is measured exactly at `p = 7, 11`. No explicit
hyperelliptic-to-hyperelliptic isogeny inside the class is implemented, so
the walk's existence along explicit steps, its length and its cost remain
open.

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
- **Observed (2026-10-14, [study](../research/smith_prop24_native_20261014/REPORT.pdf),
  protocol frozen at `31d0beec8` before data):**
  - Question 1: on 10^6 random degree-8 models over `p = 2^31 - 1` the native
    tests reproduce Smith §8 (Example 2) within ±0.01 on all four fractions:
    0.50187 [0.50062, 0.50312] with a rational tractable subgroup, 0.49999
    of subgroups pass Prop. 2, 0.50102 of those pass Prop. 4, 0.18638
    [0.18541, 0.18736] of curves have `full_success` (Theorem 2: 0.18574;
    98.75% Clopper–Pearson intervals).
  - Question 2 (exhaustive, curve + twist, classes = `L`-polynomials): the
    share of curves whose class holds a `full_success` curve is
    `C = 0.4304` at `p = 7` and `0.4750` at `p = 11`, against `D = 0.1758`
    and `0.1805` directly and the ceiling `U = 0.4897` and `0.4943`;
    `r = (C − D)/(U − D) = 0.811` and `0.938`, both "near the ceiling" under
    the frozen 0.75 threshold. Controls passed: Smith's Example 1 reproduced
    (`N`, `D`, `δ_i`, `f_i`, `#H(F_37) = 42`); the 2026-10-13 counts
    reproduced; class closure of "has a rational tractable subgroup" equals
    its direct fraction; the point count of Recillas' curve `X` agreed with
    Prop. 4 in 18,381,760 cases and disagreed in 0.
  - Post hoc, descriptive: `C` is below the independence null (`C = U` in all
    199 permutations), because 560 of 1,855 (`p = 7`) and 1,210 of 7,431
    (`p = 11`) tractable classes hold no `full_success` curve; those classes
    are smaller (median 108 vs 348 models at `p = 7`, 420 vs 1,860 at
    `p = 11`). The two conditions are neither class invariants nor
    independent across class members.
- **Missing capability:** explicit `(2,2,2)` or `(l,l,l)` isogenies between
  hyperelliptic genus-3 Jacobians, with their cost. (Native Proposition 2 and 4
  tests: done 2026-10-14.)
- **What this does not establish:** that the ceiling is reached, or the cost
  of reaching it. It also leaves out walks through non-hyperelliptic members
  of the class.

## What would unblock it

- [x] Native Proposition 2 and 4 tests, validated by reproducing Smith's
  experiment (18.58% success on random curves; Smith §8). Done 2026-10-14:
  18.638% [18.541%, 18.736%] on 10^6 curves.
- [x] Exhaustive per-class success at `p = 7, 11`: the share of curves whose
  class has a member where some tractable subgroup passes both tests. Done
  2026-10-14: 43.04% and 47.50%.
- [ ] Explicit isogenies between hyperelliptic members of a class, and the
  walk distance from each tractable non-`full_success` curve to a
  `full_success` one at `p = 7, 11`, with the cost per step in field
  operations.

## Next discriminating step

- **Done (2026-10-14):** the step recorded on 2026-10-13 (implement `M`, the
  discriminant and `s(t)`; measure class-level success) ran under a frozen
  protocol. Outcome: near the ceiling at both primes, so the walk extension
  is worth costing.
- **Action now:** build explicit isogenies between hyperelliptic members of
  one class at `p = 7` (for example Smith's trigonal construction when `X`
  is hyperelliptic, or Richelot-type `(2,2,2)` steps), and measure on the
  exhaustive census the number of steps from each tractable
  non-`full_success` curve to a `full_success` curve, with field-operation
  counts per step. Controls: the class rows of
  `research/smith_prop24_native_20261014/data/` give the reachable set
  exactly, so a walk that reaches beyond `C` or misses a class it should reach
  is a bug.
- **Narrowed (2026-10-15,
  [study](../research/smith_iso_granularity_20261015/REPORT.pdf)):** why a
  minority of classes holds no `full_success` curve. Sparsity is not the
  answer: counting `F_p`-isomorphism classes of curves instead of models
  (34,227 at `p = 7`, 324,562 at `p = 11`; `full_success` checked constant on
  every model) moves the null median of `C` from `U` to 0.4579 / 0.4864, but
  the observed 0.4304 / 0.4750 stays below all 999 permutations at both
  primes (frozen verdict: residual deficit beyond granularity). Empty
  tractable classes: 560 vs 395 expected and 1,210 vs 720, in excess at every
  class size. Successes are positively correlated within an `L`-class; the
  2-adic valuation of `#J(F_p)` does not explain it (post hoc). Still open:
  which finer invariant (candidate: the Galois module `J[4]`) is shared; and Smith's Remark 5 alternative (descending `X` when the trigonal
  map is not rational, which would raise the direct rate to 31.13% under his
  hypotheses) as a route that needs no walk.
- **Outcomes:** short walks with cheap steps make the extension a practical
  enlargement of Smith's weak set; long or expensive walks bound it as a
  counting fact only.
- **Recheck trigger:** an explicit hyperelliptic-to-hyperelliptic genus-3
  isogeny landing in the repository.

## Resolution and history

- 2026-10-13: opened with the coverage study. Unresolved.
- 2026-10-14: native Proposition 2 and 4 tests landed; Smith §8 reproduced on
  10^6 curves; class-level success near the ceiling at `p = 7` (`r = 0.811`)
  and `p = 11` (`r = 0.938`). Unblock items 1 and 2 done; the walk's cost
  stays open.
- 2026-10-15: isomorphism-class granularity tested under a frozen protocol;
  it explains 54% (`p = 7`) and 41% (`p = 11`) of the gap `U − C`, and the
  rest is a within-class correlation not yet named.
