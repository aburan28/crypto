// Findings. Numbers from data/summary_p*.json and data/controls_p*.json.
#let s7 = json("/data/summary_p7.json")
#let s11 = json("/data/summary_p11.json")
#let f(x, n: 4) = str(calc.round(x, digits: n))
#let g(s, k) = s.candidates.at(k).gap_share * 100

= Findings
- *Measured, frozen rule.* No candidate explains the within-class clustering. Conditioning the iso-level null on
  the cycle type (X1), the 2-Sylow subgroup of $J(FF_p)$ (X2) or the corrected $J[4]$ descent invariants (X3c, X3c+)
  leaves the observed $C$ below every one of the 999 permutations at both primes. The share of the gap
  "unconditioned null median $-$ $C$" that conditioning removes is at most
  #f(calc.max(g(s7, "X2_sylow2"), g(s11, "X2_sylow2")), n: 0)% (X2) and
  #f(calc.max(g(s7, "X3c_j4_profile_corrected"), g(s11, "X3c_j4_profile_corrected")), n: 0)% (X3c);
  for X1 it is negative (#f(g(s7, "X1_cycle_type"), n: 1)% / #f(g(s11, "X1_cycle_type"), n: 1)%). The verdict
  for every valid candidate is *does not explain* at $p = 7$ and $p = 11$.
- *Measured.* The invariants do move the isomorphism-level success rate: within $\#S(H) = 1$ it ranges over
  0.207–0.300 ($p = 7$) and 0.223–0.287 ($p = 11$) across cycle types, and 0.228–0.299 / 0.212–0.312 across 2-Sylow
  types with at least 200 classes (`descriptive_by_stratum_and_*` in the summaries). Rates that vary *between*
  cells but are shared by $L$-classes of every kind redistribute successes without concentrating them in
  particular $L$-classes, which is why the conditioned nulls barely move.
- *Measured, frozen rule.* X3 and X3+ (the frozen descent instrument) are *invalid*: control K4 failed on
  #s7.controls.counts.k4_mismatch ($p = 7$) and #s11.controls.counts.k4_mismatch ($p = 11$) classes, all of
  all-even cycle type, which Amendment 1 explains as the order-2 kernel of the even-degree descent map. Their nulls
  are recorded for completeness and are as far from $C$ as the valid ones.
- *Derived, not a verdict.* X3c leaves #s7.controls.counts.x3c_classes_with_ambiguous_k /
  #s11.controls.counts.x3c_classes_with_ambiguous_k classes with an undecided $k > 1$ entry. Such classes are cells
  of their own value, so the invariant is, if anything, finer than the true $J[4]$ profile; a finer conditioning
  can only shrink the null toward $C$, and it did not.
- *Interpretation (labelled hypothesis).* The clustering is not carried by the 2-adic Galois module of $J_H$ at the
  levels tested. A remaining structural candidate is the target of Smith's isogeny itself: `full_success` asks
  whether some $(2,2,2)$-quotient $J_H\/S$ with a rational trigonal map is a Jacobian over $FF_p$ rather than its
  quadratic twist (Smith Remark 4; Serre's obstruction). The quotients lie in the same $L$-class, so curves of one
  class share candidate quotients, which would correlate their outcomes. Testing this needs the isomorphism class
  of the quotient curve $X$ per subgroup (e.g. Dixmier–Ohno invariants), not a further invariant of $J_H$.
