// Findings. Numbers from data/summary_p*.json; derived quantities show their formula.
#let s7 = json("/data/summary_p7.json")
#let s11 = json("/data/summary_p11.json")
#let f(x, n: 4) = str(calc.round(x, digits: n))
#let closed(s) = (s.U - s.null_primary.median) / (s.U - s.C)

= Findings
- *Measured.* `full_success` is an isomorphism invariant on every census model (control C2), and the census models
  fall into #s7.iso_classes isomorphism classes at $p = 7$ and #s11.iso_classes at $p = 11$; a tractable $L$-class
  holds a median of #s7.median_tractable_iso_classes_per_tractable_l_class and
  #s11.median_tractable_iso_classes_per_tractable_l_class tractable isomorphism classes.
- *Measured.* With the isomorphism class as the unit, the null median of $C$ drops from $U$ to
  #f(s7.null_primary.median) ($p = 7$) and #f(s11.null_primary.median) ($p = 11$). *Derived*,
  $(U - "null median") / (U - C)$: counting curves instead of models accounts for #f(closed(s7) * 100, n: 0)% and
  #f(closed(s11) * 100, n: 0)% of the gap between $U$ and the observed $C$.
- *Measured, frozen rule.* The observed $C$ is below every one of the 999 permutations at both primes, also when
  permutations are stratified by orbit size: *#s7.verdict* at $p = 7$ and $p = 11$. Tractable $L$-classes with no
  `full_success` member: #s7.tractable_l_classes_without_full.observed observed against
  #f(s7.tractable_l_classes_without_full.expected_independent_iso, n: 0) expected ($p = 7$),
  #s11.tractable_l_classes_without_full.observed against
  #f(s11.tractable_l_classes_without_full.expected_independent_iso, n: 0) ($p = 11$). The excess appears at every
  $k$ (Figure 1b), so it is not confined to the smallest classes.
- *Post hoc, exploratory (not registered).* Within each `#S(H)` stratum the isomorphism-level success rate barely
  moves with the 2-adic valuation of $\#J(FF_p) = L(1)$: at $p = 11$ it stays within 0.239–0.271 for
  $\#S(H) = 1$ and 0.549–0.584 for $\#S(H) = 3$ (`post_hoc_by_stratum_and_v2_of_J` in the summaries). That is far
  too small to double the number of empty classes, so the 2-power part of the group order is not the hidden
  class invariant.
- *Interpretation (labelled).* Successes cluster within $L$-classes: the isomorphism classes of one $L$-class are
  positively correlated in whether Smith's two square conditions pass. Something finer than $chi_pi$ mod 2 and
  coarser than the curve is shared within a class; a candidate not tested here is the Galois structure of
  $J[4]$ or of the 2-adic Tate module, which the $L$-polynomial does not fix. For the walk lead this means that
  the reachable set is $C$, not the granularity-corrected null, and that the classes with no `full_success`
  member are a real feature of the data rather than an artefact of how models were counted.
