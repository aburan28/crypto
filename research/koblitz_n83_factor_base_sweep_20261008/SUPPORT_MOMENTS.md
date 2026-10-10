# Exact support-count screen for the retained N83 bases

This is a mathematical feasibility screen for **full smooth relations**. It does not select a factor base or estimate an index-calculus runtime. The machine-readable exact integers and fractions are in `pilot-01/support-moments.json`, produced by `support_moments.py` from the verified panel and its retained point-object headers.

Let a base contain `B` distinct points in a prime-order subgroup of order `r`. If an `m`-summand relation permits repeated points, there are exactly `M = binomial(B+m-1,m)` unordered multisets of summands. Each multiset has one fixed group sum. For a target uniform over the whole subgroup, the expected number of multiset decompositions is exactly `M/r`, regardless of collisions, orbit structure or dependence between sums.

For a target uniform over the `r-1` nonidentity points, let `Z` be the number of multisets summing to identity. The expected decomposition count is `(M-Z)/(r-1)`, hence at most `M/(r-1)`. Markov's inequality gives the rigorous ceiling `Pr(at least one full smooth decomposition) <= min(1, M/(r-1))`. We do not compute `Z`; signed bases contain cancellation pairs, so assuming it is zero only weakens the upper bound. The ceiling is **not** a probability estimate for a fixed public fixture. It also says nothing about how long a solver takes to find an existing relation.

Every retained signed-Frobenius base has `B = 166K` distinct subgroup points. The table shows the rounded nonidentity-target hit ceiling for each retained `K`; `1` means the bound is vacuous. The JSON retains exact numerators and denominators.

| Arm | K | m=2 | m=3 | m=4 | m=5 | m=6 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| a=0, 81-bit | 64 | 2.33e-17 | 8.27e-14 | 2.20e-10 | 4.67e-7 | 8.27e-4 |
| a=0, 81-bit | 256 | 3.73e-16 | 5.29e-12 | 5.62e-8 | 4.78e-4 | 1 |
| a=0, 81-bit | 600 | 2.05e-15 | 6.81e-11 | 1.70e-6 | 3.38e-2 | 1 |
| a=1, 53-bit | 64 | 6.59e-9 | 2.33e-5 | 6.20e-2 | 1 | 1 |
| a=1, 53-bit | 256 | 1.05e-7 | 1.49e-3 | 1 | 1 | 1 |
| a=1, 53-bit | 600 | 5.79e-7 | 1.92e-2 | 1 | 1 | 1 |

For this `B = 166K` construction, the first integer `K` at which the Markov ceiling becomes vacuous is:

| Arm | m=2 | m=3 | m=4 | m=5 | m=6 |
| --- | ---: | ---: | ---: | ---: | ---: |
| a=0, 81-bit | 13,247,128,046 | 1,469,216 | 16,627 | 1,182 | 209 |
| a=1, 53-bit | 788,664 | 2,241 | 129 | 25 | 9 |

This threshold is necessary only for the upper bound to reach one; it is not sufficient for a relation, rank gain, or a complete logarithm. For the primary a=0 arm, the retained K=600 four-summand base has a full-smooth hit ceiling below 2e-6 under the uniform-target model. The six-summand ceiling is vacuous at K=256 and K=600, making higher arity a more plausible full-smooth direction at these K values, while its algebraic-solving cost remains unmeasured. The diagnostic a=1 K=64 four-summand ceiling is about 0.062; its capped cold run cannot be interpreted as a rate estimate.

The declared v1 grid stops at K=900, below the primary m=5 threshold K=1,182 where this ceiling first reaches one. A future grid that tests larger full-smooth five-summand bases needs a new frozen design version and must charge the correspondingly larger base, index, relation and matrix costs. The current m=6 options remain adapter-gated rather than measured winners.

The 1,728 executed two-summand probes use fixed public targets and repeated or nested bases. Their zero-relation result is compatible with the small uniform-target ceilings, but does not validate this target model or supply independent Bernoulli trials. Solver backends, Gray enumeration, splitting and domain-preserving symmetry change search cost or duplicate representations; they do not create additional full-smooth group sums for a fixed point set and `m`. Double-large-prime partials change the accepted relation event and require separate graph/rank accounting. Future experiments must test their ordinary yield directly before claiming a runtime benefit.

Verification: `python3 -m unittest discover -s research/koblitz_n83_factor_base_sweep_20261008 -p test_support_moments.py -v` checks multiset counting, the nonidentity-target formula, Markov's ceiling on a signed small group, and exact threshold minimality. The generator checks all 54 manifest row sizes, successful generic replay, and one retained point-object header per arm and size before writing the immutable JSON receipt.
