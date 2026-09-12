# Stage 39: fixed algebraic factor-base holdout

Stage 39 removes the empirical target census from the n=41 setup. For the public parameters `n = 41`, `ell = 6`, and collection `m = 3`, it defines the abscissa set from the standard coordinate span

`W_6 = span_F2(1, x, x^2, x^3, x^4, x^5)`,

takes its Frobenius union, and closes the resulting point set under translation by the rational two-torsion point. In the workflow recipe this is `two_torsion_saturated(frobenius_union([1, 2, 4, 8, 16, 32]))`. The rule depends only on the curve representation and `ell`. It samples no subgroup target, enumerates no subgroup, and receives no discrete-log label. Constructing and materializing the predicate remains inside the measured select stage.

The holdout changes both random streams used by the tuned Stage 35 and Stage 38 runs. Relation collection uses seed 41331. The five public hash-to-curve targets use seeds 41301 through 41305, and none has a scalar constructed or supplied. Factor-base logarithms must again be derived from collected relations and certified by group multiplication.

The hosted Stage 38 sweep selected a 149-point window by the stages the parameter changes. It reduced collection plus log solving from 10.007 to 3.361 seconds, a 2.98-fold speedup, and reduced the charged third-summand lookups from 77,971,456 to 12,206,080. Collection starts with four 8,192-probe units and may extend to 64. The run retains the two-summand walked descent, sparse relation solve, same-target signed-Frobenius rho, fresh four-core build, single-CPU scientific process, process CPU time, wall time, and sampled process-tree memory.

A local same-binary tuning stream, relation seed 41131 and targets 41101 through 41105, preferred a 595-point window; that disagreement with the hosted 149-point winner is why this stage uses a new stream. At window 595 the local cell completed with 29 group-certified factor-base logarithms: algebraic construction took 0.217 seconds, total precomputation 1.022 seconds, five descents 0.284 seconds, and matched rho 0.275 seconds. The resulting algorithm-only amortized IC/rho ratio was 4.76 and the online ratio was 1.03. The production seeds above were chosen only after retiring that stream; the Linux singleton run is the untouched holdout for window 149.

This stage tests a setup constant and the stricter algebraic-definition boundary. It does not supply licensed Magma measurements, external reproduction, a full-cost crossover, a changed asymptotic exponent, or a Koblitz index-calculus SOTA result.

The first hosted attempt, run `34699355146`, completed the computation but failed independent custody replay because `predecessor_stage38.identity.path` stored the runner's absolute checkout path. Its relation seed 41231 and target seeds 41201 through 41205 are retired. The sealed output remains `completed_invalid` for scientific admission: its timings are not promoted here. The retry stores the same predecessor bytes under a repository-relative path and uses the untouched 413xx streams above.
