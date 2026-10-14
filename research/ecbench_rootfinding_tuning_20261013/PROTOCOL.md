# Prime two-sum table: frozen factor-base-size trial

Frozen on 2026-10-13 before execution. This is an incremental parameter
optimization of the exact two-summand root-finding strategy assessed in
`../ecbench_rootfinding_assessment_20261012/REPORT.md`, not a new solver or
an asymptotic claim.

## Question and prediction

On the same 20-bit prime curve, does a 64-abscissa pair table reduce the cold,
single-target *counted group-addition equivalent* cost relative to the prior
32-abscissa pair table? The previous eight-target cohort used 32 columns and
measured means of 1,056 table-setup GAE, 19,659 relation GAE, and 849
verification GAE. With `B` abscissae, a negation-folded table takes roughly
`B(B+1)` setup additions, while the idealized relation work falls roughly as
`r/B` for subgroup order `r=562427`. Balancing `B²+r/B` predicts a size near
`(r/2)^(1/3) ≈ 65`, so **64** is chosen before seeing the new targets.

The fixed reference is `rho.negation` on the same public point, and the A/A
control repeats it. The earlier direct-subtraction arm is not part of this
size comparison; its existing receipt remains a separate reference. The
generic rho floor and the 2026-10-12 rho observation remain boundary context.

## Inputs, controls, accounting

The sole spec is `specs/prime20-b32-vs-b64.json`, with the same deterministic
curve search (`bits=20, seed=59297`) as before and **eight fresh public
targets** from target seed `20261013`. Each arm gets identical target and
algorithm seeds. Both IC arms use the same exact `mitm:m=2,negation_folded=1`
oracle and `incremental-gauss-full-rank-checked` matrix policy; their factor
bases intentionally differ only in `size=32` versus `size=64`. Each has three
measured rounds after one warmup, alternating order, `max_trials=5000000`,
and a 120-second timeout per execution. All measured outcomes, including
failures, remain in the session.

Primary output: paired cold GAE per verified target, its exact ratio,
factor-base columns and digest, table setup, relation trials/work, matrix,
verification, and peak RSS. Report the matched rho rows separately. The
comparison succeeds on the counted metric if both IC arms verify all 24
measured executions, the new size has lower paired total GAE on the cohort,
and full replay agrees. Otherwise report the observed ratio and failures.
No subsequent size may be selected using this holdout cohort.

Prime factor-base lifts, Legendre tests, square roots, map lookups, modular
row operations, and some rho native work are unpriced by ecbench. Therefore
the counted comparison remains a **lower-bound operation comparison**, even
if the ratio improves. macOS runs are L0; wall values are descriptive only.
No complete-cost runtime, larger-field, or asymptotic conclusion follows.

## Execution and stop

Use the existing release ecbench binary after recording its hash and source
revision. Run the spec once into a new session directory; do not overwrite
or edit the 2026-10-12 sessions. Audit with `ecbench audit --replay-all`, save
paired comparisons, and stop after those 128 planned executions or a runner
failure. Preserve the exact command, host capsule, raw records, audit, and
visual report. Do not change the spec after an outcome is visible.
