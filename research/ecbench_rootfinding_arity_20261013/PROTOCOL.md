# Prime pair-table arity at 64 abscissae: frozen trial

Frozen before execution on 2026-10-13. This follows the 32-versus-64
two-summand size trial, whose eight targets are excluded here.

## Question and prediction

At `prime-abscissa:size=64` on the same registered 20-bit curve, does
three-summand collection lower cold single-target counted group-addition
equivalents (GAE) relative to the two-summand pair-table configuration?
The three-summand probe scans signed base points against a pair table;
its larger per-attempt cost may be offset by higher relation yield. The
candidate is selected before observing the fresh targets.

Compare `mitm:m=2,negation_folded=1` with
`mitm:m=3,negation_folded=1` at identical base size, matrix policy,
target points, and algorithm seeds. A secondary arm runs
`rootfind-subtract:m=3`: it follows the signed subtraction sweep with
batched field inversions and the same folded pair table. Its native field
work is not fully priced in GAE, so its GAE is descriptive, not a runtime
comparison. `rho.negation` is the matched reference and repeats as A/A.

## Inputs and accounting

The sole spec is `specs/prime20-b64-m2-vs-m3.json`. Use the registered
curve `icv1-fp20-t727-cd198a38` (prime subgroup order 562427), generated
by `prime_search` bits 20, seed 59297. Draw eight fresh public targets with
target seed 20261014. Run three measured rounds after one warmup, alternate
arms, one thread, `max_trials=5000000`, 120-second timeout per execution.
Preserve every failure or timeout. The planned run has 160 executions.

Primary metric: paired cold total GAE for the two `mitm` arms, with setup,
relation collection, matrix, and verification costs separately reported.
Report factor-base digest/columns, trial counts, rank, memory, and matched
rho independently. The candidate improves the counted metric only if all
24 measured runs in each primary arm verify, its paired total GAE is lower,
and full replay agrees. Otherwise report the exact observed ratio and
failures. Do not retune using these targets.

Prime factor-base lifts, Legendre tests, square roots, hash probes and row
operations, plus some rho native work, remain unpriced in ecbench.gae.
The batch-inversion arm additionally changes the field-operation mix. A
counted ratio is a partial cost comparison; no full-cost or wall-time
advantage follows. This macOS host runs at L0 unless the harness attests
otherwise. Wall observations below L2 are descriptive only.

## Execution and stop

Use the existing release ecbench binary, recording its SHA-256, source
revision, host capsule, and exact command. Run once into a new session
directory, then `ecbench verify --replay-all`, save paired comparisons and
a source-linked visual run report. Stop after the 160 planned executions or
runner failure; do not alter the frozen spec after any outcome is visible.
