# Fixed quadratic state for complete generic Boolean search

This standalone experiment continues the coefficient-aware construction studies by
measuring a complete, bounded Boolean solve. It accepts only generated synthetic
systems through the benchmark CLI. It has no curve or target interface and does
not modify the repository's production solver.

## Contract before measurement

Work in the Boolean ring over GF(2). Write an equation as

\[
f(x)=c+\sum_i l_i x_i+\sum_{i<j}q_{ij}x_ix_j.
\]

For the assignment \(x_j=b\), retain every \(q_{ik}\) between unassigned
variables and update

\[
c' = c\oplus b l_j,\qquad
l_i' = l_i\oplus b q_{ij}\quad(i\ne j).
\]

Remove variable j and its incident quadratic terms. Apply a batch of assignments
in increasing variable order. This also charges contributions of edges whose two
endpoints are assigned: the first endpoint changes a linear coefficient and the
second then changes the constant. XOR explicitly handles coefficient cancellation.
An exact residual edge count distinguishes quadratic, affine, constant and zero
equations, including degree drops after assignment.

The immutable table stores both directions of each quadratic edge. Mutable state
stores the active-variable mask, live-equation mask, residual linear coefficients,
constant coefficients and edge counts. Exact deduplication compares residual
coefficients and restricted edge sets, keeping the first equation in original
order. It never accepts hash equality as equation equality. Initial equations,
including duplicates and zeros, are preserved until the same specialization point
as the retained implementation.

This is a restricted quadratic-state experiment, not a general Macaulay schedule.
It generates no polynomial multiples. The earlier support-envelope experiment
handles newly required multipliers after generator degree drops; its result and
contract remain in `../boolean_support_envelope_20260922`. Here a degree drop is
handled directly by the existing affine propagation policy. A non-quadratic or
otherwise unsupported source takes the retained general search path. The fixed
benchmark contains only canonical quadratic systems, so that fallback is tested
but never used for its timing advantage.

The complete search policy remains unchanged: exact affine reduction, forced
assignments, most-frequent-variable selection with low-index ties, zero branch
first, and the same node cap. The trace consumer receives the same canonical
monomials in the same order without building a temporary list. Search outcomes,
models, logical counters and trace digests must match both retained search arms.
The digest is a diagnostic; exhaustive residual comparisons and independent model
evaluation supply separate correctness checks.

## Fixed reference and rejection criterion

`protocol.json` fixes four sizes, three families, two retained discovery seeds,
two fresh holdout seeds and eight rotated repetitions of all eight arms. The
reference is the pointwise fastest of **all seven** prior complete-solve methods,
including `merge_search`. The new representation must exceed 2x with a 95% paired
bootstrap lower bound on **every** family at 16, 20 and 24 variables. All nine
comparisons and all grid completions/verification are required. Timing intervals
describe the fixed fixtures and host; two holdout fixtures per group do not
establish performance across a population of systems.

Compilation, state construction, the entire search, context destruction and
result checking are charged to each cold solve. Fixture generation and preparation
of the completed reference status sit outside arm timing but inside the worker
receipt. The no-kernel arms do zero matrix work. Logical counters are invariant
work diagnostics, not calibrated operation counts or hardware-independent speed
claims. Fresh-worker RSS includes all arms; candidate-specific peak allocation
is not measured. At most 36 variables/equations fit the candidate representation.

UNKNOWN is censored and has no completion cost. No timing-only advantage can
substitute for missing verification. Full index-calculus, production-solver and
rho costs remain null. The WDSat/full-curve suite is inapplicable to this generic
standalone experiment, and no claim about those pipelines is made.

## Reproduction

```sh
python3 research/boolean_quadratic_state_20260923/run.py \
  --out research/boolean_quadratic_state_20260923/run_01
python3 research/boolean_quadratic_state_20260923/run.py \
  --protocol research/boolean_quadratic_state_20260923/confirmation_protocol.json \
  --out research/boolean_quadratic_state_20260923/run_02
python3 -m unittest discover \
  -s research/boolean_quadratic_state_20260923 -p 'test_*.py'
```

The runner requires a new output directory, freezes its source and protocol,
compiles and runs correctness tests before any measurements, retains process
receipts and hashes, and writes the complete comparison. Use another new directory
for reproduction. Never overwrite a retained run or retune against its holdouts.

`SOURCE_LINEAGE.json` pins the predecessor controls. Each new run also retains
the exact modified source, protocol, analysis, compiler and host metadata.

The confirmation protocol changes only the two held-out seeds and adds provenance
for the primary run. It uses the unchanged candidate, controls and acceptance gate.
The two runs share their 24 discovery fixtures; repetitions and repeated discovery
fixtures must not be counted as new independently generated systems.
