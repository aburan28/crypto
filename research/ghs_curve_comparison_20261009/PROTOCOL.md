# Binary-curve ECDLP comparison protocol (2026-10-09)

## Question and fixed inputs

Compare single-target ECDLP costs and binary GHS structure across registered
Koblitz models, keeping field representation, subgroup, and solver identity
explicit. This is a bounded comparison of existing curves and algorithms, not
a test of a complete higher-genus GHS index-calculus pipeline.

The frozen panel is `K_a/F_(2^n)` for `(a,n)=(0,7),(0,9),(1,9),(0,13),(0,15)`.
These five registered models give four increasing field degrees and a same-size
curve comparison at degree nine. Use the exact registry modulus for each model.
The ECDLP spec is `ecbench.spec.json`: eight independent targets per curve,
target seed 20261009, one round, no warmup, alternating order, and L0 operation
counts only. The arms are `rho.signed_frobenius` (matched operations-only
reference), `rho.negation` (generic baseline), and `bsgs.negation` (memory
tradeoff). No wall-time or scaling claim will be made from this one-round panel.

The structural companion runs `ghs_screen` with genus bound four for all five
models and `ghs_transport` on the subgroup generator only if the exact point
and subgroup order are available. A small genus or small target field alone
does not imply a cheaper ECDLP. For a certified prime subgroup of order `r`,
nonzero trace into `E(F_(2^l))` requires `r <= #E(F_(2^l))`; otherwise the
trace must kill that subgroup. Where `#E(F_(2^l))` is not counted, report the
condition as a bound or leave it open.

## Hypothesis, boundaries, and accounting

Hypothesis: same-size models can have different subgroup orders and GHS
structural rows, while generic ECDLP operation cost mainly follows subgroup
order and the available Frobenius automorphism. The generic lower-bound
reference is the harness's `sqrt(pi/(2A))` in `S = total group operations /
sqrt(r)`, with automorphism order `A` as recorded per curve. The measured
reference is the matched `rho.signed_frobenius` arm on each target. Report
both per-target correctness and per-curve paired operation counts; do not
pool targets into a single ECDLP solve.

The ECDLP interval includes setup, search, and answer verification as charged
by `ecbench`. BSGS memory is reported separately. The GHS screen is a
construction-stage diagnostic and has no complete-solver total. Neither its
wall time nor its genus is an ECDLP speed ratio. The `sect113r1` and
`F_(2^192)` results already saved under `docs/ghs-transport-evidence/` are
context only; they are not in this matched cost panel.

Success condition: all arms return independently verified scalars on all 40
targets, and the session passes exact replay audit. Any failure or timeout is
kept in the session and reported; no surviving subset is silently substituted.
Stop condition: one frozen session and its audit plus structural screens are
complete, or a concrete implementation/resource blocker is recorded. Inputs,
methods, and unit are not changed after observing results.

## Reproduction and evidence

Build and run the native `ecbench` and `ghs_screen` binaries from the recorded
commit. Preserve the `ecbench` plan, sealed session, audit receipt, structural
JSON, command log, host description, and source revision. Report raw counts and
derived ratios separately. The final report will link an editable diagram and
PDF, and will state which requested large-curve comparisons remain unmeasured.
