# Transporting quadratic coefficient row spaces

This standalone experiment tests a changed inference policy in complete generated
Boolean solves. All nine retained methods, including the current fastest
`packed_state`, remain controls. There are no curve inputs or production solver
changes. Complete index-calculus, calibrated-operation and rho costs remain null.

## Mathematical contract

Let F be a list of squarefree degree-at-most-two polynomials over GF(2), and let
S_a be substitution of a Boolean partial assignment a. Substitution is linear
on coefficient vectors, so

\[
S_a(\operatorname{span}_{\mathbb F_2} F)
=\operatorname{span}_{\mathbb F_2} S_a(F).
\]

Consequently a retained row-space basis can be specialized and reduced again
without returning to the original generator list. Row reduction preserves the
zero set, and specialization cannot increase row rank. This is coefficient
linear algebra: no polynomial-multiplier rows or ideal closure are claimed.

The global coordinate order is the original canonical order: quadratic monomials
in ascending mask order, then linear variables and the constant. At n<=36 the
domain has C(n,2)+n+1<=667 columns, stored in at most eleven u64 words per row.
Input rows number at most 36; row combinations may become dense without leaving
the domain. A variable assignment removes its incident quadratic and linear
columns, toggles surviving quadratic coefficients into linear columns when the
value is one, and toggles the current linear coefficient into the constant.
Multiple assignments compose in increasing variable order, with exact parity
cancellation and degree drops.

## Two policies, each with its own exact reference

`basis_list` and `basis_wide` maintain canonical full RREF. They extract affine
consequences, apply unit assignments and choose the most frequent variable in
the reduced basis, with low-index ties and the zero branch first. The list arm
uses explicit sorted polynomial symmetric differences. The wide arm uses packed
coefficient rows. Both must match models, every logical counter and trace digest.

`tail_list` and `tail_wide` retain quadratic rows in echelon form and canonically
reduce only the affine tail. That tail contains all affine vectors in the row
space. They additionally carry original-equation residuals in the previous packed
representation and use those residuals for the previous frequency heuristic.
Both representations are updated after every assignment, and their costs are
charged. Their stored bases are deterministic but the quadratic part is not
claimed to be canonical full RREF. The list/wide pair must again match exactly.

All four policies perform complete Boolean search when their limits permit it.
Different inference policies need not have the same tree or model as the older
controls. SAT models are checked on original equations. UNSAT requires a completed
retained search reference. UNKNOWN is censored, never UNSAT. Unsupported source
domains use the retained fallback; all benchmark fixtures fit every new domain.

## Discovery probes and frozen acceptance conditions

The two `resource_probe_*` bundles use only n24, discovery seed 17, all three
families and one repetition. Full RREF changed the branching heuristic and
increased nodes on two families. Preserving source-equation branching reduced
nodes in the second probe, but both wide policies still lost to `packed_state`.
These probes size resources and guide the mechanism. They are not holdout or
promotion evidence, and their source and results remain immutable.

`protocol.json` freezes n=12/16/20/24, three families, two discovery seeds, all four
prior regression/holdout seeds, two new holdout seeds, and thirteen rotated arm
positions. The full grid has 96 distinct generated systems and 16,224 timed
observations. No parameters change after holdout measurement.

For each wide candidate, a dramatic-gain gate requires a 95% paired-bootstrap
lower bound above 2.0 against the pointwise fastest of **all twelve other arms**
on every family at n16/20/24 in both regression and holdout splits. All eighteen
comparisons and every grid completion/verification must pass. An incremental
gate requires the same lower bounds above 1.0. A separate matched-backend gate
requires a lower bound above 1.05 against the same-policy list implementation;
that result cannot replace the strongest complete-method reference.

## Accounting and checks

Cold time includes coordinate/source compilation, state copies, row reduction,
affine propagation, tracing, every branch, destruction and validation. Fixture
generation, preparation of the independent status reference and formatting are
outside arm timings and inside fresh-process receipts. Reduction counters charge
all input rows, including zeros, and the actual union of coefficient columns.
Specialization counters charge input terms in every carried representation.
New-policy width histograms use nominal unassigned-variable count. Logical counts
are work diagnostics, not calibrated processor-operation totals.

Seventeen Rust tests include 442,368 small row-space/partial-assignment cases,
comparison with rebuilding from original generators, independent small truth
tables, dense repeated specialization across word boundaries through n36,
full-solve model/counter/trace comparisons and resource censoring. Trace hashes
are diagnostics rather than cryptographic certificates. Timing intervals concern
repeated measurements of fixed fixtures, not all Boolean systems. Whole-worker
RSS includes all arms; candidate-specific peak memory is unmeasured.

The WDSat/full-curve protocol is inapplicable to this generic mathematical driver.
No full index-calculus result or exponent change follows from a backend speedup.

## Reproduction

```sh
python3 research/boolean_basis_transport_20260923/run.py \
  --out research/boolean_basis_transport_20260923/run_01
python3 -m unittest discover \
  -s research/boolean_basis_transport_20260923 -p 'test_*.py'
```

Choose a fresh output directory. The runner freezes all source, configuration and
analysis files, compiles and runs the correctness tests, then starts the grid.

The later discovery-only census can be replayed with
`python3 research/boolean_basis_transport_20260923/affine_probe/run.py --out <NEW_DIRECTORY>`
after the full comparison finishes. The source-only copying fix in that runner
prevents a frozen-bundle replay from copying old output manifests. The original
census and corrected-wrapper replay are both retained as `affine_opportunity_01`
and `affine_opportunity_02`; their raw counts and summaries are byte-identical.
The latter records the compiler version explicitly. Neither contains performance
promotion timings. See [AFFINE_SUBSTITUTION_NEXT.md](AFFINE_SUBSTITUTION_NEXT.md)
for the unimplemented next mathematical contract.
