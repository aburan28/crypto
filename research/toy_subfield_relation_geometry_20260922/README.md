# Subfield-based factor bases: bounded relation-geometry experiment

Date: 2026-09-22. Classification: **structural diagnostic**.

The experiment tests an algorithm-design question: can a subfield-based
coordinate space make the quadratic part of a binary summation-polynomial
system narrower without destroying relation coverage? It supplies a measured
screen for subsequent toy solver work. It does not supply a discrete-log
algorithm, a GHS implementation, or an Oakley key-recovery capability.

## Why study a different field family?

[RFC 2409, sections 6.3–6.4](https://www.rfc-editor.org/rfc/rfc2409.html#section-6.3)
defines Oakley Groups 3 and 4 over fields of degrees 155 and 185. These factor
as `5*31` and `5*37`, respectively, giving proper subfields and degree-five
tower descriptions. ECC2K-130's degree 131 is prime. That structural difference
motivates a separate research family; it does not itself establish a useful
descent for either Oakley curve.

[Galbraith, Hess and Smart, EUROCRYPT 2002](https://www.iacr.org/cryptodb/archive/2002/EUROCRYPT/2290/2290.pdf)
studies Weil descent over composite-degree binary fields and explains why the
descent's usefulness depends on the particular curve. This experiment measures
factor-base geometry directly on tiny elliptic curves; it does not construct
the higher-genus curve or transport a discrete logarithm.

The existing repository's `legacy_curve_attacks.rs` contains an Oakley Group 3
demo with a planted bounded scalar and BSGS. That is a separate workload, not an
index-calculus baseline. The existing `ghs_descent.rs` also describes limitations
of its higher-genus smooth-model construction. Neither is executed here.

## Frozen question and boundaries

[contract.json](contract.json) was written before the first run. The fields are
fixed at `n = 6, 9, 10, 12`, with proper-subfield degrees `k = 2, 3, 2, 4`.
The executable accepts only these four tiny fields. For each curve and each
coordinate dimension `l = k, 2k`, compare:

- `generic_linear`: one seeded generic binary vector space.
- `subfield_linear`: a one- or two-dimensional space over the proper subfield.
- `scaled_subfield_linear`: the previous space multiplied by a fixed field element.

Each factor base contains **all** rational points whose x coordinate belongs
to its space. Equal coordinate dimensions do not imply equal point counts;
both counts are reported. Every unordered pair of **distinct** factor-base
points is enumerated. Infinity sums are excluded, repeated targets are
deduplicated, and the coverage denominator is the complete curve order.
This measures which targets have a two-point witness, not independence of
rows in a relation matrix.

Before measurement, the counting ceiling is
`distinct nonidentity targets <= min(B*(B-1)/2, #E-1)`.
The product-span bound is `dim span(V*V) <= min(n, l*(l+1)/2)`.
These are bounds on the stated geometry, not security lower bounds.
No DLP is solved, so total calibrated operations, `S`, and ratios to rho or a
DLP floor are null. The frozen WDSat suite is inapplicable to this solver-free
screen; no solver-performance comparison or equivalent-suite claim is made.

The preregistered screening criterion is a lower product-span dimension and
at least the reference's number of covered targets in three of four
full-degree holdout cells at `l = k`. Empty bases and regressions remain in
the output. Passing is a reason to test a solver, not proof that one is faster.

## Why product-span dimension is relevant

For a proper subfield `U` of binary dimension `k`, a scaled subfield line
`V = gamma U` has `span(V*V) = gamma^2 U`, also of dimension `k`.
A two-dimensional subfield space has product-span dimension at most `3k`.
These are exact structural facts; a generic binary space need not have that
closure. For a fixed target abscissa `t`, the quadratic coefficient vectors
of binary S3 are images of these products under the binary-linear map
`z -> z^2 + t*z`. Thus product-span dimension is an upper bound on their rank.

It is **not** the number of independent equations in the full system. Linear
parts and constant terms still matter, and a lower quadratic rank need not
reduce Gröbner solving degree or runtime. The current experiment measures
product-span rank, not target-dependent coefficient rank or solving degree.

## Results and follow-up control

The initial 72 variants passed the geometric screen in **4/4** holdout cells.
Inspection then found that three of those holdout groups had power-of-two
order. Their results are retained verbatim in [results.json](results.json).
For example, at `n=12, l=4` all three bases contain 15 points and cover 98
nonidentity targets, while product-span rank is 10 for the generic base and 4
for either structured base. That is an exact toy geometry comparison, but
the group-order caveat prevents treating it as evidence about a large odd
prime-order subgroup.

A second contract, [contract_followup.json](contract_followup.json), was
written before a separate run. It selects the first seeded full-degree curve
whose largest odd prime divisor has cofactor at most eight. Selection never
looks at factor-base size, rank or coverage. Rejected candidates are retained.
Its 24 variants also passed the same screen in **4/4** cells.

The secondary holdout's selected groups have orders `52=4*13`, `548=4*137`,
`1048=8*131`, and `4184=8*523`. Coverage is still measured over the whole curve;
it is not a subgroup-DLP benchmark. The full matched table is in
[RESULTS.md](RESULTS.md), with raw evidence in
[results_followup.json](results_followup.json).

In the degree-five toy extension `GF(2^10)/GF(2^2)`, at `l=4`, the generic base
has product rank 10 and reaches 86 targets with 15 base points. The structured
bases have product rank 6 and reach 187 and 204 targets, with 23 base points
each. The point-count difference is part of the observation; this is not a
same-size factor-base performance comparison.

Regressions matter: at `n=9, l=3` in the secondary holdout the generic base
reaches 18 targets, the unscaled structured base 50, and the scaled structured
base only 2. There is no universally best family. At relative degree three,
the `l=2k` product span typically fills the ambient field, losing the rank
advantage. The successful criterion permits either structured family per cell;
it does not identify one selection rule that wins on every cell.

Across both runs: **96 variants, 431,576 verified nonidentity pair witnesses**.
Every witness passes the curve equation, an inverse group identity, and the
independent S3 identity. This count includes repeated targets and is not a
count of independent index-calculus relations.

## Interpretation

The new evidence supports testing structured spaces in a bounded toy solver
comparison; it does not establish an algorithmic speedup. The next hypothesis
is whether the reduced quadratic coefficient rank helps after charging domain
encoding, failed targets, lifting and verification. All comparisons must use
identical target workloads. A higher-genus GHS route is a separate experiment.

One generic subspace and only two structured choices per cell limit
generality. Some reference bases have no nonidentity pair sums. Target coverage
is exact, but the sampled curves and bases are not an exhaustive population.
No scaling exponent, rho crossover, Oakley runtime, or cryptographic-size
security claim follows from these measurements.

## Reproduce

From the repository root, using Python 3.10 or newer with no third-party packages:

```bash
python -m unittest discover -s research/toy_subfield_relation_geometry_20260922 -p 'test_*.py' -v
python research/toy_subfield_relation_geometry_20260922/run.py --output /tmp/toy-geometry-new.json
python research/toy_subfield_relation_geometry_20260922/run_followup.py --output /tmp/toy-geometry-followup-new.json
```

Outputs refuse overwrite. The JSON records contain exact bases, curve and
field parameters, multiplicity digests, source/contract hashes, and null DLP
costs. Wall timings include verification and are diagnostic only; no timing
improvement is claimed. The seven unit tests include independent exhaustive
point enumeration on the smallest field, all nonzero-element inversion checks,
group-law checks, subfield closure, target deduplication, and size restrictions.
