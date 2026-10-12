# Parameterized support-envelope construction experiment

This is a standalone, bounded generic Boolean-polynomial experiment. It follows
the exact-support study in `../boolean_product_schedule_20260922/`, whose cache
key cannot hit a changed canonical polynomial system. No production solver
is modified or called. The experiment constructs matrices; it does not solve
systems, collect relations, accept curve inputs, or recover scalars.

## Contract

For generator slot j declare a finite envelope E_j and represent the polynomial
as `p_j(c) = sum(c[j,m] * x^m for m in E_j)` over F2, with `x_i^2 = x_i`.
Envelopes depend only on the base fixture, not future coefficient assignments.
The context is `(n, matrix degree D, active multiplier-variable mask, ordered
generator envelopes)`. Changing the coefficients is allowed. A support escape,
generator-count change, or context change rebuilds directly; it never silently
expands or replaces the envelope.

For each possible multiplier t, the compiler groups source coefficient slots
by the output monomial `m OR t`. Each group stores a bit mask. Its output
coefficient is `popcount(active_coefficients AND group_mask) mod 2`, so colliding
products cancel exactly. This is a linear transformation of the current
coefficient vector, not a stored output matrix.

If a generator currently has degree d, the constructor must visit precisely
the multipliers with degree at most D-d in numeric monomial order. The compiler
therefore covers degrees down to the lowest possible nonzero degree in E_j,
including degree zero when constants are allowed. Within a multiplier of
degree k, only source terms of degree at most D-k need routing: any selected
generator with a higher-degree active term would make that multiplier ineligible.
This proves omitted edges cannot affect an eligible product. An identically
zero generator contributes no rows. Cancellations may remove an entire row.
Only surviving rows determine the compact output columns and the current caps.

Generator slots are ordered. Reordering is accepted only if the polynomial
in each new slot lies in that slot's envelope; the result is recomputed in the
new input order. A coefficient mask has at most 32 bits. Input validation rejects
duplicate or unsorted monomials; canonicalization is never implicit.

Compilation caps and current matrix caps are different. Compilation refuses
when its plan/group/retained-storage limit is exceeded. Current output row and
column caps apply after parity, including on cache hits and direct fallbacks.
The row cap here is 512 (the preceding study used 256): an n=12 constant
generator needs 299 multipliers before counting the other generators. This is
declared before running and applies equally to all arms.

## Frozen design and falsification

`protocol.json` declares all seeds, families, limits, controls and acceptance
criteria before timing. The four arms from the preceding study remain in the
same binary; the envelope is the fifth arm. Cold batches charge one compilation,
all application/fallback work, fresh outputs, validation, and destruction.
Independent reference generation and common fixture/envelope generation are
outside arm timing, and are included in worker process receipts. Retained
capacity bytes and output payload are separate from whole-worker peak RSS.

The fixed grid has 256 cells, five arms, ten rotated repetitions, and batches
1/4/16/64 at n=6/8/10/12. Holdouts use new seeds. Each non-repeat batch contains
distinct canonical systems. The degree-cycle family exercises newly required
multipliers, constants and zero generators. The escape family must fall back
on every fourth instance. No failures or censored cells may be discarded.

Correct reuse requires every expected in-envelope hit and exact oracle equality.
Performance promotion additionally requires **all 48** declared holdout
comparisons at batches 16 and 64 to have a 95% paired-bootstrap lower bound
greater than 1.05 against direct construction, layout reuse, and packed-matrix
caching. Failure is a valid outcome. Intervals describe these fixed seeds and
repetitions only. They do not establish a population or asymptotic result.

The reference boundary is the strongest measured correct generic constructor
under the same workload. Every returned packed matrix still requires materializing
its output words. Output bytes provide a finite workload diagnostic, not a
calibrated time lower bound. No operation conversion, full pipeline or rho ratio
is available; those quantities stay null. This is an engineering/construction
diagnostic under the repository reporting convention. The WDSat/full IC suite
is inapplicable because this worker has no solver or cryptanalytic pipeline.

## Run and replay

Run from the repository root, using a new output directory each time:

```sh
python3 research/boolean_support_envelope_20260922/run.py --out research/boolean_support_envelope_20260922/run_01
python3 research/boolean_support_envelope_20260922/analyze.py research/boolean_support_envelope_20260922/run_01
python3 -m unittest discover -s research/boolean_support_envelope_20260922 -p 'test_*.py'
```

Only Python's standard library and Rust's standard library are used. The runner
freezes source and protocol before compilation, records executable hashes,
retains raw worker output and process receipts, and writes an additive manifest.
Replay in a temporary copy when checking immutable evidence. The independent
set-parity oracle is separate from the sorted-product control and the compiled
parity-mask implementation.
