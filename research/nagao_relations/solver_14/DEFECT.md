# Defect: this round stops at `d = 12`, and the stop was deliberate

`raw.jsonl` holds 14 valid records and then ends.  The run was **killed by its
author**, not by a crash or a budget.  Recorded here rather than left to be
inferred from a short file.

## What is wrong

`plantedTarget()` mints a decomposable target by enumerating every abscissa of
the subspace:

```python
for i in range(1, 1 << d):
    p = curve.pointFromX(f.fromCoords(i))
```

That is `2^d − 1` point lifts, each an inversion and a trace over `F_2^131`.
It is instant at `d ≤ 12`, about 16.7 million lifts at `d = 24`, and does not
terminate in any useful time at the `d = 32, 40, 44, 45, 48, 64` the contract's
M3 asks for.  The contract's own measurement plan therefore could not be
completed by this implementation.

The soundness gate M4 is what needs planted targets, and it is the reason the
enumeration is there at all.  Nothing about the defect affects the records that
*were* written: every `anf` and `span` row in `raw.jsonl` was produced before
the enumeration became slow, and the gate passed on all of them.

## What survives

All 14 records stand, and they already contain both headline measurements:

- **M1/M2**: the descended system has total degree 6 and refutes a
  non-decomposable target at degree 6 with the multiplier set `{1}` — that is,
  **the 131 equations alone, with no Macaulay multipliers**.  A linear
  NO-certificate.
- **M3**: the affine span of the `S₄` value set has dimension `71, 97, 123` at
  `d = 4, 5, 6` and **saturates at 131 from `d = 7` onward**, so the certificate
  exists for `d ≤ 6` and is gone at `d ≥ 7`.
- **M4**: no planted target admitted a certificate at any measured `d`.

The contract's falsifier asks whether the span stays proper at `d ≥ 45`.  These
records do not reach `d = 45`, so **this round does not answer it** — but they
show saturation from `d = 7`, and the span is monotone in `d` (enlarging `V`
enlarges the value set, hence its affine span), so saturation cannot reverse.
That argument is stated here as an argument, not as a measurement.

## What replaces it

`solver_15` re-runs M3 with planted targets minted by **sampling** admissible
abscissas instead of enumerating them, which is what the large-`d` rows need,
and carries the measurement to `d = 45` and beyond so the falsifier is answered
by evidence rather than by monotonicity.

## Why the contract did not catch it

The contract specified the `d` values to measure and the soundness gate, but
said nothing about how targets are minted.  A measurement plan that names
parameters an implementation cannot reach is under-specified; a future contract
listing `d` values should state the cost model for constructing an instance at
the largest of them.
