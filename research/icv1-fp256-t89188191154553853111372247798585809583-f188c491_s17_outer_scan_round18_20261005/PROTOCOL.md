# P-256 charged S17 outer-relation scan protocol, round 18

Date frozen: 2026-10-05

Round 17 reduced one fixed sixteen-column atom query to 39,679 charged P-256
field multiplications, but did not price how many atoms or seventeenth columns
a relation collector must visit.  This round measures that missing outer
boundary.  It is a bounded scanner and projection, not a P-256 discrete-log
claim.

## Hypothesis and two gates

For one hash-frozen sixteen-column atom, scanning a nested hash-selected prefix
of signed seventeenth factor-base columns should:

1. preserve the exact round-17 batched oracle;
2. agree on every trial with an independent table of all `2^16` signed sums of
   the atom;
3. find a length-17 planted control and produce no false positive or false
   negative on either the control or the hash-public target;
4. expose approximately one additional bit of work per added outer-prefix
   depth; and
5. retain local oracle degree 2 and the independently measured residual-degree
   envelope 3, 3, 4 for residual depths 1, 2, 3.

Those are correctness and measurement gates.  The separate attack-promotion
gate requires the projected work to collect 105% of the 131,458-column matrix
to be below `2^120` charged field multiplications, with maximum measured
structured residual degree at most 5.  Work at or above `2^128` is an explicit
negative result against the generic-rho exponent.  A low degree does not
override a failed enumeration gate.

## Frozen dependencies and identity

- curve:
  `icv1-fp256-t89188191154553853111372247798585809583-f188c491`;
- factor base: `FB1h2f8621cda105`;
- columns / signed points: 131,458 / 262,916;
- factor-base SHA-256:
  `2f8621cda105a51f8703dabe38dec87f3715ab855770acd7b2af0a14906f9e42`;
- point-set SHA-256:
  `70ab21cbda76205ff43ab6bf85c68607be9f6b89df6399e6fa98e3489b888cf1`;
- round-17 input:
  `research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_s17_batch_inversion_round17_20261005/batch-result.json`;
- required round-17 SHA-256:
  `57cf49ffbab65356b4ef69848e7824f0377bc058ddd57fe0bf6a8f7b77089ce7`;
- round-6 residual-degree input:
  `research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_dickson_residual_scaling_round6_20261004/degree-result.json`;
- required round-6 SHA-256:
  `71d63031111ba48ff831e79430bc87e6c7bac34626d7d0de4acfb4eaa65d40f4`.

Use round 17's `hash-0` packet, whose columns are
`[114616,101224,36382,129678,17773,72570,33360,118255,88179,58271,40458,117426,13544,42490,88617,31941]`
and whose exact 36-byte encoding is
`6fee18b6823879fa8e115b51b7a20941cdef561cce39f27829cab20d3a0a5fa568a47cc5`.
Hash-check both dependencies and rebuild the complete factor base before the
scan.

## Frozen targets and outer stream

The unplanted public target scalar is SHA-256 of

```text
icv1-fp256-t89188191154553853111372247798585809583-f188c491/s17-outer-scan-round18/public-target-0
```

reduced modulo the P-256 subgroup order, with zero replaced by one.  Its
preimage hash is
`77b340b378fd735bc07016816590020379a287c0df90e8e29834b13e1fd767c5`.
The scalar is retained only for independent point verification; the target is
not constructed from factor-base columns.

Visit outer column indices with the full-cycle affine permutation

```text
column(j) = (90322 + 23509*j) mod 131458
```

where the offset and stride come from the named SHA-256 preimages
`outer-offset` and `outer-stride`; `gcd(23509,131458)=1`.  Skip any of the
sixteen packet columns without replacing it.  Record cumulative checkpoints
after `2^d` stream positions for `d = 4,5,6,7,8,9,10`.

The planted target is the sum of the sixteen packet low points and the low
point at stream position 7, column 123,427.  It is a correctness control, never
a relation-yield observation.  Scan both signs of every eligible outer point
for both the planted and public target.

## Candidate and independent reference

Build the two exact eight-leaf atom images once with round 17's Montgomery
batch inversion.  For target `T` and outer low point `P`, query both `T-P` and
`T+P` through the batched quadratic image join.  Charge coefficient,
discriminant, square-root, batch-inversion, and root-construction
multiplications for every call, including misses.

Independently enumerate all `2^16` signed affine sums of the packet points into
a full-point index.  For each outer trial, look up the exact adjusted point in
that index.  This reference does not call the quadratic oracle.  Every emitted
candidate hit must have a reference witness; the planted witness must be
recovered.  Preserve duplicate sums and canonical sign masks well enough to
verify a complete 17-column relation by direct group addition.

Report, at every checkpoint and per target:

- stream positions, skipped packet columns, eligible outer columns, signed
  branch cells, oracle calls, returned roots, lookups, and hits;
- candidate/reference agreement, false negatives, false positives, and direct
  relation-verification counts;
- phase-separated field multiplications and group additions;
- batch-size histogram and scalar fallbacks;
- 36-byte persistent packet state, 8,192-byte image state, logical peak
  candidate scratch, factor-base raw point bytes, and separately excluded
  reference-table bytes; and
- canonical trial-stream and relation SHA-256 digests.

Counts are primary.  Do not headline wall time.

## Projection and sparse linear algebra

For a fixed sixteen-column atom, a complete outer scan has `131458-16`
eligible columns, two signs, and therefore `2*(131458-16)` oracle calls.  Its
exact signed candidate domain is `(131458-16)*2^17`; divide by the P-256
subgroup order and use `1-exp(-lambda_atom)` for the explicitly labelled
uniform-sum Poisson projection.  Derive expected atom scans per relation and
multiply by the measured full-outer per-atom cost.  This is an enumeration
projection; it is not a measured P-256 relation.

Project a collector with

```text
rows = ceil(1.05 * 131458)
row weight = 17
```

and report:

- target count implied by the separately labelled whole-factor-base Poisson
  probability 0.9640001820745057;
- projected relation-collection work and its base-2 logarithm;
- sparse nonzeros and a concrete CSR byte model using four-byte column indices,
  one-byte signs, and eight-byte row offsets;
- `2*columns` sparse matrix-vector products for a conservative Wiedemann work
  model, plus the separately stated Berlekamp--Massey scalar-operation model;
  and
- at least three 32-byte field vectors in the working-memory model.

Do not add unlike operation units.  Compare the collection field-
multiplication exponent with `2^120` and `2^128`; report sparse linear-algebra
additions separately.

## Stop conditions and claims

Abort on either dependency hash mismatch, factor-base identity mismatch,
packet mismatch, off-curve target, duplicate outer-stream index, oracle /
reference disagreement, missing planted witness, incorrect direct relation,
or a residual-degree receipt outside 3, 3, 4.  Stop after the depth-10
checkpoint and the frozen projections.

This round may claim a measured bounded outer-scan slope and a projected naive
atom-enumeration exponent.  It may not claim a full P-256 relation, a measured
complete collector, a measured end-to-end speedup, a logarithm, or a changed
degree of regularity for the unsplit S17 ideal.
