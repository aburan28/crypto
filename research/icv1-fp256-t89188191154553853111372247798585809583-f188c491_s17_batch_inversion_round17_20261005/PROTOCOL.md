# P-256 batched quadratic inversion protocol, round 17

Date frozen: 2026-10-05

Round 16 compressed a fixed sixteen-leaf atom to a 36-byte exact packet but
left its arithmetic unchanged.  One reconstructed target query still charges
192,640 field multiplications, of which 107,520 are denominator inversions in
the image build and query.  This round attacks that measured work rather than
the serialized state.

## Hypothesis and promotion gate

The independent denominators within each balanced image join and within one
128-quadratic target query can use Montgomery batch inversion.  One inversion
of the denominator product plus prefix and reverse products must preserve every
root while reducing counted multiplication work.

For each of round 16's four frozen P-256 packets and both frozen targets:

1. decoded columns, intermediate images, target classifications, image
   digests, and hit digests exactly match rounds 15 and 16;
2. quadratic solves, roots, and lookups remain unchanged;
3. local maximum algebraic degree remains 2, with zero linear or universal
   degeneracy;
4. warm-query field multiplications fall from 88,064 to at most 44,032; and
5. cold one-target field multiplications fall from 192,640 to at most 96,320.

Both deterministic operation gates require at least a 2x reduction.  Any
classification, digest, intermediate-image, or packet mismatch falsifies the
round.  A failed count gate rejects the engineering candidate even if it is
correct.

This is a relation-oracle arithmetic result, not by itself an end-to-end index
calculus speedup.  Relation yield, number of branch cells, matrix construction,
linear algebra, and rho remain unchanged and must not be assigned a speedup.

## Frozen dependency, curve, and factor base

- round-16 result:
  `research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_s17_packed_atom_round16_20261005/packed-result.json`;
- required SHA-256:
  `af2874ef5f0bdabec75b9bf707067232d0d6e3985f0cce4eba6834c2b5790e44`;
- curve:
  `icv1-fp256-t89188191154553853111372247798585809583-f188c491`;
- field: the exact NIST P-256 prime;
- FB1: `FB1h2f8621cda105`;
- FB1 SHA-256:
  `2f8621cda105a51f8703dabe38dec87f3715ab855770acd7b2af0a14906f9e42`;
- columns / signed points: 131,458 / 262,916;
- points SHA-256:
  `70ab21cbda76205ff43ab6bf85c68607be9f6b89df6399e6fa98e3489b888cf1`.

Hash-check round 16, rebuild and verify the complete factor base, decode and
re-encode each packet, and reproduce its selected columns before running the
candidate.

## Candidate algorithm and accounting

Coefficient and discriminant construction and P-256 square-root exponentiation
remain unchanged.  For every nondegenerate quadratic with a square
discriminant, retain its denominator `2a`.  Within one independent join:

1. build prefix products of all denominators;
2. invert the complete product once with the same counted exponentiation used
   by the reference;
3. walk backward to recover every individual inverse; and
4. construct both roots and insert them into the same exact image or lookup
   path as round 16.

The batch implementation must fall back to the scalar solver for a batch of
one, linear equations, universal equations, or any unsupported degeneracy.
Every multiplication in prefix construction, the one exponentiation, reverse
recovery, and root construction is charged.

The frozen composition groups per four samples are:

- eight two-leaf nodes of batch size 1 across both halves;
- four four-leaf nodes of batch size 4;
- two eight-leaf nodes of batch size 64; and
- one target-query batch of size 128 per target.

For a nondegenerate batch of `n > 1`, inversion work is expected to be one
384-multiplication exponentiation plus `3n - 1` product multiplications.  The
frozen expected counts are therefore:

| phase | scalar inversion muls | batched inversion muls | other muls | batched total |
|:--|---:|---:|---:|---:|
| build | 58,368 | 5,802 | 46,208 | 52,010 |
| warm query | 49,152 | 767 | 38,912 | 39,679 |
| cold one target | 107,520 | 6,569 | 85,120 | 91,689 |
| two targets, one build | 156,672 | 7,336 | 124,032 | 131,368 |

Counts, not wall time, are the primary unit.  Do not infer CPU speed from the
count ratio without a separately frozen calibrated timing experiment.

## Exact reference and boundary table

Rebuild every two-, four-, and eight-leaf image with both scalar and batched
solvers and require exact identity flags and ordered x-coordinate sets.  Fresh
exhaustive signed P-256 group addition remains the independent reference for
all intermediate images and all eight target classifications.

| sample | target | scalar / batch positive | scalar / batch muls | ratio | solves / roots / lookups | image and hit digests | FN / FP | degree | class |
|:--|:--|:--|:--|---:|:--|:--|:--|---:|:--|

Report phase-separated multiplication categories, batch sizes, scalar fallback
counts, and exact packet and result hashes.

## End-to-end interpretation

For an end-to-end IC run whose measured fraction `f` is spent in this exact
oracle, the deterministic operation-count ceiling implied by an oracle ratio
`r` is `1 / ((1-f) + f/r)`.  Leave `f` and the resulting total ratio unset
because no full P-256 collector exists yet.  The next admissible experiment is
a charged outer relation scan that measures `f`; this round only establishes
the candidate oracle cost needed by that experiment.

## Scope and stop

Stop after four packets, eight target comparisons, exact intermediate checks,
canonical compact JSON, an immediate byte-identical replay, and standard
repository validation.

Do not claim a changed 96.400018% modeled relation probability, fewer branch
cells, a full P-256 relation, lower linear-algebra cost, a logarithm, an
end-to-end `S`, or a rho speedup.
