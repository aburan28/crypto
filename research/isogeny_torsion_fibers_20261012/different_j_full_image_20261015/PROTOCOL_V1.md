# Frozen full-image factor-base discrimination protocol

**Frozen before implementation or relation measurement.** This follows
[Gate B](../NEXT_GATE_20261015.md) and is stacked on the public-only χ
result PR #1920. The hypothesis is that selecting eight factor-base columns
from the complete rational image of the explicit degree-three map changes
two-sum support or relation recovery on the named order-61 subgroup.
The full image has 1,525 rational points, while the source has 4,575;
the previous exhaustive [fiber audit](../different_j_20261015/reports/run-20261015T1738Z/REPORT.md)
found three preimages per image point. That multiplicity is an input fact,
not an assumed relation multiplier.

## Frozen input and boundary

- Field `F_(67²)` in the pinned polynomial basis; source `a=46,b=37,j=7`;
  destination `a=53,b=25,j=11`; degree-three Vélu map with kernel
  `{O,T,-T}`. The same public `G,Q` of order `r=61`, cofactor `h=75`,
  target, seeds, three rounds and 1,220-trial cap come from the
  [original V1 spec](../different_j_20261015/RUN_SPEC_V1.json).
  Original public input SHA-256
  `fa41958038b0fe01ec1664ca69acf3c44cc67af9feb38a7c8dbf6bf6a3041b65`;
  answer-bearing checker fixture SHA-256
  `3a3504b14887f20ecde8e894a7a57363f3ba9834d45f14d89568023822cbe94b`.
  The new solver fixture must derive only from the original public input.
- The generic negation-rho scale is `sqrt(pi/(2*2)) ≈ 0.8862` group
  operations divided by `sqrt(r)` under its random-walk assumptions.
  The measured reference is same-target, direct source negation rho,
  rerun in the new panel with rho A/A. The prior three-round median
  was 75.709 microseconds of current-thread CPU; it is historical context,
  not a new paired denominator. The full-image method's portable counted
  `S` remains unset until its native field, map, lookup and linear-algebra
  costs are priced in ecbench.

## Deterministic construction

Enumerate every `E_source(F_(67²))` point, including `O`, by increasing
packed `x=(x0,x1)` and each lift by increasing packed `y`. Apply the
frozen Vélu map to each. Deduplicate destination points by packed
destination key, and verify 4,575 source points, 1,525 distinct image
points and exactly three preimages for every image point. Visit the
deduplicated nonidentity image points in increasing destination key.
Project each by `[75]`; discard infinity and duplicate sign orbits under
`min(key(P),key(-P))`. Take the first **eight** distinct projected
orbits. Freeze the eight raw image points, eight projected columns and
eight first-source-preimage witnesses in a new public fixture. Do not
select or tune columns after reading relation outcomes.

## Paired solve, checking and cost

Run the existing source, isogeny-image, destination-direct, rho and rho
A/A arms plus the new full-image destination arm in three rotated rounds
on the exact same public target and seeds. Construct the image and base
inside each cold full-image arm, charging enumeration, all map calls,
deduplication, projection, pair table, failed queries, rank and linear
solve, scalar recovery, source and destination verification, field/group
counters, current-thread CPU, wall and peak RSS. Keep the 1,220-trial cap
and 512 MiB peak-RSS cap. An exhausted arm is censored. Retain stdout,
stderr, exit code, hashes and failed attempts. The independent Rust
checker must verify the base witnesses through its own field/group law,
all recovered scalars against both public targets, the original private
known answer, and one row per `(round,arm)`.

The exact support census counts distinct nonidentity sums of the 16
signed base points for each arm and independently checks every reported
relation hit by group addition. Compare relation counts, trials and rank
before CPU. Compare complete cold current-thread CPU only as shared-host
L0 descriptive data; a stage improvement or local CPU difference cannot
promote the counted-`S` frontier. If support changes, preserve both
selection and cost and repeat on a larger different-`j` order. If support
is equal, record the bounded equality and do not infer universal
equivalence. If setup dominates, separate cold and genuinely reusable
multi-target accounting in a separately frozen follow-on. All outcomes
get a report, quantitative graph, PDF and benchmark mirror.
