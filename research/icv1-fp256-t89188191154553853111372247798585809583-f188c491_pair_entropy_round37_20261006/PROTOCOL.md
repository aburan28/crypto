# P-256 cutoff-conditioned pair entropy, round 37: protocol

Date registered: 2026-10-06

## Question

Round 36 rejects uniformly random Montgomery-affine table access at the
registered 256-MiB decision depth.  Does round 33's cutoff-219 restart
condition concentrate actual signed pair accesses enough that a small hot
table can cover almost all eight setup lookups, with cold pairs constructed
on demand?

## Frozen dependencies and scope

Import round 33 by exact SHA-256
`9931bccbd9e2821f4498f465ce65f83efb400e7fb3740a239bdb40d3275ffbd8`
and round 36 by exact SHA-256
`f145231c098da83878e890d340f8ae921598519b0f2626bcc7c4cd8f9446d015`.
Keep curve
`icv1-fp256-t89188191154553853111372247798585809583-f188c491`, factor-base
width 131,458, all 17 variable columns, cutoff 219, eight signed pair
lookups, corrected base ratio `0.998569034150286`, local-oracle ratio
`0.964336477130181`, mean capacity `224.361249307511`, and remaining budget
`0.334410508423` additions per segment unchanged.

The residual capacity histogram is frozen to 6,935 columns in each class
0 through 17 and 6,628 columns in class 18.  No column, sign or rejected start
may be omitted from the conditional distribution.

## Exact signed-pair density

Let `T` be the exact number of distinct 17-column tuples whose capacity sum
is at least 219.  For every unordered capacity-class pair `(d,e)`, compute
exactly the number `C[d,e]` of retained tuples containing one fixed pair of
distinct columns from those classes.  Remove the fixed columns from the
histogram and obtain `C[d,e]` from the coefficient sum for 15 remaining
columns at capacity at least `219-d-e`.

For any one of the four sign assignments of that fixed column pair, its
conditional presence probability is

```text
C[d,e] / (4*T).
```

The signed entry population of a class block is four times its unordered
column-pair population.  Rank all 190 class blocks by exact per-entry density.
For a table of `K` signed entries, allow a fractional final block and fill the
`K` highest-density entries.  This is an optimistic oracle cache: it may pick
arbitrary individual entries from a class block.

For a retained signed 17-column start, let `E_K` be the number of cached
signed unordered pairs present among all 136 possible pairs.  Regardless of
how the selector pairs 16 columns after seeing the whole start, its number of
cache hits is at most `min(8,E_K)`.  Therefore

```text
expected hot hits <= min(8, E[E_K]),
cold constructions >= max(0, 8-E[E_K]).
```

This deliberately ignores matching conflicts and grants adaptive global
pairing, so it is a candidate-favouring upper bound, not a proposed selector.
Verify that the full signed table has expected co-presence exactly 136.

## Registered table depths and cost lower bound

Evaluate `K` equal to the 64-byte Montgomery-affine capacities of 2 MiB,
64 MiB and 256 MiB, plus the smallest optimistic `K` whose cumulative
co-presence reaches `8-0.334410508423`.  Also record the complete
34,562,148,612-entry table.

Grant every hot lookup zero time and zero memory cost.  Grant every cold pair
exactly one extra group addition to construct its sum, ignoring the two
factor-base reads and allocation.  Convert the resulting cold-construction
lower bound with

```text
projected/rho >= 0.998569034150286
               + 0.964336477130181*cold/(224.361249307511+1).
```

This is stricter in the candidate's favour than round 36's measured access
cost.  A cache depth fails if its optimistic projected lower bound exceeds
rho.  The threshold depth is storage accounting only; it is not claimed to
fit any particular cache or achieve the granted zero-cost access.

## Exact controls

- Recompute `T` independently from the frozen histogram and require it to
  match round 33 exactly.
- On at least two deterministic tractable histograms, exhaust every retained
  tuple and compare every fixed-pair completion count with the dynamic
  program.
- Require zero count, ordering, population, normalization, monotonicity,
  false-positive and false-negative failures.
- Hash the ordered class rows and frontier.  Emit deterministic JSON with the
  exact integer counts represented as decimal strings.

## Gates and stop condition

Advance a hot-tier implementation only if an admissible measured table depth
has an optimistic cold-construction lower bound at or below `0.334410508423`
additions per segment.  Passing this screen would not establish measured
time, achieved P-256 coverage, structured degree at most 5, collection below
`2^120`, per-row cost below `2^103`, or non-genericity.

If all registered depths fail, publish the exact minimum entry/byte threshold
and reject capacity-only hot tiering.  Do not benchmark a hot/cold selector or
attempt a full-depth unplanted relation after a failed screen.  Never call the
oracle cache an exhaustive search or count a discarded branch as searched.
