# P-256 pair-table representation audit, round 34: protocol

Date registered: 2026-10-06

## Requested scope and evidence plan

| continuing requirement | round-34 evidence |
|:--|:--|
| pursue parity with rho | test whether the round-33 cutoff-219 stage retains its `0.9985690342` rho estimate after the omitted table-access work is represented honestly |
| retain all 17 variable columns | import the exact round-33 cutoff and its without-replacement 17-column distribution unchanged |
| exact group replay | validate every decoded table entry on P-256 and compare every eight-entry accumulated point with an independently constructed reference |
| respect storage gate | price compressed, affine and projective representations over all `34,562,148,612` signed-pair entries |
| preserve promotion gates | keep measured local timing, projected full-table behavior, relation probability, degree, collection and non-generic status separate |

## Question

Round 33 stored each signed pair sum as a 33-byte compressed P-256 point but
charged only eight table lookups and eight additions per selected segment.  A
compressed point is not directly addable: its `y` coordinate must first be
recovered.  Does any complete representation of the pair table fit below
`2^50` bytes and keep the cutoff-219 stage at or below rho after lookup and
decoding work is expressed in measured P-256-addition equivalents?

This experiment audits only the implementation cost omitted by round 33.  It
does not convert that round's occurrence upper bound into a proved relation
probability.

## Frozen dependency and parity budget

Import the canonical round-33 JSON by exact SHA-256.  Require cutoff 219,
mean capacity `D=224.36124930751055`, eight pair entries and eight setup
additions per segment, one target-correction addition, and the corrected stage
ratio `r0=0.9985690341502859`.

Let `r_local=0.964336477130181` be the frozen local-oracle ratio.  The maximum
additional work that can be admitted without crossing rho is

```text
extra additions / sample = (1-r0)/r_local,
extra additions / segment = (D+1)*(1-r0)/r_local.
```

Record this threshold before timing.  A representation passes the local
timing screen only when the upper confidence endpoint of its incremental
median cost is below that segment budget.  Negative incremental values caused
by timing noise are clamped to zero.

## Registered representations

Evaluate all of the following complete layouts:

1. `compressed33`: SEC1 parity byte plus canonical `x`.  Each lookup must load
   33 bytes, compute `rhs=x^3-3x+b`, recover `y=rhs^((p+1)/4)`, select the
   registered parity, validate `y^2=rhs`, construct a projective point and add.
2. `affine64`: canonical `x || y`.  Each lookup loads 64 bytes, constructs a
   projective point with `z=1`, validates a deterministic sample on-curve and
   adds.  The trusted materialized table need not repeat the on-curve test in
   the timed loop.
3. `projective96`: the implementation's three 32-byte Montgomery coordinates.
   Each lookup loads the directly addable 96-byte projective point and adds.
4. `coefficient32`: because the round-31 base exposes known coefficients,
   load eight 32-byte pair coefficients, include the seventeenth coefficient,
   reduce their signed sum modulo the group order and perform one constant-time
   256-bit scalar multiplication.  Report this generic representation but do
   not treat it as a non-generic index-calculus step.

For every layout report complete projected bytes, `log2(bytes)`, and the
`2^50` storage gate.  No layout may use the 33-byte projection while timing a
larger already-decoded object.

## Deterministic timing design

Use the repository's constant-time P-256 field and complete projective point
formulas in one release binary.  Generate valid table entries deterministically
from hash-selected nonzero scalars and verify reference equality before timing.

Use table working sets of `2^10`, `2^15` and `2^20` entries for every layout,
plus the largest layout-specific power of two that stays below 512 MiB.  Visit
eight independently hash-selected indices per segment.  Precompute the index
stream outside timed regions.  Benchmark these loops separately:

- index-stream/control baseline;
- eight directly supplied projective additions without table lookup;
- layout-specific lookup, decode and eight-point accumulation;
- one constant-time scalar multiplication for `coefficient32`.

Run one untimed warmup and nine timed repetitions in an interleaved,
deterministically rotated order.  Each repetition must run enough segments to
exceed 50 ms for the direct-add baseline where practical.  Record wall time,
process CPU time, checksum, compiler version, target, CPU model, table size,
bytes touched, and all raw repetitions.  Use medians; estimate a deterministic
nonparametric 95% interval for the median by sorted repetition ranks 1 and 7
(zero-indexed) out of nine.  Report both total and incremental costs relative
to the direct-add loop.  The decisive screen uses the fastest observed
working-set median as a candidate-favouring lower estimate and its registered
upper endpoint as the pass criterion.  Larger unmaterialized tables are never
claimed to be faster than the measured minimum.

## Correctness controls

- Check `size_of::<P256ProjectivePoint>() == 96`; otherwise record the actual
  size and use it in storage projections.
- Round-trip at least 4,096 deterministic valid points through compressed and
  affine layouts; reject invalid parity, square-root, curve or equality cases.
- For at least 4,096 deterministic eight-entry segments, require compressed,
  affine, projective and coefficient reconstructions to equal the independent
  reference point.
- Hash the canonical inputs, index stream, decoded outputs and raw timing rows.
- Re-run the release command and require byte-identical semantic output after
  removing explicitly identified volatile host timing fields; timing rows
  themselves are retained, not fabricated as deterministic artifacts.

## Gates and stop condition

The round-33 pair-table stage survives only if at least one faithful layout:

- has zero correctness failures;
- projects below `2^50` materialized bytes;
- has a registered timing upper endpoint within the frozen parity budget; and
- makes no assumption that an unmeasured 1.14--3.32 TB table is faster than the
  fastest measured cache/DRAM working set.

Even if this local screen passes, end-to-end promotion still requires proved
P-256 usable-relation probability, structured residual degree at most 5,
collection below `2^120`, per-row cost below `2^103`, and a demonstrably
non-generic algorithm.  Attempt no full-depth unplanted relation unless all
gates pass.  If every faithful layout fails, publish the measured overhead and
identify representation/decoding as the dominant obstruction for this
round-33 selector, without generalizing to all factor bases or selectors.
