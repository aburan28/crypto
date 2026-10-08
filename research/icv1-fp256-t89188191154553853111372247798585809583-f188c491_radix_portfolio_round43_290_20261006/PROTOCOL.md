# P-256 exact radix portfolio, rounds 43–290: protocol

Date frozen: 2026-10-06

This protocol defines 248 distinct selector-routing hypotheses in one frozen
Cartesian sweep. It does not claim 248 adaptive discoveries: each numbered
round is one preregistered cell, and the complete manifest must preserve
negative and projection-only outcomes.

Round 41 compared only 12- and 18-bit radix digits. Round 42 showed that
fixed-universe direct addressing passes the storage gate but loses badly to a
second complete request-generation pass. The remaining exact one-generation
table family is stable radix routing with other digit widths. This sweep
closes that finite width/depth design space.

## Frozen identity and dependencies

- curve:
  `icv1-fp256-t89188191154553853111372247798585809583-f188c491`;
- measured family: round-33 known-log low-delta/scalar selector;
- comparison factor base: `FB1h2f8621cda105`, 131,458 columns;
- pair-key width: 36 bits;
- pair universe: 34,562,148,612 signed entries;
- requests per start: 8;
- request record: 16 bytes;
- streaming accumulator: 16 bytes;
- round-42 result SHA-256:
  `82adf1c47c60a9004a83ac061d9a75ce908f3bb9c28f537260b43b7d5016aa2a`;
- round-42 isolation SHA-256:
  `7a6d11512161cee3020f599e8e705f6ea8fa68ba4e14aada73e0716230d4dad8`.

The binary must verify the exact round-42 schema, curve, dependency hashes,
request digest, native-control digest, zero-failure gates, full-universe byte
total, and non-promotion decision before producing evidence.

`FB1h2f8621cda105` remains a comparison boundary. Its degree-4 residual
does not transfer to the measured known-log family. Relation collection,
per-row cost, and sparse linear algebra remain unset unless a separately
proved composition is supplied.

## Round map

Let the radix digit width be `w in {6, 7, ..., 36}` and let the start-depth
index select:

    depth_index  0   1   2   3   4   5   6   7
    start bits  12  14  16  18  20  24  32  40

The round number is:

    round = 43 + 31 * depth_index + (w - 6)

Thus rounds 43–73 use `2^12` starts, rounds 74–104 use `2^14`,
rounds 105–135 use `2^16`, rounds 136–166 use `2^18`, rounds
167–197 use `2^20`, rounds 198–228 use `2^24`, rounds 229–259 use
`2^32`, and rounds 260–290 use `2^40`.

Every round must receive a separate deterministic JSON artifact. A manifest
must bind all paths and SHA-256 hashes and reject any missing or duplicate
round number or parameter pair.

## Exact router

Use the same deterministic SHA-256 rejection-sampled request stream as rounds
40–42. For width `w`, perform `p = ceil(36 / w)` stable
least-significant-digit radix passes. The final pass uses only the remaining
key bits, so its histogram has `2^(36 - w(p-1))` buckets rather than
`2^w`.

Every executed cell must:

- route all eight requests from every start;
- compare the complete sorted request sequence with a canonical
  `(pair, owner_slot)` reference;
- compare every reconstructed accumulator;
- compare distinct-pair counts and deterministic digests;
- regenerate every reference accumulator independently;
- apply one-unit negative controls;
- report false positives, false negatives, operation counts, traffic, peak
  algorithmic bytes, and wall time.

No probabilistic branch may be discarded or credited as exhaustive.

## Execution boundary

A cell is `complete_checked` only when both:

- start depth is at most `2^18`; and
- the largest histogram has at most `2^24` counters.

All 76 such cells (`w = 6..24` at four depths) must execute. Other cells
are `projection_only`: compute exact pass counts, counter sizes, record
traffic, and peak-memory formulas without allocating or routing them.

A projection-only artifact must set:

- `complete_checked=false`;
- false-positive and false-negative fields to null;
- wall time and addition equivalents to null;
- `applicable_memory_tier_measured=false`.

It must never inherit exactness from another depth or width.

## Registered accounting

For `M=8S` requests and pass histograms `B_j`:

    record traffic = M * (16 generation + 32 per radix pass
                          + 16 sorted scan + 32 accumulator read/write)
    histogram traffic = sum_j 24 * B_j
    peak bytes = 2 * 16M + 16S + 8 * max_j(B_j)
    field/group operations = 0 routing operations

Also record comparisons, histogram increments, scatter writes, sorted scans,
and accumulator updates. Synthetic pair-value work is shared across widths
and is not a native group-operation claim.

## Timing and fitted exponents

Build the release binary before execution and run it through
`tools/isolated_bench.py run`, one thread pinned to one CPU. The binary
must:

- warm representative widths;
- perform an A/A native P-256 addition control;
- retain the one-shot wall time for every complete checked cell;
- select the five fastest widths at `2^16` by one-shot time, breaking ties
  by smaller peak bytes and then smaller width;
- execute seven balanced rotating-order repetitions of those five widths on
  the identical frozen `2^16` request vector;
- report medians, minima, checksums, native-addition equivalents per start,
  and ratios to the 0.334410508423-addition routing allowance;
- fit time exponents for every executed width across start depths
  `2^12, 2^14, 2^16, 2^18`;
- fit exact memory exponents for all widths.

The shortlist timing excludes request generation because the frozen vector is
shared. Import round 42's complete one-generation cost as the complete-route
control; do not add or subtract independently timed medians to manufacture a
new end-to-end measurement.

## Gates

No round is promoted unless all of:

- the exact cell, not another parameter cell, has zero false positives and
  false negatives;
- its structured residual degree is at most 5 for the same measured family;
- projected relation collection is below `2^120`;
- projected cost per usable relation is below `2^103`;
- peak materialised storage is below `2^50` bytes;
- a complete selector is measured below rho in the applicable memory tier;
- no discarded probabilistic branch is counted as exhaustive.

Do not attempt a full-depth unplanted relation unless every gate passes. If
no cell passes, publish the complete 248-round ledger and identify the best
measured routing kernel, the best formula-only width, and the dominant
end-to-end obstruction.
