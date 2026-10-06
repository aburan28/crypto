# P-256 two-pass independent router, round 41: protocol

Date frozen: 2026-10-06

Round 40 measured a complete three-pass independent-start router at 3.43258
native P-256 addition equivalents per start, 10.2646 times its 0.33441 budget.
That row includes SHA-256 generation of the experimental uniform request
stream.  In an actual selector, pair identifiers are outputs of the upstream
start construction.  This round separates those costs and tests whether a
two-pass 18-bit radix router can fit the routing headroom on a pre-generated
request stream.

This is a routing-kernel stage experiment.  A pre-generated-input crossover is
not an end-to-end speedup, and a RAM result may not be transferred to the
multi-terabyte projected tier.

## Frozen dependencies and scope

- curve:
  `icv1-fp256-t89188191154553853111372247798585809583-f188c491`;
- measured family: round-33 known-log low-delta/scalar selector;
- pair universe: 34,562,148,612 signed entries;
- requests per independent start: 8;
- round-40 result SHA-256:
  `52154eb455ea7edd187d83944ba4b6d953d92014d4f05919e4bd8b7bc097ecbf`;
- round-40 isolation SHA-256:
  `c04812b438ea26ed6b0f94fee6a7425730fd49b395c9f92c70db0c61915b59a9`.

Import and verify the exact round-40 constants, depth-`2^16` request digest,
control digest, zero-failure gates, pair-only threshold, selected storage row,
and measured complete-router result.  Reject a mismatch.

`FB1h2f8621cda105` remains only a separate degree-4 comparison base.  The
measured family's structured residual degree is unknown.  No round-41 timing
may be used for the 138,031-row Dickson relation projection.

## Registered variants

Use the identical deterministic request stream and 16-byte request records.
All variants perform an exact sorted distinct-pair scan and update one 16-byte
accumulator per start.

1. **complete 3x12 baseline**: generate requests with round 40's SHA-256
   rejection sampler, then perform three stable 12-bit radix passes;
2. **pre-generated 3x12 control**: clone the frozen request vector inside the
   timed interval, perform the same three passes, scan, and reconstruct;
3. **pre-generated 2x18 candidate**: clone the same vector, perform two stable
   18-bit radix passes, scan, and reconstruct.

The pre-generated rows exclude upstream request derivation and must say so.
Their counted traffic still includes the initial 16-byte record write/clone.
Registered logical traffic is 160 bytes/request for 3x12 and 128 bytes/request
for 2x18.

## Correctness and depth ladder

At start depths `2^12`, `2^14`, `2^16`, and `2^18`:

- compare both radix outputs byte for byte after sorting;
- compare distinct counts and every reconstructed accumulator;
- regenerate every expected accumulator independently;
- apply one-unit negative controls;
- record false positives, false negatives, radix counts, traffic, RAM, and
  deterministic digests.

Repeat complete exhaustive toy references over pair universes 17, 31, and 63.
The native group replay remains round 40's frozen 4,096-control result; this
round changes routing only and must verify the imported digest.

## Timing

Build the release binary before timing.  Run it through
`tools/isolated_bench.py run`, one thread pinned to one CPU.  The binary must:

- warm every variant;
- perform an A/A native-addition control;
- execute at least seven balanced rotating-order repetitions at `2^16`
  starts;
- record raw wall times and identical checksums;
- report medians, minimums, 3x12/2x18 ratio, native-addition equivalents per
  start, and ratio to the 0.334410508423-addition routing budget.

Run a second timing ladder at `2^14`, `2^16`, and `2^18` with at least five
repetitions per variant.  Fit the observed time exponent against starts.  Wall
time is secondary to exact traffic and correctness; contended runs are not
pooled.

## Projection and promotion gates

Combine the measured router equivalents only with round 40's pair-construction
rows as an explicitly optimistic RAM-stage projection.  Report the resulting
ratio at `2^38`, `2^40`, and the largest batch below `2^50` bytes.  For every
row record whether its memory tier was measured; projected external-memory rows
must remain unverified.

Promotion requires all of:

- zero false positives and false negatives;
- exact equality with the three-pass reference on complete checked instances;
- structured residual degree at most 5 for the measured family;
- relation collection below `2^120` and per usable relation below `2^103`;
- peak materialised storage below `2^50` bytes;
- measured complete selector time below rho in the applicable memory tier;
- no discarded branch counted as exhaustive.

Do not attempt a full-depth relation unless every gate passes.  If 2x18 fits
the RAM routing budget but its projected rows require unmeasured external
memory, publish that split result without claiming parity.

