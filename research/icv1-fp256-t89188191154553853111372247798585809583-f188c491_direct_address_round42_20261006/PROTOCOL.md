# P-256 fixed-universe direct-address replay, round 42: protocol

Date frozen: 2026-10-06

Round 41 reduced exact radix traffic by 20% but retained storage linear in the
number of independent starts.  Its optimistic sub-rho rows require 299 TB to
1.126 PB of materialised records in an unmeasured tier.  This round tests the
next representation boundary: exploit the fixed 34,562,148,612-entry signed
pair universe directly, store one value per pair, and regenerate the start
stream rather than storing every request.

This is a representation-stage experiment.  Scaled dense-table timing is not
timing of the 557-GB P-256 table, and deterministic request replay is not free.
No stage-only result may be called end-to-end parity.

## Frozen dependencies and scope

- curve:
  `icv1-fp256-t89188191154553853111372247798585809583-f188c491`;
- measured family: round-33 known-log low-delta/scalar selector;
- pair universe: 34,562,148,612 signed entries;
- requests per independent start: 8;
- value width: 16 bytes;
- seen bitmap: one bit per pair;
- round-41 result SHA-256:
  `8bf3dd1d2e7942858c6cbf8486170cb6275aa0ec5e5380e686daeebcd9fd26a6`;
- round-41 isolation SHA-256:
  `6e6dedb3f9367229808b1d3eb7c9978c9445fd01643a72078a8489215d9ebf32`.

Import and verify the round-41 schema, curve, request digest, exactness gates,
native-control digest, and non-promotion decision.  Reject any mismatch.

`FB1h2f8621cda105` remains only a separate degree-4 comparison base.  This
experiment does not transport its structured degree to the measured family
and does not project a 138,031-row Dickson relation collection.

## Candidate and matched reference

For a registered power-of-two universe `U`, reduce the same SHA-256
rejection-sampled request stream modulo `U`.

1. **materialised 3x12 reference**: generate every 16-byte request record
   once, perform three stable 12-bit radix passes, construct each distinct
   pair value once, and reconstruct every start accumulator;
2. **direct-address replay candidate**:
   - first pass: regenerate every request, test a one-bit seen map, and on the
     first occurrence construct and write the 16-byte pair value;
   - second pass: regenerate every request, read its direct-addressed value,
     and reconstruct every start accumulator.

The candidate must not retain request records or discard requests.  Its
materialised state is exactly `ceil(U/8) + 16U + 16S` bytes for universe
`U` and `S` starts, plus fixed vector metadata.  Report this formula
separately from allocator RSS.

Use a deterministic 128-bit pair-value function independent of request order.
Compare all accumulators, distinct counts, and deterministic digests with the
materialised reference.

## Correctness and scaling

Run complete exhaustive toy controls over universes 17, 31, and 63 and start
counts 4, 8, and 16.

Run a matched exact ladder at load one (eight requests per start and
`8S = U`) for `U = 2^12, 2^14, 2^16, 2^18, 2^20, 2^22`.  At every cell:

- compare every reconstructed accumulator;
- compare exact distinct-pair counts and expected first-occurrence counts;
- regenerate selected values independently;
- apply one-unit negative controls;
- record false positives, false negatives, candidate/reference digests,
  pair-value constructions, peak materialised bytes, and logical traffic.

Run a saturation ladder at `U = 2^18` with request loads 1/4, 1, 4, and 16.
This measures whether state becomes universe-bounded while every request is
still replayed exactly.

## Timing

Build the release binary before timing.  Run it through
`tools/isolated_bench.py run`, one thread pinned to one CPU.  The binary
must:

- warm both variants;
- perform an A/A native P-256 addition control;
- execute at least seven balanced rotating-order repetitions at
  `U = 2^20`, `S = 2^17`;
- time complete request generation, allocation, pair-value construction, and
  reconstruction for both variants;
- record raw wall times, medians, minima, identical checksums, native-addition
  equivalents per start, and ratio to the frozen 0.334410508423-addition
  routing allowance;
- run a five-repetition load-one timing ladder at
  `U = 2^16, 2^18, 2^20, 2^22` and fit time and memory exponents.

No timing from a table smaller than the full pair universe is an applicable
external-memory measurement.

## Traffic and projection

Record two traffic models:

1. semantic bytes: bitmap bits, first-occurrence value writes, second-pass
   16-byte value reads, and 16-byte accumulator read/write;
2. a 64-byte cache-line model that charges a full random line for each bitmap
   probe and direct value lookup, plus write allocation/writeback for new
   values.

Project the exact full-universe state and both traffic models at `2^38`,
`2^40`, and the largest start count below the `2^50` storage gate.  Combine
only measured arithmetic equivalents with round 40's pair-construction model,
label the result optimistic, and keep `applicable_memory_tier_measured=false`.

## Promotion gates

Promotion requires all of:

- zero false positives and false negatives;
- exact equality with the materialised reference on complete checked cells;
- structured residual degree at most 5 for the measured family;
- relation collection below `2^120` and cost per usable relation below
  `2^103`;
- peak materialised storage below `2^50` bytes;
- measured complete selector time below rho in the applicable memory tier;
- no discarded branch counted as exhaustive.

Do not attempt a full-depth relation unless every gate passes.  If the fixed
table passes the storage formula but only scaled RAM timing exists, publish
that split result without claiming parity.
