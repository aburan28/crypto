# P-256 independent-start batch routing, round 40: protocol

Date frozen: 2026-10-06

Round 38 rejected correlated codebooks because hot pairs removed independent
start entropy.  This round tests the remaining construction: preserve every
independent start, batch its eight signed-pair requests, deduplicate the
requests globally, construct each distinct pair once, and route the result
back to its start.  The candidate may approach rho only if pair reuse saves
more than exact routing, reconstruction, and storage cost.

This is a selector-stage experiment.  It is not an end-to-end speedup unless
every promotion gate below passes.

## Amendment 1: factor-base scope correction (before execution)

The initial protocol text identified the 131,458-column pair universe as
`FB1h2f8621cda105`.  That is not a valid composition.  The near-rho constants
imported from round 36 descend from round 33's known-log low-delta/scalar
family; round 19's structured residual maximum of 4 belongs only to
`FB1h2f8621cda105`.

This amendment was committed before implementation or execution.  The
experiment therefore has two explicitly separate readings:

1. the radix router is measured for the round-33/36 131,458-column known-log
   selector family, whose structured residual degree remains unknown;
2. `FB1h2f8621cda105` is a comparison boundary only.  Its degree-4 residual is
   preserved, but no round-33 rho constant is transferred to it.

Consequently this round cannot discharge the original 138,031-row Dickson-base
relation-collection gate.  Any sub-rho row is a routing-stage result for the
known-log selector and must leave the Dickson-base end-to-end cost unset.

## Frozen curve, factor base, and dependencies

- curve:
  `icv1-fp256-t89188191154553853111372247798585809583-f188c491`;
- measured selector family: round-33 known-log low-delta/scalar construction,
  modelled with 131,458 columns;
- comparison factor base: `FB1h2f8621cda105`, 131,458 columns, not composed
  with the round-33 cost constants;
- relation model: 17 distinct variable columns, signed balanced `8+9`;
- signed pair universe: 34,562,148,612 entries;
- pair requests per independent start: 8;
- useful-row target: 138,031 independent rows including allowance;
- round-19 selector dependency SHA-256:
  `3096540621408e4a48cfa18963ad01da9686d3527ee26776c8b6cf6f45a71114`;
- round-36 Montgomery-affine dependency SHA-256:
  `f145231c098da83878e890d340f8ae921598519b0f2626bcc7c4cd8f9446d015`;
- round-38 support dependency SHA-256:
  `4297328eee07822331c68bddd65d8b0155ba1c72c6e0dda044dc97d42aba5cce`.

The binary must reject changed hashes, schemas, curve names, factor-base
identity, pair-universe size, base ratio, local-oracle ratio, mean capacity,
or allowed-addition budget.

Frozen comparison constants imported from round 36 are:

```text
base ratio to rho                 0.998569034150286
local oracle ratio to rho         0.964336477130181
mean segment capacity           224.36124930751055
allowed extra additions/start     0.33441050842298764
```

The pair-construction-only projection is

```text
ratio = base_ratio
      + local_ratio * distinct_pairs_per_start / (mean_capacity + 1).
```

Routing cost is additional and may not be silently assigned zero.

## Exact independent-start stream

For start `s` and slot `j`, derive a uniform candidate pair identifier by
SHA-256 domain separation and rejection sampling into `[0,N)`.  The stream is
fixed before execution.  All eight requests remain present; no start or branch
is discarded for deduplication.

Registered depths are `2^12`, `2^14`, `2^16`, `2^18`, and `2^20` starts.
For each depth record:

- request and distinct-pair counts;
- observed reuse and the exact uniform-occupancy expectation
  `N*(1-(1-1/N)^(8S))`;
- request survival ratio and deviation from expectation;
- logical record bytes, peak materialised bytes, and routed accumulator bytes;
- radix passes, bytes touched, route updates, and replay failures;
- candidate false positives and false negatives.

Provide complete exhaustive references for small pair universes and compare
the radix result with a direct ordered-set implementation.

## Structured radix router

Use a fixed 16-byte record containing a 64-bit pair identifier and a packed
start/slot identifier.  Perform exactly three stable 12-bit LSD radix passes,
then one sorted scan.  Construct a deterministic contribution once per
distinct pair and update the owning start accumulator for every request.

Regenerate the original request stream independently and verify every start's
eight-slot completion mask and contribution checksum.  For 4,096 selected
requests, independently replay the routed pair contribution in the native
P-256 group against direct construction from deterministically derived endpoint
points.  These are routing controls; round 34/36 remain the evidence for the
real factor-base pair representation.

Count logical traffic explicitly:

```text
request generation write       16 bytes/request
three radix read/write passes  96 bytes/request
sorted scan read               16 bytes/request
accumulator read/write         32 bytes/request
total                         160 bytes/request
```

Report projected external traffic separately from peak storage.  Byte counts
are not group operations until a measured conversion is supplied.

## Timing protocol

Counted requests, pair constructions, radix passes, traffic, and group
operations are primary.  Wall time is secondary.

Run the release binary through `tools/isolated_bench.py run` with one thread.
The binary performs one warm-up, an A/A timing control, then at least five
interleaved baseline/candidate repetitions.  The baseline is native P-256
group addition on a deterministic point stream; the candidate is the complete
RAM radix route and reconstruction at the registered timing depth.  Record
host, compiler, pinned CPU, raw timings, checksums, context-switch/pressure
metadata from the isolation wrapper, medians, and the A/A spread.

No RAM timing may be extrapolated as disk timing.  If the projected batch
exceeds measured RAM, leave complete time unset unless an external-memory run
is completed under the same isolation rules.

## Full-depth projection and gates

Solve for the smallest batch whose exact uniform-occupancy expectation makes
the pair-construction-only row reach rho.  Also sweep powers of two through the
largest batch whose materialised request buffers plus accumulators fit below
`2^50` bytes.  For each projected row report:

- starts, requests, expected distinct pairs and additions/start;
- pair-construction ratio to rho;
- record bits and bytes, peak storage, radix traffic, and traffic/start;
- required sustained routing bytes per native group-addition time to stay
  inside the frozen rho headroom;
- whether that throughput is measured in the applicable memory tier;
- projected number of batches and traffic for the complete relation search;
- 138,031-row collection, duplicate/rank allowance, and sparse linear algebra.

Promotion requires all of:

- zero false positives and false negatives on complete checked instances;
- exact group replay of every reported native routing control;
- structured residual degree no greater than 5 for the measured round-33
  family; the comparison base's degree-4 result is not transferable;
- projected relation collection below `2^120` field/group-equivalent
  operations and below `2^103` per usable relation;
- peak materialised storage below `2^50` bytes;
- measured complete time below matched rho in the projected memory tier;
- no discarded probabilistic branch counted as exhaustive coverage.

Attempt no full-depth unplanted P-256 relation unless all gates pass.  If pair
construction crosses rho only after unmeasured multi-terabyte routing, publish
that split result and keep parity unestablished.
