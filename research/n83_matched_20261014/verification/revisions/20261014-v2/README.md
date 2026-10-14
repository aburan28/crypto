# Verification-only revision, 14 October 2026

Measured worker, controller and fixture source freeze remains
`c6e8daafb5a59d27e1be88a2f61f0bd6ee83ef03`. The parent authorized this separate
checker revision before main timing. No measurement worker, controller,
fixture, policy, seed, order or statistical estimator changed.

The read-only review found no matrix, schema or failed-row false success.
It identified two standalone checker coverage gaps that the controller already
covered: the two derived preparation split fields were not rebound to raw
summary values, and rho's millisecond total was not rebound to its integer
nanosecond total. A genuine launch failure's null CPU/RSS values were also
classified as invalid evidence rather than a retained incomplete no-child run.

The revised checker:

- Binds shared and per-target preparation to raw summary timings. IC uses
  total preparation and zero per-target preparation; rho uses total setup
  minus all G-only jumps, and all G-only jumps divided by the declared count.
- Binds rho target and online milliseconds to the integer nanosecond interval;
  all integer phase/residual identities remain exact.
- Accepts null CPU/RSS/exit metrics only for `LAUNCH_FAILED` with no semantic
  answers, records `child_launched=false`, and retains the incomplete outcome.
  Such a run can never pass the complete-panel gate. Actually launched children
  still require finite nonnegative CPU and exact RSS-unit conversion.

Floating bindings permit `8 * f64_epsilon * max(1, |a|, |b|)` after decimal
round-trip, conversion and subtraction/division. At 2,100 seconds this is below
0.004 nanoseconds in millisecond units; a one-nanosecond clock discrepancy is
rejected. The tolerance does not relax integer clock identities or permit a
missing timer field.

[native_self_test_receipt.json](native_self_test_receipt.json) passed no-child
and actual-child cases, four malformed/false-launched negative cases, and the
float-tolerance boundary. This is a small native checker test, not a solver
measurement. Build and test outputs are retained here; the build is warning-free.

[old_checker_hashes.sha256](old_checker_hashes.sha256) pins the preceding
source and executable. The preceding source bytes are retained in
[verify_matched_campaign_v1.cpp.snapshot](verify_matched_campaign_v1.cpp.snapshot),
and its executable remains `/workspace/n83_matched/verification/verify_matched_campaign`
(also copied to `verify_matched_campaign_v1`).
[new_checker_hashes.sha256](new_checker_hashes.sha256) pins the revised source
and executable `/workspace/n83_matched/verification/verify_matched_campaign_v2`.
Neither checker has yet replayed final campaign results. Use the v2 executable
for final replay, metered outside solver timing.

```sh
taskset -c 0 g++ -std=c++17 -O2 -Wall -Wextra \
  research/n83_matched_20261014/verification/verify_matched_campaign.cpp \
  -o /workspace/n83_matched/verification/verify_matched_campaign_v2
taskset -c 0 /workspace/n83_matched/verification/verify_matched_campaign_v2 \
  --self-test \
  research/n83_matched_20261014/verification/revisions/20261014-v2/native_self_test_receipt.json
```

The self-test receipt is create-new; reproduction needs a new output path.
