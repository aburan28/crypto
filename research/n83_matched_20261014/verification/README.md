# Native fixture and answer verification

The frozen fresh panel contains 32 synthetic public targets on
`icv1-f2m83-t6151469093347-9eba6cd5`: `a=b=1`, field polynomial
`z^83+z^7+z^4+z^2+1`, subgroup order `8569786107849059`, generator
`(6614244651398184860292015,8868260410029249391539529)`.
The files are for separate one-target experiments. The combined file is an
ordered container and does not authorize shared-target rho state.

[generate_targets.cpp](generate_targets.cpp) and
[check_fixtures.cpp](check_fixtures.cpp) use standalone C++ bitwise field and
affine group arithmetic rather than the Rust workers' arithmetic. The checker
reproduces all draws and validates every accepted point, its subgroup membership,
its validation-only scalar, its exact bytes and its signed Frobenius orbit.
Both programs reuse the historical independently implemented arithmetic in
[native_arithmetic.hpp](native_arithmetic.hpp),
[native_helpers.hpp](native_helpers.hpp) and
[fresh_fixture_helpers.hpp](fresh_fixture_helpers.hpp); those files were copied
without altering the historical studies. This is independent solver-answer
verification, not a second independently implemented fixture generator.

The fixed SplitMix64 seed is `202610146301`. A 64-bit rejection step produces
unbiased nonzero scalar proposals, followed by rejection of excluded signed
Frobenius orbits. The initial 617 excluded orbits are the 600 K600 representative
orbits, the original public target and all 16 targets from the two preceding
optimization studies. Both preceding eight-target draw streams are reproduced
before exclusion. Each new accepted orbit is also excluded from later proposals.
The accepted point distribution is conditioned on these exclusions. A fixed
seed supplies reproducibility, not a proof of randomness.

All 32 proposals were accepted, with zero range or orbit rejections. The
[generator receipt](fixtures/generator_receipt.json) binds its inputs, point
files, sidecar, executable and source. The
[fixture replay receipt](fixture_replay_receipt.json) passed all checks,
including all 99,600 K600 signed points and labels.
[matched_fixture_helpers.hpp](matched_fixture_helpers.hpp) contains the new
617-orbit exclusion and 32-target replay logic.

Only `fixtures/target_01.jsonl` through `fixtures/target_32.jsonl` are one-target
worker inputs. Every file contains exactly one coordinate pair, encoded as
decimal strings. `holdout_32_targets.jsonl` contains the same point-only rows in
01-through-32 order. The validation sidecar contains the expected scalars and
the complete draw/rejection trace; it is checker-only and must never be passed
to a measured child. No scalar appears in a point file. Output directories,
generated files and verification receipts use exclusive creation.

The [four negative checks](negative_checks/receipt.json) rejected an existing
output directory, an existing verification receipt, changed point bytes and an
incorrect scalar. For the incorrect-scalar case the sidecar's recorded hash was
updated first, so rejection came from the independent arithmetic scalar replay,
not simply a stale file hash. The retained shell script only prepares file
copies and launches native programs; it performs no research arithmetic or
statistical calculation.

Exact build and generation commands, run from the repository root:

```sh
taskset -c 0 g++ -std=c++17 -O2 -Wall -Wextra \
  research/n83_matched_20261014/verification/generate_targets.cpp \
  -o /workspace/n83_matched/verification/generate_targets
taskset -c 0 g++ -std=c++17 -O2 -Wall -Wextra \
  research/n83_matched_20261014/verification/check_fixtures.cpp \
  -o /workspace/n83_matched/verification/check_fixtures
taskset -c 0 /workspace/n83_matched/verification/generate_targets \
  experiments/koblitz-single-target-n83-20261006 \
  research/n83_evaluation_20261012/optimization/verification/fixtures \
  research/n83_speed_20261012/verification/fixtures \
  research/n83_matched_20261014/verification/fixtures
taskset -c 0 /workspace/n83_matched/verification/check_fixtures \
  experiments/koblitz-single-target-n83-20261006 \
  research/n83_evaluation_20261012/optimization/verification/fixtures \
  research/n83_speed_20261012/verification/fixtures \
  research/n83_matched_20261014/verification/fixtures \
  research/n83_matched_20261014/verification/fixture_replay_receipt.json
taskset -c 0 bash \
  research/n83_matched_20261014/verification/run_fixture_negative_checks.sh
```

The commands above have completed. Regeneration and replay require a **new**
output directory and receipt path. Executables are local scratch artifacts;
their hashes are retained in [native_executable_hashes.sha256](native_executable_hashes.sha256).
These preliminary receipts establish fixtures only. They do not establish a
solver result, performance comparison, independent-host replay or common
operation accounting. Final measured-answer and file replay will have separate
receipts after the campaign finishes.

## Final campaign checker, implemented and compiled

[verify_matched_campaign.cpp](verify_matched_campaign.cpp) reconstructs the
prospective execution matrix independently: 32 processes, 218 target rows,
104 IC scalar answers, 114 rho scalar answers, and 104 four-point IC witnesses.
It binds each ordered target, seed, point file, configuration and executable,
checks the recorded clean environment and common resource envelope, compares
observations with their raw records and costs, and checks actual file hashes.
Rho's three embedded compiled-source bindings and any relevant-source bindings
explicitly frozen in the protocol are checked without claiming a complete
dependency closure. Post-run snapshot hashes are distinguished from bindings
recorded before or during execution.

For each of the five A/A pairs the checker requires identical points, walk
and jump seeds, policy, recovered scalar, recorded outcome and preparation
charges. This is equality of exported outcomes and counted work, not a proof
that every intermediate state is identical. The success receipt includes
`checked_aa_semantic_pairs=5`, `aa_semantic_gate=true`, `checked_target_rows=218`
and `invalid_rows=0`, bound to the campaign manifest and frozen protocol.

Missing, failing or incomplete process/target rows remain in the report. A
complete `PASS` requires all planned processes and answers. Returned scalars
must be canonical nonzero subgroup residues and pass independent C++ group
replay before the validation sidecar is checked. IC decomposition witnesses,
abscissae and pinned intermediate points are also independently replayed.
This verifies returned answers, not the complete IC rank and factor-log
transcript. IC's five reported timing fields are checked with their sum and
residual retained; the current input-validation mismatch is not erased.

After all reserved solver timing has ended, run this checker through a
separately metered native isolated-bench invocation:

`verify_matched_campaign ARCHIVE PRIOR_OPTIMIZATION_FIXTURES PRIOR_SPEED_FIXTURES MATCHED_FIXTURES CAMPAIGN FROZEN_PROTOCOL NEW_RECEIPT`

The standalone executable is
`/workspace/n83_matched/verification/verify_matched_campaign_v2`.
[The separate verification revision](revisions/20261014-v2/README.md) explains
the strengthened preparation/timer bindings and launch-failure classification,
with preserved old/new source and executable hashes. The old
[final_checker_hashes.sha256](final_checker_hashes.sha256) remains the original
pre-revision pin; the measured workers and fixtures remain unchanged. The final
compile is warning-free; the earlier warning and correction remain in the build
logs. Final answer/file replay has not yet
run. Its external audit wall, CPU and RSS belong to a separate receipt and
are outside the already charged built-in solver verification intervals.

## Singleton rho reshape preflight

[check_rho_preflight.cpp](check_rho_preflight.cpp) is compiled and ready to
compare the old singleton worker with the new reused-setup worker on the fixed
historical speed-study `holdout_1`, using seed `202610126241` and the selected
Fold/reference/256/DP6/escape-off policy. It compares point, scalar, seeds,
policy, all nine common outcome fields, preparation charges and total charges;
it derives the old table-payload lower bound from the new table-entry count.
Both returned scalars also pass independent native arithmetic replay.

`check_rho_preflight OLD_RECEIPT NEW_TARGET_JSONL NEW_RECEIPT` runs the check.
Adding `--self-test` before those arguments requires rejection of copied
records with a changed walk counter and a changed scalar. Neither invocation
launches a solver. The actual old/new solver preflight and semantic receipt
are still pending; their measurements are excluded from the fresh-panel
statistics. This check establishes exported singleton semantics only, not
every trajectory state, every target or a performance result.
