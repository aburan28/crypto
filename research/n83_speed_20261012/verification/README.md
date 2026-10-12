# Independent native fixture and answer verification

Generation, fixture replay and all 139 answer replays are complete. The generator and checker use
standalone C++ field and group arithmetic, not the Rust solver's arithmetic.
All input files are public synthetic examples. This directory makes no timing
or speed claim; the campaign controller supplies performance evidence separately.
The [final verification report](FINAL_VERIFICATION_REPORT.md) links the complete
[answer receipt](final_receipt.json), [file provenance audit](file_provenance_receipt.json)
and controlled negative checks. All answer, count and nine recorded-semantic
pair gates passed; all 139 file and resource records passed the separate audit.

The exact curve is `icv1-f2m83-t6151469093347-9eba6cd5`, with
`y²+xy=x³+x²+1` over `F₂[z]/(z⁸³+z⁷+z⁴+z²+1)` and subgroup order
`8569786107849059`. Its generator is
`(6614244651398184860292015,8868260410029249391539529)`.
Neither the primary a=0 model nor the other degree-83 modulus is substituted.

The reusable arithmetic comes from the previous optimization's independent
checker. [native_arithmetic.hpp](native_arithmetic.hpp) contains only its field,
group, JSON, primality and modulus helpers; the previous checker's main function
and its historical metadata are excluded. [native_helpers.hpp](native_helpers.hpp)
checks the exact archived fixture and all 600 representatives, 99,600 signed
Frobenius orbit points and labels. [fresh_fixture_helpers.hpp](fresh_fixture_helpers.hpp)
then reproduces the previous eight optimization targets, verifies their file and
sidecar hashes, and excludes their entire signed Frobenius orbits.

[generate_targets.cpp](generate_targets.cpp) uses SplitMix64 seed `202610126201`,
unbiased modulo conversion of 64-bit pseudorandom proposals to `[1,r−1]`, and
rejection of every excluded orbit. The frozen stream accepted eight proposals
in eight draws, with zero range or orbit rejections. The resulting three training
and five holdout files are in [fixtures/](fixtures/). The accepted point law is
conditioned on exclusion of the original target, K600, previous optimization
targets and earlier accepted targets; it is not an unconditional point sample.
The seed makes the experiment reproducible, not a proof of randomness.

Only the eight one-point `*.jsonl` files are solver inputs. The combined
`holdout_5_targets.jsonl` is an archival convenience. The
`validation_sidecar.json` contains expected scalars for independent validation
and must never be passed to a measured solver. The generator uses exclusive
creation for the output directory and every output file; the checker refuses
to overwrite a receipt.

The [generation report](GENERATION_REPORT.md) and
[exclusion diagram](generation.svg) explain the invariants. The
[fixture replay receipt](fixture_replay_receipt.json) passed curve, subgroup,
scalar replay, signed-orbit disjointness, exact draw-stream reproduction and
file binding for all eight fresh targets and all eight old targets.
The generator receipt SHA-256 is
`2b2a4987a5f93e751f5a28eb7856d31f8833a5b4d55f9baa92eb2c5ec427e0ae`.

Build and fixture replay, from the repository root:

```sh
mkdir -p /tmp/n83-speed-check
g++ -std=c++17 -O2 -Wall -Wextra \
  research/n83_speed_20261012/verification/generate_targets.cpp \
  -o /tmp/n83-speed-check/generate_targets
g++ -std=c++17 -O2 -Wall -Wextra \
  research/n83_speed_20261012/verification/check_speed.cpp \
  -o /tmp/n83-speed-check/check_speed
/tmp/n83-speed-check/check_speed --fixtures-only \
  experiments/koblitz-single-target-n83-20261006 \
  research/n83_evaluation_20261012/optimization/fixtures \
  research/n83_speed_20261012/verification/fixtures \
  /tmp/n83-speed-check/new-fixture-receipt.json
```

The already built scratch binaries are
`/workspace/n83_speed/verification/{generate_targets,check_speed}`. Their hashes
are recorded in the generation and replay receipts. Binaries are not committed.
Regeneration needs a new output directory; use the same first two inputs and
pass that new directory as the generator's third argument. Rebuilt generator
receipts have different executable/path-dependent source hashes; point and
sidecar bytes remain determined by the fixed inputs and seed.

Full panel replay after the 139 processes complete:

```sh
/tmp/n83-speed-check/check_speed \
  experiments/koblitz-single-target-n83-20261006 \
  research/n83_evaluation_20261012/optimization/fixtures \
  research/n83_speed_20261012/verification/fixtures \
  /workspace/n83_speed/panel 139 \
  /tmp/n83-speed-check/new-panel-receipt.json
```

[check_speed.cpp](check_speed.cpp) recursively checks each `execution.json`
against its `producer.json`. It binds each one-point input, seed and policy
file; checks curve, generator, subgroup, lambda, scalar-input flags and group
verification flags; and replays every returned scalar with independent C++
arithmetic. It checks the complete 10 A/A, 18 arithmetic, 81 policy and 30
holdout process/answer counts and per-fixture counts. It also checks the nine
arithmetic pairs for equal answer, jump seed, outcome counters, charges and
rho policy. These exported equalities do not establish every intermediate state.

All missing receipts, no-answers, caps, nonzero exits, invalid answers and
unpaired outputs remain in its report. A full success requires every planned
process to complete and return an independently verified scalar. Failure and
cap preservation never converts them to successful answers. Independent replay
is outside all measured solver intervals. The checker does not establish
independent timing reproduction, a common calibrated operation unit, host
isolation, independent-host reproducibility or an admission gate.

The separate [file provenance auditor](verify_file_provenance.cpp) checks actual
bytes against recorded SHA-256 bindings after answer replay. It binds the solver
executable, every Rust source named in the compiled-source receipts, the original frozen protocol, point
and configuration files, observations, producer output, stdout and stderr. It
also checks all four frozen stage matrices, immutable training-selection inputs,
recorded clean environment/resource settings, CPU component sums and RSS unit
conversion. It snapshots every retained panel file, including heartbeat and
isolation logs, and distinguishes new snapshot hashes from pre-existing recorded
bindings. JSON stage protocol copies must match the frozen protocol structurally;
they need not preserve its key order or original serialization bytes.

```sh
g++ -std=c++17 -O2 -Wall -Wextra \
  research/n83_speed_20261012/verification/verify_file_provenance.cpp \
  -o /tmp/n83-speed-check/verify_file_provenance
/tmp/n83-speed-check/verify_file_provenance \
  /tmp/n83-speed-check/new-panel-receipt.json \
  /workspace/n83_speed/protocol.json \
  f5aa7d7fe2e58ddcfb8177dcc9d9907f4c5aeb859cadea983448b872b1f05d7c \
  /tmp/n83-speed-check/new-file-provenance-audit.json
```

This file audit verifies byte bindings and recorded metadata, not independent
timing reproduction or the truth of unobservable physical-host conditions.
