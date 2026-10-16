# Validation outcomes and scope

The graph package's final source hash is frozen in every final plan and prefix
replay receipt. `native-tests.txt` contains six passing unit tests and two
passing binary tests. No failed outcomes were reclassified as passes.

`isogeny-cli-tests.txt` records the separate installed-tool suite: 95 passing
tests and one pre-existing ignored test. The current `*-layout-*.json` receipts
show `crypto isogeny class` constructing a complete mixed-degree control and
`isogeny class` replaying it. Both command adapters and the standalone command
produce the same final source fingerprint for the fixed large input plan.
Earlier `*-entry-*.json` receipts are development snapshots from before the
storage metadata was renamed to describe its exact three-u64 layout; their
older source hashes are preserved and are not final-source receipts.

## Whole-library check remains a combined-run failure

The required `cargo test --release --lib` invocation ran to completion:
**4,166 passed, 2 failed, 117 ignored**. The two failures are in the unchanged
`src/cryptanalysis/jv_cover.rs`:

- `nagao_finds_planted_six_sums_and_agrees_with_the_oracle`: its recorded first
  trial exceeded the 120-second F4 budget (`timed_out: true`).
- `replayed_f4_agrees_with_the_full_solver_and_the_oracle`: the count of eligible
  completed cases was zero; the assertion requires at least 30.

Both tests pass individually with the same compiled test executable and unchanged
source. The complete failing run and both successful reruns are retained here.
The source file's SHA-256 is
`a6d69fc93e10329b9608b0af6898fbd67867ace64d149e716d6830475424ca42`,
identical to the base revision's file. This supports a concurrency-sensitive
failure diagnosis but does not establish a green combined suite or a separately
executed clean-base suite. No unrelated solver code, budget or ignore attribute
was changed. The full run preceded the final catalogue/metadata refresh;
the refreshed catalogue is checked separately by the native cover and view checks.

These are correctness-control logs. Their printed internal durations/costs are
not admitted performance measurements, and are not used for a speed comparison.

The global curve-name lint reports five retired-name occurrences in two
unchanged archival result files. The changed-line check reports zero problems;
both outputs are retained. No frozen archive was rewritten to clear that gate.

## Catalogue checks

The native cover generator and independent `--check` invocation both report
2,120 verified models, zero invalid inputs and zero unsupported models. The
browser and leaderboard `--check` outputs are retained after materializing
their required frozen metadata in this sparse checkout. Model registration used
the isolated native metadata importer; it does not enable order-preserving JSON
maps in the cryptographic crates or change their canonical certificate bytes.

Conductor is unavailable (exit 127), as recorded in the protocol and problem.
The protocol was written in the working tree before these runs; its first Git
commit is this delivery, not an earlier preregistration commit. These runs make
no preregistered performance-improvement claim.
