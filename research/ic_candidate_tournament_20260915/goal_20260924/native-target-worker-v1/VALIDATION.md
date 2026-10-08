# Native target worker source and file controls

This is a source-integration result for the disclosed educational n17 fixture,
based on PR #1362 at `30c84566f782e6503901d187348709dfaad05cc2`.
No scientific preparation panel, F5/SAT search, target worker `run`, registration
freeze, target exposure or comparative benchmark was executed. Native execution
admission, natural relation yield, fresh qualification and speedup remain unknown.
The three old confirmation sets and consumed control registrations remain closed.

## What changed

`prepared_target_worker` connects the bounded F5 recording API to new exclusive
start/completion files. It requires its own consumed target claim, exact original
worker path and byte hash, compiled source/build identity, immutable inventory,
bounded configuration and frozen single-thread environment. Before online timing,
it invokes the original preparation capsule's frozen native auditor and binds its
complete rank-29 F5 input to the registered producer and mathematical hashes.
It never runs an archive restoration or an old preparation worker. This production
integration has compiled but has not run against a scientific capsule.

The attempt journal fixes target, backend family, seeded coefficients and query
cap. It syncs the new directory and exclusive files; completions link to their
start. Validation or write failure permanently stops that writer. Bounded
inspection rejects gaps, orphans, altered links and unknown members, and preserves
an unfinished start without admitting a result. It accepts the F5 and SAT record
shapes; this does not implement a complete SAT worker or verify a SAT model.
Actual record writes stay inside the existing F5 online interval. Setup, input
loading and original preparation audit stay outside it. A producer without a
verified scalar exits unsuccessfully after retaining its report.

## Retained local validation

| Evidence | Outcome | Meaning |
| --- | --- | --- |
| `source-controls-v1.log` | Compile failure, E0603 | Private SHA helper import; no tests ran. Fixed by importing the public helper. |
| `source-controls-v2.log` | 23 passed | Initial successful controls; original files retained. |
| `source-controls-v3.log` | 23 passed, exit 0 | Final source controls, including the replayable watchdog publication layout. |
| `clippy-v1.log` | Exit 0 | No warning attributed to the new worker files; inherited crate warnings and an unknown configured lint remain visible. |
| `worker-build-v1.log` | Exit 0 | Final non-test development binary built offline with the locked dependencies. |
| `development-worker-identity.json` | Source/build identities null | This binary cannot satisfy scientific dispatch identity. |
| `prefix-cli-replay-v1.json` | Byte-identical to `watchdog-control-v3/prefix.json` | Actual non-test `inspect` CLI reproduces the published structural result. |
| Explicit five-file rustfmt check; workflow Ruby YAML parse; git whitespace check | Passed | Formatting and workflow syntax only. |

The 23 tests comprise 11 meaningful new contract/journal controls, one inactive
child helper in the parent suite, seven reused native file/transport/asset controls
and four reused preparation-contract controls. They are not 23 solver executions.
Asset extraction reads accepted archives as data and launches no archived binary.
Both families' journal payloads and the preparation-receipt control are synthetic.
Receipt syntax alone cannot establish native or mathematical admission.

The final watchdog control ran a file-only child. It wrote and synced
`attempts/start-000.json`, flushed `DURABLE_START_READY`, then reached the 2,000 ms
deadline. The retained original receipt reports `timed_out: true`, identical
before/after test executable hashes, and confirmed process-group drain. Inspection
retains trial 0 as pending and zero completions. Original receipt, stdout, stderr,
record bytes, inventory and structural output are in `watchdog-control-v3/`.
`watchdog-control-v2/` preserves the preceding control's different publication
layout. This tests abrupt process termination; it is not power-loss testing or an
interrupted cryptographic run. Cryptographic searches executed by the child: zero.

Builds, tests, clippy, CLI inspection and parsing used the shared
`/Volumes/SSD990/llm/tmp/ic-native-busy-20261002 busy --` wrapper. Rust/Cargo were
Homebrew 1.93.1 on macOS ARM64. Incremental compilation and development/test debug
information were disabled, and `IC_TARGET_SOURCE_MANIFEST_SHA256` and
`IC_TARGET_BUILD_SHA256` were explicitly unset for the final controls and binary.
The ignored development Cargo.lock was unchanged from the dependency worktree;
it is not added to Git. Busy serialization does not establish host isolation.
Linux controls are configured in the existing disclosed replay CI job; local
Mac results do not prove a Linux pass.

For a separate control replay, use a new receipt directory and retain its output:

```sh
/Volumes/SSD990/llm/tmp/ic-native-busy-20261002 busy -- \
  env -u IC_TARGET_SOURCE_MANIFEST_SHA256 -u IC_TARGET_BUILD_SHA256 \
  TMPDIR=/Volumes/SSD990/llm/tmp \
  CARGO_TARGET_DIR=/Volumes/SSD990/llm/tmp/ic-native-preparation-build-cache-20261004 \
  CARGO_INCREMENTAL=0 CARGO_PROFILE_DEV_DEBUG=0 CARGO_PROFILE_TEST_DEBUG=0 \
  IC_TARGET_RECORD_RECEIPTS=/Volumes/SSD990/llm/tmp/NEW-target-file-control \
  /opt/homebrew/bin/cargo test --locked --offline \
  --bin prepared_target_worker -- --test-threads=1
```

This command runs file/transport controls only. Do not substitute worker `run`,
an old scientific registration, or an archive executable. Prefix CLI replay used
`inspect --execution <watchdog-control-v3> --registration-sha256 <64 ones>` and a
new output outside `attempts/`; that synthetic seal grants no scientific authority.

## Exact final source and evidence digests

All hashes below are SHA-256, over the retained original bytes.

| File | Digest |
| --- | --- |
| `src/bin/prepared_target_worker.rs` | `a793e1c77e4e8dfd06ce0ffae46501684909c42a7f1443db5ebf9e00fabe51f8` |
| `src/bin/prepared_target_worker/contract.rs` | `dd8e489565204f07df781b586f0a5bf3f5da40ca3f440f14a5aa9e068bae3797` |
| `src/bin/prepared_target_worker/journal.rs` | `b5a89fc84ddfb9bf803c3a5738021ffd48a2a37fc606f3aa4f4b412f5222e735` |
| `src/bin/prepared_target_worker/contract_tests.rs` | `7dbaebd6313f2031bc3a34617872fec6a3e80e76d4a2a6495b1cedaa913a1aac` |
| `src/bin/prepared_target_worker/journal_tests.rs` | `41c6df45c3ba93a0e851785cd72bb7f7371b918756ff08be76647f67ab7e5982` |
| `source-controls-v3.log` | `f970de7c1fe6235919187695e88393456f4cdcb9355e4714a15274f2ae510810` |
| `clippy-v1.log` | `89455ee1abca81fd28972ddb84bd7898bb3144045a6e0f9a6ec8429ed0e0041b` |
| `worker-build-v1.log` | `a80e3f85680b6ad52ab075e8afaed2039f1024bc07930ebdc74ace93f19715ac` |
| `watchdog-control-v3/control.receipt.json` | `e14dd8b19fc9f15ce5f43247e83c74c509271bb266e698b5e3fc2735cdad9157` |
| `prefix-cli-replay-v1.json` | `311946e221d186401cdf7881627e54dafd4d0dbc0f9b078308c29f5aeb6cc1e1` |
| Development non-test binary, recorded metadata only | `d44fb8e95d421fdd61b012737c92ae31ea49ce152d4d9ee6fec285389e5911fe` |
| Ignored development Cargo.lock | `f99127c279e83c5fd8a474a604e5457a96d0ed16fb4c84758c78f51f1291ce79` |

The development binary is not published as a scientific executable. Its identity
and file hash document only the build and data-reader checks.

## Admission still pending

The new target controller/freezer and independent complete-run auditor are not
implemented by this PR. New native source-admitted ordinary preparation for F5
and CryptoMiniSat, their public external seals and one-use execution receipts,
target exposure/custody and independent failed-attempt accounting remain required.
SAT also needs accepted preinitialized transport; a cold per-query CLI launch or
subtracted startup interval cannot silently answer the one-target online question.
Canonical candidate/workload/run identities, fresh paired incumbent and strong
same-point rho, operation/memory accounting, and host isolation/noise receipts
remain required. This is engineering progress; full goal complete is false and
online speedup is null. Dependency CI failures are retained in
[CI_SOURCE_CLOSURE.md](CI_SOURCE_CLOSURE.md).
