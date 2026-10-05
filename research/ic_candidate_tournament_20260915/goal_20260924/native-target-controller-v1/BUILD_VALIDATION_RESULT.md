# Validation-only target build and portable byte custody

This is a build-custody control on the disclosed educational n17 fixture. It
does **not** register or execute a scientific target solver. Its preparation
binding is explicitly synthetic and non-executable. The frozen worker was
invoked only for `build-identity`; publication and replay ran data checks.
The three closed confirmation sets remain closed. No new natural-relation
yield, one-target recovery, paired rho result, speedup or host-isolated CPU
measurement follows from this record.

The final source freeze used committed source
`5614c9a429becb77bceda469980ff1981225d1a1`, the checked-in
`control-v2/config.json`, the synthetic preparation pins from
`control-v2/registration.json`, Homebrew Cargo/Rust 1.93.1 on local macOS
ARM64, and an uncalibrated host-context record. The source was copied in full,
dependencies were vendored with `--locked --offline`, and the copied source
built `prepared_target_worker` and `icprog` in release mode. Original
vendor, compiler, build, worker-identity, environment and timing receipts are
in [build-validation-v1](build-validation-v1). The full original capsule and
outer command log are retained locally outside Git; the publication contains
every sealed immutable file and all small sidecars.
[Freeze result](BUILD_VALIDATION_FREEZE.json) and
[independent repository replay result](BUILD_VALIDATION_REPLAY.json) are
retained as machine-readable summaries.

| Evidence | Exact value |
| --- | --- |
| External validation registration seal | `571ec54d81e02ba68a6238236a51808631758a60f0dc997da684d398e6436662` |
| Source manifest SHA-256 | `4185dea51aef339b1fcb424fd6985c37e92ca483f74de9b5e6301cf71a48a44e` |
| Frozen worker SHA-256 | `5523d08f4d92e461640e5dd9c36a112c4f98c12d5d8b7c9409fd72d59b86b182` |
| Frozen checker/publisher SHA-256 | `82b61f35394500e22eb5cb3fecc3875d39b92c4c3af271873a1e276a125d6170` |
| Published Gzip archive SHA-256 | `d9f8b4db6cb6b9b9cad746f4b11fe439f19c476eaecd4c5df2c7cb64bf392ad3` |
| Sealed immutable files / complete archive regular files | 6,944 / 6,948 |
| Vendor/build/identity process exit codes | 0 / 0 / 0 |
| Validation-only / scientific worker calls / archived binaries executed | true / 0 / 0 |
| Source-bound execution admitted / headline speedup | false / null |

The frozen checker published and immediately replayed the full archive. A
second invocation replayed the copied repository artifact independently and
again returned `PASS_DATA_ONLY_TARGET_BUILD_CUSTODY`. Replay checks the exact
external seal, registration/config/host bytes, complete sealed source and
binary inventory, original build receipt set, tool and invocation paths,
environment, archive hash and every bounded USTAR member. It extracts or
executes no archived binary. In separate retained local negative controls, a
one-byte archive change failed with `target custody archive bytes differ`;
removing the worker-identity receipt and refreshing the outer descriptor still
failed with `target receipt set differs from original sealed inventory`.
The source tests additionally reject a modified host sidecar and attempts to
turn a validation publication into an execution or speedup claim.

Local targeted regression checks passed: two target-custody tests, 24 target
controller tests, 23 worker/contract tests and four prior ordinary archive
tests. The focused watchdog tests passed with debug info disabled. Changed
Rust files passed rustfmt; the workflow parsed as YAML; Clippy correctness
checking exited successfully with inherited warnings outside the new paths.
These are controls, not independent IC trials. The prior Linux replay failure
and passing corrected job are documented in
[CI_LINUX_DEBUG_BINARY.md](CI_LINUX_DEBUG_BINARY.md).

CI replays the checked-in archive as data on Ubuntu and macOS with the
external seal above and retains the replay JSON. It does not run its macOS
release binaries on Linux. The publication remains unconsumed and is never a
scientific execution registration. A future one-target study still needs
original native F5 and SAT preparation/target runs, canonical candidate and
workload manifests, natural failed-attempt/yield accounting, fresh matched rho,
and the physical-host isolation receipt before a controlled CPU speedup claim.
