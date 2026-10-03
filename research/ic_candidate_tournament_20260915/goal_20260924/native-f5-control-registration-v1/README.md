# Native prepared F5 control preregistration

**Status: registered locally, not dispatched; source acceptance and this
preregistration's final-head CI remain pending.** This is the actual frozen
capsule for the disclosed synthetic n17 known-input control. It is separate from
the non-executable compact build-validation capsule. No scientific invocation
claim exists, no new query has run and no target has been recovered by this
registration. This historical stage label must not authorize later archive
restoration or retry after consumption.

The [controller](../native-f5-control-v1/PROTOCOL.md) is in
[PR #1284](https://github.com/aburan28/crypto/pull/1284); its
[full custody implementation](../native-f5-custody-v1/PROTOCOL.md) is in
[PR #1294](https://github.com/aburan28/crypto/pull/1294). The exact clean source
snapshot is `3eb99b4d2d934059834174a4d75789bf410b916a`. Both implementation PRs
must be accepted, and all applicable exact-head checks on this preregistration
must pass, before its sole scientific dispatch. Old registrations and the three
historical confirmation sets stay consumed and closed.

| Frozen item | Value |
| --- | --- |
| Canonical registration SHA-256 | `5a50e01c02fa7e14f227dea413cb9e569f73d68d684a68541e49fd357660f060` |
| Worker SHA-256 | `c4a2f160faeb57a2451171a2f105c05c77d57abad76e13eb8112f57ac47c20b0` |
| Frozen dispatcher / independent checker SHA-256 | `5055542e85dc805c6a0c8ee9e09548d53f8bf2a84892e15e02695d580221e8b1` |
| Published archive | 28,451,993 bytes; SHA-256 `4cf00be22f9a5eb5084e4374cf8e211dedf2a1953ba84dabf84350579d46c5d0` |
| Source/build identities | `fd62f42c05b548d76af959a7d88b84535b20d3b0a2d7acd3dad8957f6dcc6b03` / `d52c569b427027e01e7485e628b09b3f71e25561e2ca55dd00f1f0d465a7d832` |
| Synthetic curve / supplied point | `EC1N17Ckb1hbbe2b5b6b1e6` / `[52411,72106]`, previously disclosed |
| Seed / maximum target queries / outer deadline | `2026100310` / 8 / 900 seconds |
| PDP / internal kernel | Three-summand bounded MatrixF5, degree 3, reduction budget 8192, dense Boolean rows |
| Base / relation columns | 63 geometric points, 62 subgroup-usable points, 29 folded columns |
| Immutable files / full archive regular files | 5,954 / 5,960 |
| New ordinary-query executions | 0 |
| New complete IC1 candidate/workload/run identifiers | Pending; this is a disclosed control registration |
| Execution/runtime/fresh qualification, headline and promotion | False |
| Paired rho/incumbent costs, online time and speedup | Unknown |

The [original registration](publication/registration.json),
[seal](publication/seal.json), [config](publication/config.json),
[whole preparation](publication/preparation.json),
[mathematical-only worker input](publication/job.json) and
[physical host context](publication/host-context.json) are copied without
alteration. [Publication metadata](publication/PUBLICATION.json) binds them to
the full [capsule archive](publication/capsule.tar.gz). The archive retains Rust
sources, the actual worker example, exact Cargo manifest/lockfile, compiled-in
resources, offline vendored dependencies, both release binaries, original build
logs/receipts and compiler identities. Build caches are excluded.

The offline build used native macOS ARM64 Rust/Cargo 1.93.1 with an empty Cargo
home, explicit compiler and fresh build target. Source/build identities were
compiled into the worker and verified by its only call, `--build-identity`, with
no job input. No scientific solver ran. The host is the physical Apple M4 Pro
machine in the retained host record; this is not a calibrated performance
environment. The build log retains the xcrun temporary-directory warning and
successful exit rather than hiding it.

Native preparation replay reconstructs all historical 216 ordinary attempts,
61 witnessed relations, 155 complete geometric negatives, rank 29 and all column
logs. Original Python preparation provenance remains explicit. These retained
records are not new native natural yield. The engine is bounded Macaulay
MatrixF5 with its F5 row criterion, not a complete incremental Gröbner algorithm.
The worker receives geometry, logs, public point and limits; it receives no
historical target attempts, witness, recovered scalar or certificate proof.

## Archive custody and preserved failure

The first macOS stdout archive had zero padding outside its gzip stream. The
native checker rejected it with `invalid gzip header`; system gzip reported
trailing garbage. [The full first archive and original sidecars](failed-packaging-v1/FAILURE.txt)
remain preserved. A new named-file archive passed unchanged native custody
verification, retained in [the local receipt](local-macos-arm64-custody-v2.json).
No source, seed, target, registration or claim changed, and no gzip check was
relaxed. Linux/macOS CI independently checks the accepted-format archive and
requires rejection of the first packaging attempt. This corrects publication
formatting, not a mathematical failure or performance result.

Custody checks the external registration seal, exact sidecars and every archived
file, with bounded decompression and strict USTAR/path/member/trailer rules. It
extracts and executes nothing. A portable custody pass verifies bytes only;
it does not validate the Mac binaries on Linux or admit a scientific execution.

```sh
/tmp/ic-native-busy busy -- target/debug/icprog f5-control-replay-custody \
  --publication research/ic_candidate_tournament_20260915/goal_20260924/native-f5-control-registration-v1/publication \
  --registration-sha256 5a50e01c02fa7e14f227dea413cb9e569f73d68d684a68541e49fd357660f060 \
  --out /tmp/native-f5-registration-custody-NEW.json
```

## Sole execution after acceptance

The live original capsule is
`/private/tmp/ic-native-f5-control-20261003-preregister-v1`. Before dispatch,
confirm accepted source PRs, durable publication, passing applicable final-head
CI/review and absence of a live invocation claim. No new path, restored copy,
modified budget or replacement seed authorizes another invocation.

```sh
/Volumes/SSD990/llm/tmp/ic-native-busy-20261002 busy -- \
  /private/tmp/ic-native-f5-control-20261003-preregister-v1/immutable/bin/icprog f5-control-execute \
  --capsule /private/tmp/ic-native-f5-control-20261003-preregister-v1 \
  --execution /private/tmp/ic-native-f5-control-20261003-execution-v1 \
  --registration-sha256 5a50e01c02fa7e14f227dea413cb9e569f73d68d684a68541e49fd357660f060
/Volumes/SSD990/llm/tmp/ic-native-busy-20261002 busy -- \
  /private/tmp/ic-native-f5-control-20261003-preregister-v1/immutable/bin/icprog f5-control-audit \
  --capsule /private/tmp/ic-native-f5-control-20261003-preregister-v1 \
  --execution /private/tmp/ic-native-f5-control-20261003-execution-v1 \
  --registration-sha256 5a50e01c02fa7e14f227dea413cb9e569f73d68d684a68541e49fd357660f060 \
  --out /private/tmp/ic-native-f5-control-20261003-independent-audit-v1.json
```

The synchronized create-only claim consumes the invocation before launch,
including setup error, timeout, crash, interruption or failed recovery. Preserve
all raw output, delivery/watchdog/PID/drain receipts and terminal failure; never
resume, retry or extend this registration. Original admission requires the
preexecution-frozen checker and the same externally published seal. It starts
no solver and independently checks every attempt, negative, witness and scalar.

All producer wall values will remain uncalibrated diagnostics. External audit
time is separate; no matched reference, resource/noise calibration or fresh
target is registered here. A successful known-input control cannot establish
the primary single-target speedup or a tournament winner. The persistent goal
still requires new native natural-yield evidence with failed attempts, an
exposure census, exact IC1 method/workload identities and a new frozen paired
one-target comparison with the incumbent and matched rho. This preregistration
neither generates fresh targets nor reopens prior confirmation data.
