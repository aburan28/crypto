# Native prepared SAT control registration

**Terminal status: consumed and closed.** [RESULT.md](RESULT.md) records the
sole native execution, two independently proved negatives, source-verified
witness, recovered scalar and original frozen audit. The preregistration
description and commands below remain historical evidence; never execute,
retry, resume or extend this registration, including an archive restoration.

This is the executable preregistration for the bounded disclosed-input control
specified in [the accepted protocol](../native-prepared-sat-control-v1/PROTOCOL.md).
The native controller was accepted in [PR #1267](https://github.com/aburan28/crypto/pull/1267).
The actual clean source snapshot is merged commit
`d972396a0b476ab63889edcc4dd57745fc665bd2`, which also includes the accepted
ecbench online-window extension. The SAT controller sources are unchanged from
PR #1267. All three historical confirmation sets and all older consumed
registrations remain closed.

**Preregistration stage: no execution claim has been acquired and no exporter
or SAT solver has been launched.** This is a historical before-execution label,
not permission to replay the registration later. Inspect the external claim
and published terminal result before any dispatch; a restored archive is for
audit and must never substitute for a consumed invocation.

| Frozen item | Value |
| --- | --- |
| Canonical registration SHA-256 | `b098ebbe135abd4706e62d6286f513fe1c4f07fad61f7250487bd2d86ad317e4` |
| Native producer SHA-256 | `ec3208e308435ebb11a6b1e78c1589e825d8ee581cfc04f198f288668977aca1` |
| Frozen dispatcher / independent checker SHA-256 | `c5338662f688c15f52e2b02debf5e7f82ee9c5450a8de7b1ad0961225a27bf41` |
| Supplied synthetic target | `[52411,72106]`, previously disclosed |
| Query seed / exporter nonce | `2026100302` / `2026100303` |
| Target-query cap / conflicts per query | 8 / 1,000,000 |
| Exporter / solver / outer deadlines | 30 seconds / 60 seconds / 900 seconds |
| Immutable source, asset, binary and receipt files | 5,959 |
| Preparation | 63 geometric points, 62 usable points, 29 folded columns; rank 29 |
| New ordinary-query executions | 0 |
| Candidate, workload and run IDs for a new comparison | Pending; this is a control registration |
| Fresh qualification, headline and promotion | False |
| Paired rho/incumbent costs and speedup | Unknown |

The [actual registration](publication/registration.json) inventories every
immutable file; [the seal](publication/seal.json) fixes its canonical digest.
The [configuration](publication/config.json) and
[whole preparation](publication/preparation.json) are copied without alteration.
The [publication metadata](publication/PUBLICATION.json) binds the sidecar bytes
and complete compressed [capsule](publication/capsule.tar.gz).

The capsule retains the complete Rust `src` tree, manifest and lockfile,
compiled-in resources, offline vendored dependencies, both release binaries,
accepted native exporter/CryptoMiniSat assets and original build sources,
compiler/Cargo identities, environment, exact commands, build logs and
before/after executable digests. Release builds ran from the retained snapshot
with an empty Cargo home, explicit compiler and offline dependency resolution.
Build output/cache directories are excluded from publication; their immutable
receipts and final executable bytes are retained.

The accepted historical preparation is independently reconstructed before
registration and again before execution: all 149 ordinary records, 37 verified
relations, 106 complete geometric negatives, six inconclusive outcomes, rank
trajectory and 29 column logs. Its Python provenance and original failures stay
visible. These are existing ordinary-yield records, not new native yield.

## Publication verification

The native [custody auditor](../../../../examples/prepared_sat_registration_audit.rs)
checks the external registration seal, input/configuration hashes, executable
pins and every archive member. It rejects missing, duplicate, modified,
unregistered, linked and unsafe-path members, invalid headers and truncated
trailers. It reads the archive as data and extracts or executes nothing.
The [CI workflow](../../../../.github/workflows/ic-native-sat-registration.yml)
runs these controls and the complete publication audit on Linux x86-64 and
macOS ARM64. A custody pass proves retained bytes, not solver execution.

```sh
./target/debug/examples/prepared_sat_registration_audit \
  --publication research/ic_candidate_tournament_20260915/goal_20260924/native-sat-control-registration-v1/publication \
  --registration-sha256 b098ebbe135abd4706e62d6286f513fe1c4f07fad61f7250487bd2d86ad317e4 \
  --out /tmp/native-sat-registration-custody-new.json
```

Registration/build setup first encountered local storage exhaustion while
creating a separate checkout. Git removed that incomplete new checkout; no
capsule or scientific claim existed at that point. The clean existing checkout
was used instead, and only disposable task-owned incremental/build caches were
removed. Final capsule construction succeeded on the separate system volume.
The earlier capsule-v2 source-freeze validation remains unexecuted and is not
this registration; it predates final launch/admission guards and must not be
dispatched.

## Sole execution and independent audit

Execute only after this preregistration PR is accepted and all applicable
exact-head checks pass. The live original capsule is
`/private/tmp/ic-native-sat-control-20261003-final`. No new copy, restored archive,
changed point, seed, budget or output path authorizes a second invocation.
Use the exact frozen dispatcher under the native busy wrapper, with the
explicit seal and a new output directory:

```sh
/Volumes/SSD990/llm/tmp/ic-native-busy-20261002 busy -- \
  /private/tmp/ic-native-sat-control-20261003-final/immutable/bin/icprog sat-control-execute \
  --capsule /private/tmp/ic-native-sat-control-20261003-final \
  --execution /private/tmp/ic-native-sat-control-20261003-execution-v1 \
  --registration-sha256 b098ebbe135abd4706e62d6286f513fe1c4f07fad61f7250487bd2d86ad317e4
```

The create-only claim is synchronized before any child launch. A source gate
failure, timeout, crash or interruption consumes the invocation. Preserve the
claim, partial progress, all original source/process receipts and the terminal
outcome; never resume, retry or extend it. The frozen checker then audits into
a separate new report, without launching search:

```sh
/Volumes/SSD990/llm/tmp/ic-native-busy-20261002 busy -- \
  /private/tmp/ic-native-sat-control-20261003-final/immutable/bin/icprog sat-control-audit \
  --capsule /private/tmp/ic-native-sat-control-20261003-final \
  --execution /private/tmp/ic-native-sat-control-20261003-execution-v1 \
  --out /private/tmp/ic-native-sat-control-20261003-independent-audit-v1.json
```

Actual native solver execution is admitted only on macOS ARM64. Other hosts
may perform custody/math controls; their passes do not validate these native
assets on different hardware. The busy lock serializes heavy work; it is not
L2 timing isolation. All wall values are uncalibrated diagnostics, with external
independent replay outside the producer interval. A complete audited recovery
qualifies only this source-bound disclosed-input control. It is not a primary
online speedup, fresh relation yield or tournament winner.

The full F4/F5-plus-SAT goal remains active. It still requires complete native
F4/F5 control execution, new natural-yield evidence with failed attempts,
reconciled exposure exclusions, exact IC1 candidate/workload identities and a
fresh frozen paired comparison with strong compatible incumbent/rho references.
The merged ecbench online-window implementation is available for review and
integration; its public unplanted-target and IC1 identity gates are still open.
No historical confirmation set is reopened.
