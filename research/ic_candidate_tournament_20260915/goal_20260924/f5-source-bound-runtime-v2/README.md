# F5 v2 adapter: interface and publication prerequisites

The [new protocol](PROTOCOL.md) fixes one bounded n17 F5 development invocation
with prospective seed 2026093032. It has **not been dispatched**. The consumed
[v1 invocation](../f5-source-bound-runtime-v1/README.md) is closed after its input
schema failure; its raw evidence is preserved in merged
[PR #1058](https://github.com/aburan28/crypto/pull/1058), merge
`5885d629d8c30261ee8c83157fa51792406545c0`.

[PR #1060](https://github.com/aburan28/crypto/pull/1060) corrects two independent
adapter defects. It sends `target_seeds: []` to the actual native `Vec<u64>`
input while retaining `[null]` in output fixture provenance. It also publishes
finite raw resource seconds as roundtrip decimal strings and retains/hashes the
original metrics bytes. Native scientific clocks and phase costs remain integer
nanoseconds; resource diagnostics do not become calibrated performance evidence.

The separately preregistered [native interface control](interface-control/README.md)
passed. Two valid inputs parsed with every expected field/default, and all six
invalid inputs rejected. It imports the actual private `Job`/`Config` types from
the exact pinned native source, exposes only a parser accessor, and never calls
`execute_job`. Archive/source/input/output replay passed without native execution.
The live v2 registrar requires the exact proof archive hash before freezing a
new invocation and binds that hash in the canonical candidate. This parser
control's direct Python source inventory is not complete execution attestation;
the future IC job separately requires complete controller/package/interpreter
and native pre/post gates.

Regression validation also found a shared supervisor cleanup bug: after a
successful watchdog group kill and leader reap, a redundant second signal can
produce macOS `EPERM` when only orphaned zombies remain. Cleanup now skips the
second signal after a successful kill or an already-exited group. Normal and
failed-controller cleanup still kills the inherited group, and errors on the
first kill still propagate. The actual inherited-native-child watchdog and a
post-reap permission-error regression passed. No old registration, source
manifest, phase record or sealed mathematical result was rewritten.

[controls/FINAL-FOCUSED.log](controls/FINAL-FOCUSED.log) retains 28 passing checks,
including actual watchdog/RSS controls, v2 source/identity/publication/proof
controls, and closed F5 failure and complete SAT archive replay. The initial
[sandboxed full regression failure](controls/FULL-SANDBOXED-FAILED.log) retains
two `/bin/ps` permission errors. The subsequent
[post-reap cleanup failure](controls/FOCUSED-POSTREAP-FAILED.log) is also retained.
These failures were not relabeled as algorithmic failures or silently dropped.
[controls/FINAL-FULL.log](controls/FINAL-FULL.log) records the final complete
regression: 370 tests, 367 passed and three declared macOS skips (strict frozen
Linux statistics and two Linux affinity/instruction-runner controls). The
repository-local skill also validates. Exact-head Linux CI remains a separate
merge gate.

The v2 controller remains development-only on the admitted physical macOS ARM64
worker. It has no incumbent or rho arm and no fresh-target qualification.
Natural relation yield, final rank, complete F4/F5 recovery, scientific online
cost, normalized S/floor and speedup for this v2 registration are unknown. The
next result PR must freeze the accepted snapshot and commit the external
invocation receipt before the single IC attempt, then preserve and independently
replay every terminal outcome. A parser pass or a merged adapter does not meet
the active goal: complete F4/F5 admission and a fresh paired incumbent/rho
comparison are still required. All three old confirmation sets remain closed.
