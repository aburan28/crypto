# Canonical admission in both campaign drivers

The existing native screen and calibrated tournament now admit the reviewed
public-point producer before measurement. Each admitted arm/case freezes the
actual source/build/kernel, independent usable-base census, complete method,
canonical candidate and one-target workload. Successful runs reconstruct those
identities and independently verify the relations, rank, target descent, scalar
replay, native online interval and (when profiled) exclusive instruction ledger.
A frozen random run-number namespace prevents identifier collisions when the same
workload is executed again in another campaign. Early local controls retain their
original identifiers under their sealed contracts. Ignored configuration flags
cannot create a second IC competitor with the same canonical method.
Failed executions keep the same candidate/workload key and null verified online
timing and complete costs. Rho uses a distinct reference identity.

The per-target `single_target_online` table contains the public point, canonical
IDs, IC/rho online milliseconds, correctness and speedup. Cold process time and
instructions remain supplementary, separate measurements. Native screens retain
all reference widths. The instruction driver continues using the existing
portfolio, selection and replay stages, but admission controls alone cannot
enable promotion: comparative reference qualification and the frozen familywise
protocol are still pending. This work consumes zero improvement rounds.

The [protocol](PROTOCOL.md) was declared before the controls. The initial macOS
build failed because Cargo merges configuration arrays; the corrected driver
records an explicit host target and empty encoded-rustflags override. The next
preparation stopped because a strict canonical writer rejected a floating-point
process duration. Process measurements now have a separate raw JSON artifact.
Both failures remain in the archive.

Two subsequent fixed-input screens each verified 12/12 jobs. A deliberately
one-trial candidate then failed to collect enough relations; its six-job screen
retains five verified jobs and one failed canonical row. That row has no verified
online timing, instruction total or speedup. All 30 receipts revalidated after
fresh archive extraction with their own frozen evaluators. These are laptop
integration controls, not a comparative performance qualification.

Local evidence is in
[`ic-driver-native-controls-20260925.tar.zst`](../../evidence/ic-driver-native-controls-20260925.tar.zst),
SHA-256 `93808969fe84a1a3eaf169908f8295b1f763d09b32abb4bc3389d544845621be`.
It contains 2,698 files (48,193,193 expanded bytes), including exact source,
executables, evaluators, build logs, raw runs and failures. Build caches are
excluded. [History](native-history.json) and
[fresh extraction audits](native-archive-validation.json) are retained separately.

```sh
python3.11 research/ic_candidate_tournament_20260915/evidence/restore.py \
  --archive ic-driver-native-controls-20260925 --out /tmp/ic-driver-replay
python3.11 /tmp/ic-driver-replay/ic-driver-native-controls/native-v4/evaluator/autolab.py \
  verify --round /tmp/ic-driver-replay/ic-driver-native-controls/native-v4
```

`driver_smoke.py` is a bounded integration test of these same drivers, not another
tournament implementation. The producer workflow runs it on Linux amd64 against
the prepared scaled source and retains its artifacts, including failures. It
executes native screening and the tournament's A/A and smoke stages; development,
selection, confirmation and replay are frozen but unexecuted in this control.

## Merged implementation and durable Linux replay

[PR 755](https://github.com/aburan28/crypto/pull/755) merged at
`0de6d004d81fb8c2d3e86ecef24c738d86fec853`, following successful applicable
checks on reviewed head `04ca67fca90d2fd6f6c6dc4d4b34ff7ffac469c5`.
The final Linux workflow passed the 100-test harness suite and release arithmetic
controls. Its measured integration census is:

| Control | Verified | Retained failures | Scope |
|---|---:|---:|---|
| Three optimized producer sources | 39 IC + 13 rho native/profile pairs | 0 | Public-point input, certificates, exclusive phase and online timing closure |
| Existing native driver | 17 jobs | 1 intentional one-trial failure | Canonical admission, repeated-run identity and null failed-row costs |
| Existing instruction tournament | 15 native/profile pairs | 0 | Six A/A pairs and nine smoke pairs |

The tournament's development, selection, confirmation and replay inputs were
prepared but those stages were not executed. These controls consume **zero
improvement rounds** and do not qualify an incumbent or support a speedup claim.
Repeated controls across code revisions are not new independent target samples.

The [Linux history](linux-history.json) retains four successful workflow versions,
including the initial producer bundle's omitted transitive dependency and the
older Python 3.12 floating-point summaries. The initial bundle can be replayed
jointly with the exact missing dependency from the same CI run's sealed native
evaluator; its provenance is retained and none of the original artifacts were
rewritten. Later bundles include that dependency. Integer nanosecond aggregation
replaced floating-point process-clock summation, so the final driver artifacts
replay identically under Python 3.11 and 3.12.

[`ic-driver-linux-controls-20260925.tar.zst`](../../evidence/ic-driver-linux-controls-20260925.tar.zst)
is the durable bundle: SHA-256
`972f1b6e16cc9a09960fccc80ba5f64331c3d12a5a312581521c53cea53d334c`,
9,324,450 compressed bytes, 83,113 files and 372,870,950 expanded file bytes.
It preserves exact source, binaries, frozen evaluators, raw jobs, profiles,
certificates, logs, failures, earlier audits and final laptop controls. Build
caches are excluded; original compressed profile bytes are retained unchanged.
All **28 read-only audits** passed after fresh hash-checked extraction; see
[the replay receipt](linux-archive-validation.json). The cross-version checks in
that receipt apply to the final native, failure and tournament driver bundles.

```sh
python3.12 research/ic_candidate_tournament_20260915/evidence/restore.py \
  --archive ic-driver-linux-controls-20260925 --out /tmp/ic-driver-linux-replay
python3.12 research/ic_candidate_tournament_20260915/goal_20260924/driver-admission/audit_linux_archive.py \
  --bundle /tmp/ic-driver-linux-replay/ic-driver-linux-controls \
  --out /tmp/ic-driver-linux-replay/audit.json --also-python311
```

The optional cross-version flag requires `python3.11` on `PATH`. Without it,
25 audits run under Python 3.12. Replay executes the archived independent
checkers, not new performance measurements. The archive remains in Git after
GitHub Actions artifact retention expires.
