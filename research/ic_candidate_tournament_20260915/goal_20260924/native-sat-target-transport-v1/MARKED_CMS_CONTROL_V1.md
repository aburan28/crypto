# Marked CryptoMiniSat stdin control, version 1

Status: frozen *design* for a future disclosed-input compatibility run. No
marked CMS binary, result, SAT target registration or speed claim exists yet.
Do not launch these roles while the separately frozen natural SAT preparation
is running. This control cannot change, refill, or reopen that preparation or
the three closed confirmation sets.

## Inputs fixed before building the marked binary

Use the three already disclosed n17 source instances in
[CONTROL_PLAN.md](CONTROL_PLAN.md), in order `query-00`, `query-01`,
`query-02`. The exact point and manifest/ANF/CNF/Magma SHA-256 values in that
table are required inputs, not outputs to rediscover. The accepted CMS binary
must hash to
`6c509f09622f103d8a3ad90afc151e1c4031275c052d7f8af481465b4f27f2af`.
The marked variant is built from the retained CMS source commit
`7ae1b4a74259cdce223a584281fb8f090bbd3eed` and the exact
[post-buffer patch](cms-stdin-ready-postbuffer.patch). The earlier
[one-line sketch](cms-stdin-ready.patch) is not an eligible primary-online
readiness implementation: `StreamBuffer` allocation would still occur after
its marker. The accepted source archive has
SHA-256 `467b1c3d00a7d6e893332b4d8b42c6326301974d22885aa745f1f926da050323`.
Retain that archive, patch, toolchain, dependencies, build command/log,
binary and all hashes before the first solver role. It is a **new** candidate binary, never the
accepted CMS executable under a changed label. Review that the marker is
flushed after stdin stream, DIMACS parser and stream-buffer construction,
before the first CNF read. If source review or build provenance fails, stop with that failure
and do not launch a control solver.

For each source instance run these three roles once, in this order:

1. Accepted CMS, file argument naming the frozen CNF.
2. Accepted CMS, launched without a file, observed at its original stdin
   reader line before feeding the same CNF.
3. Marked CMS, launched without a file, observed at the exact new
   `c PREPARED_STDIN_READY_v1` line before feeding the same CNF.

Use one CMS thread, `--verb 1 --random 1 --maxsol 1 --maxconfl 1000000`,
and a **120,000 ms wall watchdog per role**, fixed for every point and mode.
The old 60-second `query-00` file timeout remains a failed historical control;
this 120-second panel is separate and cannot retroactively pass it. Use a new
create-only output tree. Retain raw stdout/stderr, exact stdin bytes and
write-count/hash, startup/READY timestamps, PID and group-drain receipts,
exit code, timeout flag and status for every role, including failures.
Perform all solver work after the ongoing natural panel releases the host.
The checked [build recipe](build-marked-cms.sh) also waits for that panel's
original terminal and a passing, full-rank original audit. It creates a new
output directory, extracts only the accepted archived source/dependencies,
applies the pinned post-buffer patch, and retains compiler output and binary
hashes. A successful build receipt says `BUILT_UNVALIDATED`; only the nine
disclosed solver roles and a separate data-only replay can pass transport
parity. The build recipe is never a target execution or speed measurement.
The source-committed `ic_marked_cms_control` Cargo example runs the three
roles in the fixed order for each disclosed point. It uses create-once output,
fixed source hashes and expected classes, bounded process groups, exact stdin
READY markers, and retains every role's raw output and watchdog/drain receipt.
It can report provisional parity only; its result explicitly requires a
separate data-only replay before transport admission.
The `icprog marked-cms-transport-audit` command is that data-only replay. It
checks the original nine role receipts, stdout/stderr and PID ledgers, exact
READY-before-stdin records, disclosed source hashes, native exit/status,
source-valid SAT models, and full-point geometry. It runs no solver. A pass
admits only the disclosed transport parity; it still sets target execution,
natural yield and online speedup to unadmitted or null.

The run passes transport parity only if all nine roles finish within their
fixed watchdogs, all three roles for each point return the same recognized
status, the marked READY occurs before any CNF byte is written, and all
process groups drain. On `SAT_MODEL`, parse each original model separately,
verify it against the frozen ANF and CNF, and independently lift it to a
full-point witness on the declared query. Models need not have identical
bytes. On `SOURCE_UNSAT`, independently enumerate the geometric three-sum
class and require it to be empty. An incomplete, timeout or native error row
is retained and makes this parity gate fail; no rerun or increased limit is
allowed under this version. Verify source files and binary hashes again after
each role. The checker must run data-only against the original execution and
must not call CMS.

The separate prepared-exporter controls already passed on these three points
for their earlier source snapshot. A target capsule must repeat or replay
that parity with the **exact exporter build** it pins. Even complete transport
parity is only a prerequisite: the full natural SAT preparation, one-use
target capsule, fresh card, original execution audit, same-point references
and isolated-host receipt remain separate gates.
