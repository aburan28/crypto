# Native F5 target mathematics and bounded job transport

This is a verification foundation for the pending native F5 controller. It
replays retained data from the consumed prepared F5 v3 control and executes no
decomposition solver, ordinary query, target generation or archived executable.
It is not an executable registration. The original control stays closed.

## Fixed question and inputs

Can independent native code reconstruct every target query, distinguish a
complete geometric negative from an exhausted solver, and verify recovery and
the exclusive phase ledger without trusting the producer's scalar?

The only mathematical instance is the disclosed synthetic n17 control:
`EC1N17Ckb1hbbe2b5b6b1e6`, field modulus 131081, curve a=b=1,
subgroup order 65587, cofactor 2, G=(43693,23339), and known public target
Q=(52411,72106). The retained query seed is 2026093032 and the cap is eight.
These are historical replay inputs, not seeds for another dispatch.

The whole preparation certificate is pinned to canonical SHA-256
`d8c5d6679fe89561606154785fdfba763835bf261dceb08fae193eff9f5636ae`.
The raw target report `prepared-f5-v3-control-v1/controls/pipeline.stdout`
is pinned to SHA-256
`57cb2e8e885c5da4912abd32b767b26ac6263175c5ec98f8f20b93fbc8d250ca`.
Neither historical file is modified. The preparation verifier independently
reconstructs the ordinary matrix and logs before target verification.

## Independent checks

The verifier reconstructs the seeded query law and each aG+bQ. It proves a
claimed three-summand negative by exhaustive finite-group pair-complement
checking over all 63 geometric points, including killed torsion, repetitions
and identity pair sums. A positive witness must group-readd. Recovery uses
cofactor projection and the 29 independently verified column logs, followed
by an independent scalar multiplication to Q. The actual usable base count
is 62; it is never replaced with the 63-point geometric count or 29 columns.

The producer's bounded Macaulay MatrixF5 engine, degree-three limit and 8192
reduction budget are checked. This engine uses an F5 row criterion; this does
not claim implementation of a complete incremental F5 Gröbner basis algorithm.
An exhausted or unsupported attempt cannot become a proved negative. Attempts
must appear in query order, stay within the cap and stop after recovery.

All five exclusive target phase slots must be present. A completed target
requires all five integer costs and exact closure to the recorded online
interval. An incomplete target retains null unperformed phases and cannot
claim a scalar or complete online cost. The replay does not newly measure
these historical costs, establish a calibrated host, or add an external
verification interval to the old timer.

## Native launch foundation

The shared watchdog accepts an explicit bounded stdin job (at most 1 MiB) and
declared environment. It clears ambient environment variables, rejects
duplicate keys and a conflicting helper handshake, hashes the executable
before and after, preserves raw output and creates a durable launch receipt.
Stdin writes are nonblocking and advance inside the deadline loop; a worker
that never reads the job cannot block timeout enforcement. Receipts retain
the requested input hash, delivered byte count and any delivery failure.
The watchdog terminates the process group and confirms drain. A delivery
failure cannot admit a successful launch merely because the worker exits zero.

The transport tests use shell/cat/sleep correctness fixtures, not cryptographic
producer calls. They compare an exact job larger than pipe capacity, verify
the declared environment, and require a nonreading worker to time out and
retain its partial-input receipt. Historical handshake/watchdog tests remain.

## Commands and evidence

Run builds and tests under the repository's native busy lock. This lock
serializes work; it does not establish a quiet or calibrated measurement host.

```sh
/tmp/ic-native-busy busy -- cargo test --locked --bin icprog f5_target -- --test-threads=1
/tmp/ic-native-busy busy -- cargo test --locked --bin icprog sat_control::native -- --test-threads=1
/tmp/ic-native-busy busy -- cargo build --locked --bin icprog
/tmp/ic-native-busy busy -- target/debug/icprog f5-target-replay --root . --out /tmp/NEW-f5-target-replay.json
```

The output is create-only. CI on Linux and macOS tests the verifier and shared
watchdog and uploads each platform's native replay receipt. No archived macOS
binary executes on Linux. The receipt retains checker source/binary identity,
the old native producer build identity and historical Python orchestration
provenance. This is postexecution analysis and does not replace the original
preexecution-frozen auditor or retroactively repair its source coverage.
Checker self-identification allows a bounded 512 MiB executable so Linux debug
test builds can be identified. It retains the regular-file/symlink gate and
reports the path, file size and limit on read failure. This is replay provenance,
not a relaxation of the scientific worker's native launch limits.

## Gates still open

The native F5 controller still needs a complete source/dependency/build freeze,
its mathematical-only stdin fixture, exact native environment, one-use claim,
separately published registration and preexecution-frozen independent auditor.
A new disclosed control must use its own frozen seed and stop conditions;
never retry an old registration or choose a seed after searching for success.
This PR supplies no such registration and launches no scientific worker.

New ordinary yield for both families and a new frozen paired comparison remain
pending. That comparison needs the exposure census, complete IC1 method and
workload identities, native ecbench admission, matched rho/incumbent,
calibration, resources, arm order and uncertainty analysis. The primary
question is one previously unseen target's independently verified online wall
time. No new candidate, workload or run result is minted by this replay.
Source-bound runtime admission, fresh qualification, headline eligibility and
promotion remain false; reference, operation calibration and speedup are unknown.
The three consumed confirmation sets remain closed.
