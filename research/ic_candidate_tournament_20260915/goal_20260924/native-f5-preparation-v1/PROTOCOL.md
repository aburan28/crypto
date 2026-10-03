# Native retained F5 preparation replay

This gate replaces the historical Python preparation verifier with Rust for
the disclosed synthetic fixture. It starts no solver, new relation query or
target solve. Historical F5 and SAT registrations and all three tournament
confirmation sets stay closed. The next complete native F5 controller needs
its own before-execution source freeze and bounded registration.

The hypothesis is that independent native reconstruction reproduces all
ordinary inputs and the complete matrix/log state of the accepted F5
certificate without modifying its bytes, vocabulary or provenance. The
reference is the retained certificate, not a timing or matched-rho arm:

- Whole canonical certificate: `d8c5d6679fe89561606154785fdfba763835bf261dceb08fae193eff9f5636ae`.
- Ordinary input seal: `6678815b78b9414d81f94b595d3e10ab38f388511fa3a6831891f685e96f833d`.
- Mathematical state: `edbff76da6442b9f2e5e8235682c9bf1052465f310c765ba8a37d018a60bf107`.
- Exact curve representation: `EC1N17Ckb1hbbe2b5b6b1e6`, field modulus 131081,
  subgroup order 65587, cofactor 2, generator `[43693,23339]`.

Success requires all 216 original ordinary queries in order, 61 group-verified
witness rows, 155 independently proved geometric negatives, full rank 29,
32 dependent rows, no duplicate relations and all 29 independently recovered
and scalar-verified logs. Keep the 63 geometric points, including the torsion
point: their subgroup images contain 62 usable points and 29 folded columns.
The original outcomes remain `witness` and `proved_unsat`; do not normalize
the original certificate to the SAT vocabulary.

Replay must reject an altered external seal, target contamination, missing or
reordered geometry, false witnesses or negatives, a failed query entering the
matrix, changed chronology, projection, rank trajectory, column logs,
provenance or promotion flags. Stop on any discrepancy. Do not repair or
regenerate the historical certificate, dispatch a solver, or reinterpret
rank-stopped query counts as fixed-sample yield.

The native reader uses bounded, no-follow input reads and create-only output.
It checks input and checker executable bytes before and after reconstruction.
The receipt binds checker, independent curve oracle and orbit helper source
digests, binary digest and host architecture/OS. The existing SAT mathematics
module supplies the orbit and exhaustive pair-sum helpers without changing its
source, because that exact module is part of the closed SAT frozen audit.
F5 matrix elimination and retained outcome checks remain explicit in the new
F5 module. The production Rust solver is not used for this independent replay.

All replay costs are offline verification. No online arithmetic or calibrated
timing is measured; online time, operation calibration and speedup stay null.
Fresh/headline/promotion eligibility and complete native scientific runtime
admission remain false. Original Python preparation provenance is retained.
Passing this gate supplies neither new native natural yield nor a fresh
tournament comparison.

```sh
rustc --edition 2021 src/bin/isolated_bench/busy_unix.rs -o /tmp/native-f5-busy
/tmp/native-f5-busy busy -- cargo build --locked --bin icprog
/tmp/native-f5-busy busy -- target/debug/icprog f5-preparation-replay \
  --root . --out /tmp/native-f5-preparation-new-replay.json
```

Linux x86-64 and macOS ARM64 CI independently rebuild and replay these data,
including rejection controls. Never execute an archived Mac binary on Linux.
Replays are data/math checks, not hardware speedup tests. Keep native F5
controller/preparation production, new ordinary-yield evidence, exposure
reconciliation and fresh paired ecbench qualification marked pending.
