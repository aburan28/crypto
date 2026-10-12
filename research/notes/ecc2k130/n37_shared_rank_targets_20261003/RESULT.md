# n37 shared-rank point-only target recovery

Decision: **`SHARED_RANK_TARGET16_RECOVERY_ADMITTED_N37_SOURCE`** as a
correctness and cost-accounting gate. The [protocol](PROTOCOL.md) was
committed as `611cf9db` and [PR #1309](https://github.com/aburan28/crypto/pull/1309)
opened before target generation. The point-only input, separate scalar
fixture and independent input receipt were committed as `5a8f7fff` before
the target producer ran. The producer and replay source were committed as
`c4219948`; the verifier's subsequent correction to include mandatory
native lookup units in its GAE audit is preserved in this PR.

The native fixture preparer inventoried 126 tracked historical point files
with 55,368 rows and excluded 13,443 n37 signed-Frobenius orbits after
adding the 3,108 frozen base points and 55 rank probes. A separate replay
regenerated all candidate scalars with the general binary-curve group law,
checked the inventory and every `[d]G=Q`, and verified that the 16 accepted
point-only targets were orbit-disjoint from this frozen set. No candidate
was rejected in this deterministic selection. The target producer opened
only [`TARGETS.points.jsonl`](TARGETS.points.jsonl), never the
[`FIXTURE.json`](FIXTURE.json) containing the scalars.

All **16/16** public points were recovered on the direct three-summand
attempt, with **16/16** independent full-point scalar replays. No shifted
residual or failed PDP attempt occurred, so this run does **not** empirically
qualify the 1–63 shifted-attempt fallback. The independent target verifier
used the archived point-defined source base and general binary-curve group
law to check each witness sum, folded row, cofactor adjustment, recovered
scalar and `[d]G=Q`. It also compared the new rank setup against the
independently verified #1302 transcript after removing only timing fields.
One altered witness and one altered recovered log each produced a `FAIL`
receipt. A second producer run agreed with the committed run on every
non-timing field.

| Exclusive charged phase | GAE lower bound | Notes |
| --- | ---: | --- |
| Shared base, table, rank and checks | 194,826 | Built once without Q; rank 42/42 |
| Target query and subgroup checks | 736 | 16 point-only Q; target-dependent |
| Target PDP | 1,220 | 602 lookups and 618 additions |
| Target witness checks | 48 | All 16 three-point sums checked |
| Target log combination | 0 | Modular arithmetic remains unpriced |
| Target scalar replays | 624 | 16 full-point checks |
| **All 16 targets** | **2,628** | Sum of exclusive online phases |
| **Shared setup plus 16 targets** | **197,454** | Incomplete native-work price |

The target online intervals sum to 1.266917 ms on an **unisolated** macOS
arm64 host. This wall observation is a diagnostic only; it is not a
one-target online speed result or a multi-target speedup. Native field
operations, hash work, allocations, modular log combination and peak table
memory do not have a common calibrated price in this record. The GAE total
is a lower bound, not a fully charged `S`. There is no same-Q strong rho
arm, denominator, confidence interval, crossover or n131 extrapolation.
Class: **accounting / correctness diagnostic**. No algorithmic boundary
ratio is measured here.

The [evidence manifest](EVIDENCE.json) pins source, binaries, lockfile,
inputs, [raw transcript](RAW.json), [independent target receipt](REPLAY.json)
and both mutation results. The next measurement must first pair one
previously unseen point with a fully verified strong signed-Frobenius rho
online solve under one resource envelope and include all target-dependent
work. A separately frozen 1,024-point same-Q batch can then test shared
setup amortization; it must charge the setup once, all descents, native
work and memory, and use an isolated A/A-controlled timing protocol. Neither
comparison is supplied by this correctness gate.

Reproduce the frozen input and transcript replay from the repository root:

```sh
cp research/notes/ecc2k130/n37_native_m6_mitm_20261002/Cargo.lock Cargo.lock
cargo build --release --locked --example n37_shared_rank_targets --example n37_shared_rank_target_replay
evidence=research/notes/ecc2k130/n37_shared_rank_targets_20261003
target/release/examples/n37_shared_rank_target_replay "$evidence/TARGETS.points.jsonl" "$evidence/FIXTURE.json" /tmp/n37-shared-target-input-replay.json
cmp /tmp/n37-shared-target-input-replay.json "$evidence/INPUT_REPLAY.json"
target/release/examples/n37_shared_rank_target_replay "$evidence/TARGETS.points.jsonl" "$evidence/FIXTURE.json" "$evidence/RAW.json" /tmp/n37-shared-target-replay.json
cmp /tmp/n37-shared-target-replay.json "$evidence/REPLAY.json"
```
