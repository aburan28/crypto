# n37 target-blind shared-rank gate

Decision: **`TARGET_BLIND_RANK42_ADMITTED_N37_SOURCE`** for the exact frozen
42-column source base and counted signed-Frobenius three-summand oracle. The
[protocol](PROTOCOL.md) was committed as
`5d13c68ea7d9ac89df91d44e2f418a32ccdcc7f7` and [PR #1302](https://github.com/aburan28/crypto/pull/1302)
opened before measurement. On the registered `icv1-f2m37-tm534059-32aad96b`
curve, **55 target-blind `[a]G` probes produced 55 valid relations, 42
independent rows and 13 dependent rows**. Rank reached **42/42** well before
the 1,000,000-trial cap. No public target Q or target scalar was supplied to
this rank setup.

The separate Rust replay reconstructed the ordered 3,108-point source base
from the frozen `SOURCE42.jsonl`, recomputed all 55 probes through the general
binary-curve group law, verified every witness sum and dense row, and solved
every accepted-row prefix with its own Gauss-Jordan code. It checked all
**42 solved base-column logarithms** by full-point multiplication, including
the orbit coefficient and subgroup cofactor. Both final raw runs report
`verified=true`, `exhausted=false`; both independent receipts say `PASS`.
Two deliberate alterations, a witness index and a table-addition count,
produced `FAIL` receipts. The Linux workflow replays the archived runs,
regenerates non-timing evidence, and repeats those rejection controls.

| Exclusive target-blind setup phase | Counted GAE | Important native or group work |
| --- | ---: | --- |
| Source base | 4,260 | 3,108 signed points, 42 columns; 105 abscissae scanned |
| Folded pair table | 178,710 | 66,822 additions; 64,467 retained entries; 111,888 Frobenius maps and lookups |
| Rank search, including all trials | 8,498 | 55 scalar probes, 55 oracle hits; 3,164 canonicalisations |
| Dense elimination | 1,110 | 55 rows, 42 independent, 13 dependent; 1,110 row operations |
| Witness and base-log verification | 2,248 | 55 witnesses and 42 full-point column checks |
| **Shared setup total** | **194,826** | No target descent or rho run in this gate |

The five phase charges sum exactly to the reported total using
`Calibration::default()`. This is a **charged-operation lower bound**, not a
fully priced native cost: field operations, hash work, allocations and full
memory use remain unpriced. The table is 91.7% of the charged setup. The
base record gives a 149,184-byte factor-base vector floor; it excludes the
table's hash allocation and allocator overhead. The source builder's
`wall_ns` was initially left at zero even though its operation ledger was
complete. We preserved those first two raw runs as `RAW_INITIAL*.json`,
committed the timer and replay-ledger correction as
`14aef9e7de84e353f9a59153498c8dca712ecd31`, and reran exactly the same
frozen workload. The four runs match on every non-timing field. Final A/B
base times are 7.576/4.664 ms and table times 85.358/61.284 ms on this
unisolated macOS arm64 host; they are diagnostics, not admitted speed data.

The exact [evidence manifest](EVIDENCE.json) records raw, replay, binary,
lockfile, source and normalized-record SHA-256 hashes. The two final raw
records are [A](RAW.json) and [B](RAW_B.json), with independent receipts
[A](REPLAY.json) and [B](REPLAY_B.json). The earlier timer-incomplete runs
and receipts remain [here](RAW_INITIAL.json), [here](RAW_INITIAL_B.json),
[here](REPLAY_INITIAL.json) and [here](REPLAY_INITIAL_B.json). Final A, final
B and the initial records share SHA-256
`aa650e78c32a66e0496dc8fd1c0fbda31bf8d5ce8b96e139a78d1d520475f321`
after deleting only the five `wall_ns` fields and applying `jq -S`.
The ordered point-and-label digest is
`8460ac4c28515db701c3897a03b4ce0f28abf7cd56fd759ad095f98436a76dcf`;
the source file digest is
`0a32de24a5680ff46baf9543e8bf8e32447323491e2a2adedce2548c14a25f75`.

This closes the **reusable-rank correctness gate only**. The earlier n37
single-target rank-43 result carried a target unknown and cost 22.6 times
its matched strong rho in the bounded cold charged metric. This 42-column
system can instead be built before any Q is supplied, but it has not yet
recovered a fresh Q with that precomputation. There is **no complete IC
cost, matched batch-rho denominator, speedup or crossover** here. The result
does not transfer to n41/n53, to a degree-263 descendant, or to ECC2K-130
at n131 by itself, and does not alter a boundary exponent.

The next registered comparison should freeze a new point-only Q cohort and
a residual policy before measurement, charge this rank setup once, charge
every target descent and verification, and solve the same cohort with a
strong signed-Frobenius batched rho under matched resource and memory limits.
That is the first test of whether shared rank improves the complete batch
economics; an online single-target claim still needs its own paired online
measurement.

Reproduce from the final source revision with the archived dependency lock:

```sh
cp research/notes/ecc2k130/n37_native_m6_mitm_20261002/Cargo.lock Cargo.lock
cargo build --release --locked --example n37_shared_rank_folded --example n37_shared_rank_replay
target/release/examples/n37_shared_rank_replay research/notes/ecc2k130/shared_rank_folded_table_20261003/RAW.json /tmp/n37-shared-rank-replay.json
cmp /tmp/n37-shared-rank-replay.json research/notes/ecc2k130/shared_rank_folded_table_20261003/REPLAY.json
target/release/examples/n37_shared_rank_folded /tmp/n37-shared-rank-new.json
target/release/examples/n37_shared_rank_replay /tmp/n37-shared-rank-new.json /tmp/n37-shared-rank-new-replay.json
```
