# V7B same-input repeat after controlled release suite

Run: 2026-10-15 18:22 UTC. Host: Darwin arm64, macOS 26.6,
shared L0. The controlled `cargo test --release --lib -- --test-threads=2`
run had exited with 4,083 passed, zero failed and 117 ignored before this
repeat began. The [frozen V7B protocol](../../REPEAT_V7B.md) kept the same
source SHA-256 `55e23e6c7c361ce295ea7c3313f82ad4734ef900eb7a71bdafc9ac78dbbcea6e`,
binary SHA-256 `d4eb7405211a1a433f570663334c413f8abdce6195a6f43c7c3f31fa087a67ab`,
66 public order-67 targets, H selection, nine-entry source table,
rotated order, caps and three seeds as V7. The solver inputs contained no
answer scalars. The source and binary hashes were checked immediately
before execution.

Each of pairing, direct `F_401` negation rho, and direct source-group
BSGS returned all 66 scalars in all three rounds. A separate arithmetic
checker replayed **594/594** `[k]G=Q` outputs. Common input reading and
subgroup validation are added to every arm at each prefix. Each arm's
CPU includes cold setup, failed work, target solves and final checks.
The BSGS table was rebuilt in every round and reused within that round;
its full-list ledger has nine entries, 283 lookups and zero false
candidates. Its final 66-target operation counts match V7: 488
additions, 552 doublings, 132 scalar multiplications, 3,726 base-field
multiplications, 1,586 squarings and 1,032 inversions. Pairing and rho
primitive ledgers, attempts, wall time, and process peak RSS are retained
in every raw round.

| Arm at 64 distinct targets | Round 0 cold CPU ms | Round 1 | Round 2 | Median |
| --- | ---: | ---: | ---: | ---: |
| Direct source-group BSGS | 0.378 | 0.497 | 0.276 | 0.378 |
| Direct prime-field negation rho | 0.938 | 0.959 | 0.548 | 0.938 |
| Pairing with reused H and finite-field table | 1.219 | 1.837 | 1.014 | 1.219 |

The paired pairing/BSGS ratios at `N=64` were **3.224, 3.696,
3.674** (median 3.674). Direct rho/BSGS ratios were 2.482,
1.930, 1.986 (median 1.986). Pairing/direct rho ratios were 1.299,
1.915, 1.850. Direct BSGS was the lower complete CPU arm at all five
prefixes in all three rounds. The median paired pairing/BSGS ratio was
1.995 at `N=1` and 3.708 at `N=66`; the graph gives the whole prefix
series and observed three-round ranges. V7's concurrent-test `N=64`
pairing/BSGS median was 5.190. The change between V7 and V7B is a
shared-host timing observation, not an inferred effect size from test
content or a portable operation-count comparison.

![V7B complete CPU ratios to source BSGS](GRAPH.svg)

These are whole cold current-thread CPU measurements on this fixed
order-67 subgroup. Host sharing and its residual variability are visible
in the three rounds. A counted `S` conversion, isolated interleaved
timing, second-host replay and increasing subgroup orders remain open.
The canonical ECDLP frontier is unchanged.

[Combined V7B raw](../../COMPACT_BSGS_401_V7B.jsonl), SHA-256
`9208af97ef9e2c0db554c887e324249a2c69f8991f02543d4c7f2e722b26bf70`;
[independent replay](../../COMPACT_BSGS_REPLAY_401_V7B.json), SHA-256
`68aa6cc8e9c3ee48ce1d7135da7326c08ac2345dd94d647124606d3c6376994d`.
Individual round JSONL, stderr and the controlled library-test stdout/stderr
are retained in the study directory. The test stdout is in
`../../../validation_release_lib_controlled_20261015.stdout`.
