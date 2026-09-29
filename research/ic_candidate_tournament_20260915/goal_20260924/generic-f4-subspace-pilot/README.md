# Disclosed-point F4/F5 standard-subspace recovery pilot

Status: **dispatch diagnostic measured; recovery under max_trials=1 failed; not a qualification**.

See [RESULT.md](RESULT.md). Encoder-feasible `standard_subspace` d=6 lets pinned
`f4`/`f5` enter MatrixF4/MatrixF5 with `unsupported: false` on all five
disclosed inventory points. Zero relations/solutions at `max_trials=1`.

| Gate | State |
| --- | --- |
| Feasibility gate + d6 inventory control | Merged via PR [#952](https://github.com/aburan28/crypto/pull/952) |
| Live v2 measurement | Do not retry seed `2026092902` |
| This pilot (dispatch series) | Complete — [RESULT.md](RESULT.md) |
| Fresh competitive F4/F5 registration | Blocked until recovery mechanism + v2 exposure census |

No scoreboard competitive row is updated from this directory until a later
evidence PR lands measured, verified receipts.
