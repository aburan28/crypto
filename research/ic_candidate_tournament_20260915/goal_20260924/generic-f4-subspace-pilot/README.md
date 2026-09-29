# Disclosed-point F4/F5 standard-subspace recovery pilot

Status: **dispatch diagnostic measured; recovery under max_trials=1 failed; not a qualification**.
`promotion_eligible=false`. Canonical receipts are [RESULT.md](RESULT.md).
[`run-20260929`](run-20260929/RESULT.md) is a same-budget replay and is not a new registration.

See [RESULT.md](RESULT.md). Encoder-feasible `standard_subspace` d=6 lets pinned
`f4`/`f5` enter MatrixF4/MatrixF5 with `unsupported: false` on all five
disclosed inventory points. Zero relations/solutions at `max_trials=1`.

| Gate | State |
| --- | --- |
| Feasibility gate + d6 inventory control | Merged via PR [#952](https://github.com/aburan28/crypto/pull/952) |
| Live v2 measurement | Do not retry seed `2026092902` |
| This pilot (dispatch series) | Complete — [RESULT.md](RESULT.md) (PR #961) |
| Same-budget replay | Retained in [`run-20260929`](run-20260929/RESULT.md); do not rerun |
| Fresh competitive F4/F5 registration | Blocked until recovery mechanism + v2 exposure census |

No scoreboard competitive row is updated from this directory. The retained
run is a disclosed-point diagnostic with unknown complete cost.
