# Decision: generic degree-263 full-point addition edge

**Gate: PASS, limited to representation correctness.** The preregistered
[protocol](PROTOCOL.md) was committed before the candidate, and
[FROZEN.json](FROZEN.json) locked source commit
`31997e6b9598b5bdad08d3ca36b6dbfef805bd3a` plus 12 source/input hashes
before the held run. The complete raw [evidence](evidence/final/receipt.json)
records a producer and a separate independent replay, each exiting 0 within
the 600-second child and 1 GiB RSS acceptance limits.

| Panel | Measured result |
| --- | ---: |
| Toy curves over GF(2²), GF(2³), all `a∈{0,1}` and all `b≠0` | 20 |
| Canonical `(P,Q,R)` rows | 14,456 |
| Exhaustive valid and invalid point/λ circuit evaluations | 152,432 |
| Archived raw rows, including invalid encodings and leaf controls | 19,974 |
| Reference/DAG mismatches | 0 |
| Exact degree-263 leaves `[1,0]`, `[1,4]` | 2/2 pass, 40 controls |
| Exact leaf DAG variables / packed-model limbs | 920 / 15 each |
| Exact leaf DAG nodes | 357,429 / 357,438 |
| Producer wall / peak RSS | 6.96 s / 174,931,968 B |
| Independent verifier wall / peak RSS | 7.31 s / 173,473,792 B |

The full-point edge covers canonical infinity, inverses, doubling,
different-x addition, and the x=0 self-inverse point for both actual leaf
coefficients. The separate verifier used polynomial multiplication and long
division, recomputed every archived row, and independently checked both saved
leaf curves. A fresh archive replay returned `ARCHIVE_REPLAY_PASS` with
`rows_sha256=543c8c283d87b36de6701f692ac00bc566a9109e8ee3512ffbf54171a82e370a`.

This gate does **not** establish an implicit m≥3 native PDP producer, a
solver-checkable UNSAT certificate, a leaf relation yield, or a speedup over
rho. In particular, no solver has consumed the 357k-node edge at n=131; the
node count is a construction size, not a solve time. The 16 saved native points
per leaf are test inputs, not a useful-sized factor base. Keep common-unit
full-DLP cost and n131 crossover **unset**.

The next evidence-ranked implementation target is a three-summand native leaf
constraint that binds each point to a preregistered implicit factor domain and
keeps intermediate points fully labelled. Export the resulting exact Boolean
instance and re-check every positive model by an independent point sum. Admit
a solver comparison only after a necessary physical-support/rank screen and
explicit SAT, checked-UNSAT, UNKNOWN, and resource-exit semantics. Then compare
native, transported, pullback, and original bases at equal useful size on held
PDP targets. The current edge is reusable across these arms but cannot decide
which arm wins.
