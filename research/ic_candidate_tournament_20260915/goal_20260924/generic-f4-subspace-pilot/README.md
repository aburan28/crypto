# Disclosed-point F4/F5 standard-subspace recovery pilot

Status: **one local execution retained** in [`run-20260929`](run-20260929/RESULT.md). `promotion_eligible=false`.

The ambient-orbit F4/F5 arms in
[generic-backend-qualification-v2](../generic-backend-qualification-v2/README.md)
are statically encoder-infeasible (`4n > 64` variables) per the post-registration
audit in that directory's `STATIC-FEASIBILITY` files (PR #952). SAT arms are
out of scope for that static check. The live v2 campaign
([run 36580669479](https://github.com/aburan28/crypto/actions/runs/36580669479))
remains the sole empirical record for seed `2026092902` and must not be
redispatched.

This pilot answers the next bounded question on **already disclosed** points:
does `standard_subspace` dimension 6 let the pinned worker actually dispatch
F4/F5 and, under a declared budget, recover any verified one-target logs?
See [PROTOCOL.md](PROTOCOL.md).

| Gate | State |
| --- | --- |
| Feasibility gate + d6 inventory control | Merged in PR #952; this run's preflight passed |
| Live v2 measurement | Sole dispatch 36580669479; do not retry |
| This pilot execution | Done once: 10/10 jobs entered MatrixF4/MatrixF5, node budget exhausted, zero relations |
| Fresh competitive F4/F5 registration | Not authorized by this pilot |

No scoreboard competitive row is updated from this directory. The retained
run is a disclosed-point diagnostic with unknown complete cost.
