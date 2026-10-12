# Preregistered full-rank continuation for the n37 compact-orbit pipeline

Status: **preregistered, unmeasured**. Push this protocol and open its PR
before running the frozen spec. This follows the merged [cold-source
measurement](../cold_compact_orbit_ic_20261003/RESULT.md), whose one-target
scalars verified but whose relation matrices stopped at rank 19–38 of 43.

## Boundary and hypothesis

The registered public-challenge ladder rung is `E_0: y²+xy=x³+1` over
`GF(2^37)`, ICV1 `icv1-f2m37-tm534059-32aad96b`, prime subgroup order
`r=230603167`, signed-Frobenius order `A=74`. Use
`S=total charged group-addition equivalents / sqrt(r)`; the generic-group
floor is `sqrt(pi/(2A))=0.1456948091`. The same-Q reference is
`rho.signed_frobenius_strong`. The fixed 42-column base and eager folded
three-summand table require 66,822 group additions before any query;
their addition-only setup floor is `S=4.4003461`.

Hypothesis H1: continuing the *same deterministic relation stream* after
the target is pinned reaches rank **43 of 43** and independently verifies
the target scalar on every new Q within the fixed cap. The control stops at
the first pinned target; its scalar and its entire relation prefix must
agree with the full-rank arm. H2: full-rank continuation costs additional
online work, but the eager-table cold single-target lower bound remains
above the measured strong-rho cost on these Q. These are finite n37
questions, not an ECC2K-130 or asymptotic crossover claim.

## Frozen method and correctness gates

Use the already audited native
`compact-orbit-scan:columns=42,raw_x_cap=1000000` factor base and
`mitm-frobenius-counted:m=3` oracle. Their source inventory must retain
the 42 columns, 3,108 signed points, source SHA-256
`0a32de24a5680ff46baf9543e8bf8e32447323491e2a2adedce2548c14a25f75`
and base BLAKE3
`8423b135df3515b0e284126d9eeb95f71e3901ae12908b3cad8eb03d1f31d4bc`.
Keep `Cargo.lock` SHA-256
`b28d3c2d81146a40d00df85f45c6460d9e2bc5c26875307a06cfd15efff99365`.

Implement `incremental-gauss-full-rank` as a new `linalg` choice. It uses
the identical dense incremental matrix and relation stream as
`incremental-gauss`, but continues until all 43 columns are independent.
Only then read the target scalar, charge `[d]G=Q`, and report success.
Record the first rank/trial at which the target was pinned, the final
rank, independent and dependent relation counts, and exhaustion. If the
cap is reached after early pinning but before rank 43, report an exhausted
full-rank attempt, never a successful full-rank solve. Reject an
inconsistent row or wrong scalar. The old `ic.pipeline` default and its
historical replay/serialization must remain unchanged.

The two IC arms use exactly the same `factor_base`, `oracle`, `targets=walk`,
`max_trials=1000000` and algorithm seed per paired run; they differ only
in `linalg`. Before measurement, a deterministic unit test must show the
early-pin/full-rank distinction and cap behavior. On measured Q, the
first-pin checkpoint and relation prefix must agree; if the session format
does not expose individual relations, compare prefix counters and scalar,
and label that weaker check accurately.

## Frozen new public-Q measurement

Use [SPEC.json](SPEC.json) with `target_kind=public`,
`targets_per_curve=8`, **new** `target_seed=202610030741` (disjoint from
the earlier `202610030641` set), alternating three-arm order, one warmup,
two measured rounds, session seed `202610030742`, L0 operation-count
isolation, and 600-second per-run timeout. `ecbench plan` must show all
workloads registered; save its exact output before running. This is 72
executions including warmups. The solver never receives planted scalars.
No support, seed, cap, unit, arm or admission rule may change after any
target outcome is seen.

Run `verify --replay-all`. H1 passes only if every one of 16 measured
full-rank runs reaches rank 43, returns a scalar independently checked by
`[d]G=Q`, the 16 early-pin controls and 16 strong-rho runs also verify,
and all 48 measured replay records are identical. Report every failure,
timeout, rank and target individually, with any run that fails a gate
left in the session. H2's single-target cost rejection passes only if the
minimum full-rank IC **charged lower bound** exceeds the maximum strong-rho
charged value over the same 16 Q-round pairs. Save exact per-Q paired
costs, phase splits, the generic-floor and strong-rho ratios, source and
input hashes, binary and host identity, and the audit receipt in this PR.
L0 wall times are descriptive. Native field and hash work still unpriced
must remain explicit; neither arm earns a full-speed claim from charged
lower bounds.

If H1 fails, the decision is a failed full-rank gate at this cap, with the
observed rank distribution and bottleneck, not an invitation to raise the
cap on these Q. A new support, arity or cap would require a new protocol
and disjoint Q. A passing n37 gate still does not imply a useful batch
amortisation: shared IC setup must face an actual automorphism-aware
**batched** rho reference, full rank, memory and all native work. n41/n53
counting capacity and degree-263 descendant transport remain separate
follow-up gates; no n131 transfer follows from this measurement.
