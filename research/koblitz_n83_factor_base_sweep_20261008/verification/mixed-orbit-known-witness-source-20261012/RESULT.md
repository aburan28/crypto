# Archived N83 known-answer control, source preparation only

The frozen [fixture](fixture.json) gives an exact 3-small/2-residual relation
in the previously imported K=1,182/K=2,048 signed-Frobenius pair. It selects
the phase-zero positive point from small columns 0, 1, 2 and residual columns
1182, 1183, then sets `a=0`, `b=1` and `T` to their group sum:

`T = (6595116073687692399614168, 9473620351503970875544348)`.

[known_witness.py](known_witness.py) reopens the exact local copies of the two
S3 objects, checks their compressed SHA-256 values against the prior download
receipt, compares all 196,212 nested-prefix rows, checks the complete archive
row counts, verifies the selected points and `T` on the curve and in the
primary subgroup, and recomputes the group sum. The prior retained-search
audit is bound by SHA-256
`32420a62bb6ccf8f9df7f03dea9040e518606740968110dc8c0ac47cb50cafd5`.
This construction has a guaranteed relation by design. It measures neither
natural relation yield nor success on the pinned public target.

The Rust `koblitz_n83_mixed_orbit_worker` now has `verify-control` and
`control-run` modes. `verify-control` independently imports the exact archive
pair, compares the passed fixture to compile-time pinned bytes, checks every
selected index against the imported Rust bases, and recomputes the group sum.
This source-only validation **passed** for the retained objects with
`solver_search_executed=false`. The new `control-run` mode pins `a=0,b=1`,
uses the same 3+2 domains and first-pair ordering as the natural worker, and
requires a 4 GiB zero-swap Linux cgroup. It writes separate control-only
configuration and result receipts, including solver stats and a group replay
of any returned candidate. It is prepared but **has not been executed**.

[control_supervisor.py](control_supervisor.py) now makes the proposed next gate
concrete. Its authorization schema is `n83.additional-local-budget/v3` with
`prior_total_charged_seconds="15000"`,
`purpose="one_mixed_orbit_known_witness_control"`, exact fresh build/run paths,
and a new explicit user-approval reference. It rejects the expired v2 schema
and any allowance other than 1,500 seconds. An offline, pinned-image Linux
build has a 300-second conservative charge; the read-only import/control
verification has 120 seconds; one 5,000-conflict/32-model SAT call has 900
seconds; cleanup and audit reserve 180 seconds. Launched stages charge their
whole caps even if they fail early. Docker stages use no network, the exact
S3-object downloads, and verified container removal. The search runs at one
CPU, 4 GiB, zero swap and a worker wall cap of 840 seconds.

The supervisor freezes the clean source commit and SHA-256 values for code,
fixture, manifests, objects, build receipt and binary. It validates the prior
15,000-second audit hash before launch; rechecks the source and inputs after
execution; retains raw stdout/stderr, container state, cleanup, an outer
disposition and the conservative budget ledger; and independently replays any
successful candidate against the archived point domains and group equation.
Its [four source-only tests](test_control_supervisor.py) verify stale and
underfunded grants are rejected before launch, a bare success label cannot
pass without a candidate, a cap cannot be called a success, and the worker wall
cap needs a matching receipt. Preparing and testing this supervisor has not
launched a container or SAT search.

The [constructed replay check](replay-check.json) opens the same compressed
archive bytes and passes the five-point witness through the supervisor's
candidate validator. It then changes a small-point index, replaces a residual
point with a small-base point, and changes the fixture digest; each mutation
is rejected. [replay_check.py](replay_check.py) reproduces this record. These
are deliberately constructed receipts for auditing the validator, not output
from a solver. The supervisor now also requires the worker's fixture and both
object BLAKE3 values to match the pinned identities.

The next controlled experiment is exactly one `control-run` call under a
separate authorized local budget and this external container supervisor, with the
two archive hashes and fixture hash frozen before launch. A bounded failure
would leave solver access unknown; a valid candidate would establish only
known-answer solver access. A completed refutation would be a producer
failure, since the fixture contains a checked witness. Neither outcome can
rank factor bases by total index-calculus runtime. Natural relation yield,
graph closure, rank, sparse linear algebra, individual-log work and matched
rho remain unmeasured for this pair. The previous 15,000-second cumulative
charge is unchanged by this source preparation and no additional search
authorization is implied.
