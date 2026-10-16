# First retained-N83 mixed-orbit partial-search gate

**Result: `UNKNOWN_search_cap`.** The authorized first trial imported the exact
K=1,182 small signed-Frobenius base and nested K=2,048 envelope, then built a
three-small/two-residual SAT search for the pinned public target. Its one solver
call used the full 5,000-conflict cap with **zero models**. Consequently the
worker emitted zero partial candidates and completed rows. The graph and rank
paths received no natural input. This is a capped relation-stage observation,
not a complete cold index-calculus runtime or factor-base comparison.

| Stage | Descriptive wall | Receipt | Conservative charge |
| --- | ---: | --- | ---: |
| Offline Linux arm64 build | 170.639 s | `PASS_build_only`; binary SHA-256 `264df6474b2acc4865091a3eef8da97424c694affa0dc994971de7d01df4e346` | 600 s |
| S3-backed nested-base import | 8.080 s | `PASS`: 196,212 small points/1,182 orbits; 143,756 residual points/866 orbits | 360 s |
| One mixed-orbit search | 434.368 s | `UNKNOWN_search_cap`: one call, 5,000 conflicts, zero models, zero partials | 1,200 s |
| Cleanup, frozen-input audit and receipts | Included in 444.473 s run-supervisor wall | Both runtime containers exited without OOM and were removed; event replay passed | 240 s |
| **Additional gate** | Build plus run supervisor: 615.111 s | No complete-runtime value or selected base | **2,400 s** |

The full user-approved additional 2,400 seconds are charged under the frozen
stage-cap protocol, even though observed wall was lower. Combined with the
prior 12,600-second conservative ledger, the cumulative local charge is
**15,000 seconds** and this additional grant has **zero seconds remaining**.
The stage walls are L0 diagnostics; there is no matched, fully charged
single-target index-calculus result from which to select the minimum-runtime
factor base.

The search event records 101,636 peak SAT variables, 465,738 domain clauses,
714,072 peak clauses, 4,980 XOR rows and a 468,664,320-byte cgroup memory
peak. These are one capped trial's construction and solver diagnostics. The
solver did not refute the equation. In particular, zero models at a conflict
cap is not evidence that the target has no decomposition in this base.

The source was clean commit `f0a1e55f7cb6d965b0791ef68cb636ec8b8a3a1a`.
The downloaded small and envelope objects matched SHA-256
`eb342b574e244e1e67e8e78c805f33ee88edcdb4f3a14d4114e6f413beede380`
and `aacb62784dee0c9f5e06328dbb85a4b6a1c7243e778328da82900e7284f3a74a`.
Their exact [S3 locations, downloaded sizes and 2026-10-12 head responses](evidence-01/objects.json)
are retained. The heads establish object presence and length; downloaded
SHA-256 values establish the bytes used by this run.
The [authorization](evidence-01/authorization.json), [build receipt](evidence-01/build/build.json),
[run outer receipt](evidence-01/run/outer.json), [budget](evidence-01/run/budget.json),
[worker result](evidence-01/run/run/result.json), [trial event](evidence-01/run/run/events.jsonl),
[compiled worker](evidence-01/build/worker) and raw stage logs are preserved.
The separate [audit script](audit.py) passed and wrote [audit.json](audit.json):
it rechecks committed source bytes, the retained binary, object downloads when
present, stage commands and log hashes, the first-trial coefficient schedule,
search counters, container exits, cleanup, and cumulative accounting.

Validation after the evidence edit passed: `cargo test --release --lib -q`
(4,000 passed, 113 ignored),
`cargo test --release --example koblitz_n83_mixed_orbit_worker -q` (3 passed), the study Python suite
(46 passed), and the boundary autolab Python suite (43 passed). An initial
sandboxed library run was interrupted after local-network and solver-process
tests failed under restricted permissions; the completed host-permission run
is the release-suite result above.

The next meaningful comparison still needs natural relations, verified graph
completion and independent rank gain, sparse linear algebra, one verified
individual logarithm, and a matched rho control with complete cold costs.
This gate grants no further retained-N83 execution.
