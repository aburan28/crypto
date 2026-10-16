# K=600 mixed-orbit natural-search gate, source preparation

`k600_mixed_orbit_supervisor.py` prepares an exact, opt-in natural-target
comparison using the already stored K=600 small base and K=2,048 envelope.
It leaves the historical K=1,182 supervisor and its receipts unchanged. The
new gate invokes the worker's `verify-alt` and `run-alt` modes under the same
public corpus, first trial coefficients, three-small/two-residual shape,
first-pair ordering, 5,000-conflict cap, 32-model cap, one CPU, 4 GiB and zero
swap as the earlier K=1,182 retained search. It independently replays any
reported partial or completed relation against the archived point domains,
source equations, graph coefficients and group law. A cap, missing receipt,
resource exit or inconsistent replay cannot be promoted to relation yield.
Each event must also report the source-derived 415,686 domain clauses for the
K=600 three-small/two-residual shape, so archive or domain-selection drift
cannot masquerade as the matched experiment.

The source-only gate is anchored to the [prior search audit](../mixed-orbit-search-15000-20261012/audit.json)
at cumulative conservative charge 15,000 seconds, SHA-256
`32420a62bb6ccf8f9df7f03dea9040e518606740968110dc8c0ac47cb50cafd5`,
and the [K=600 Rust import receipt](../k600-mixed-orbit-worker-source-20261012/import-receipts.json),
SHA-256 `9d1d8336ba3941ee03b37d67e4da8a447f2abe0c2613090f098d83ed4b936de0`.
The small compressed archive is pinned to SHA-256
`3e85c9e419fcddd3c9faed3c84470bf96fcad9178e48e1d73629369335841ecf`
and BLAKE3
`b4dda45e6c3edfca673db6f3f82c16cce8c860f5e62c092155258df399a0363e`.
The K=2,048 envelope remains pinned to its prior hashes. Build and run
require a clean committed checkout, frozen source files and archives, a
content-hashed authorization record, a locally pinned no-network Linux image,
fresh external output paths and container cleanup receipts.

The authorization schema is `n83.additional-local-budget/v4`, with purpose
`mixed_orbit_k600_two_residual_first_gate`, exact prior-audit and import-receipt
digests, exact build and output paths, an explicit user approval reference,
and at least **2,400 fresh local seconds**. The full conservative charge is
600 build, 360 verification, 1,200 search and 240 cleanup/audit seconds;
each launched stage is charged its declared cap even if it exits early. The
worker wall is 1,150 seconds within the search cap. This ledger anchor is
valid only while no other retained N83 run has changed the 15,000-second
starting total; an intervening run requires a new gate and authorization.
Preparing and testing this source did not launch Docker, SAT, or another
retained N83 search, and did not consume that fresh allowance.

The focused source-only tests exercise authorization rejection before Docker,
exact prior-receipt binding, forged import counts, wrong worker receipts,
stage command and cgroup settings, event counters, point equations, residual
domain membership, source provenance, completed rows, and a two-partial graph
cycle. They use synthetic public points and mocked stage launches. A passing
synthetic control is not a measured K=600 relation result. The next search
would still be a bounded first trial; both bases need natural relation yield,
graph/rank growth, sparse linear algebra, individual-log and matched rho
measurements before total runtime can be ranked.

Validation on the prepared source: the full release library suite passed
4,000 tests with 113 ignored; the N83 study Python suite passed 51 tests and
the boundary autolab suite passed 43. The new five-test supervisor control
is included in the study suite. Python compilation, CLI help and whitespace
checks passed. These checks did not construct or search a retained N83 SAT
instance.
