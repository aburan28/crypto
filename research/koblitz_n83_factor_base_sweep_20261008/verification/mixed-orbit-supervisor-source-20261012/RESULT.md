# Mixed-orbit supervisor source gate

The supervisor source prepares a first retained-N83 two-residual partial
producer gate using the already stored K=1,182 and K=2,048 prefix-related
bases. A separate run-specific authorization for at least 2,400 additional
local seconds is required before its `build` or `run` mode can launch Docker.
The fully charged prior allowance remains 11,700 seconds. This source gate
did not launch either mode, import the N83 objects, or search for relations.

The supervisor pins the existing Linux arm64 image, both object SHA-256
digests, committed source, matched public corpus, one trial and its hard
resource limits. It retains stage commands, exit/timeout/OOM/cleanup states,
raw logs, frozen input hashes and conservative charges. It validates the
worker's import and search receipts, checks per-trial counters against the
event ledger, and independently replays any partial and completed row using
binary-curve group arithmetic. Censored solver, graph and resource outcomes
remain separate from producer failures. The frozen [first-run protocol](PROTOCOL.md)
sets the exact gate and its accounting; it is not a runtime ranking.

The source tests exercise authorization rejection before Docker, command
construction and stage handoff under mocks, malformed-receipt rejection, and
positive/negative independent group replay. The study Python suite passed
43/43 tests. The first four-thread full release-library run passed 3,999
tests and failed one existing WDSat mock-child timing check: its nominal
2,000 ms deadline expired under suite contention. That exact test passed
1/1 in isolation. The full release-library suite then passed serially with
4,000 tests, zero failures and 113 ignored. The [machine receipt](validation.json)
binds final source hashes and all four raw logs, including the failed run.
No retained-N83 relation yield, independent rank or complete cold runtime is
measured here.
