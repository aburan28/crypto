# Frozen first cross-stream N83 gate

This directory defines a single retained K600 public-X hash base versus a
K600 public-X sequential pool. The base is the already audited object
`b4dda45e6c3edfca673db6f3f82c16cce8c860f5e62c092155258df399a0363e`;
the pool is
`e5600aa76a88246a4c3f73fce83ba73479e79c34e88d6205fffdb91ce09e0a99`.
Both are 600-column, 99,600-point signed-Frobenius pilot archives with seed
2026100801. The pool HEAD was observed at 3,020,387 compressed bytes with
ETag `9e41c76de859046eb295aec39e635c3c`; these metadata alone do not prove
archive contents. The importer verifies the manifest BLAKE3 and complete point
and group structure. The supervisor independently replays the archive labels,
exact signed-orbit set difference, bounded result identity, and any candidate
group equation.

The frozen sequence is an offline four-CPU/eight-GiB Linux arm64 build, then
one-CPU/four-GiB/no-swap import and one-trial natural-target partial search.
The conservative stage charges are 600 seconds build, 360 import, 1,200 search,
and 240 cleanup, totaling 2,400 additional local seconds. Earlier K600 work
is charged 17,400 seconds; a complete use of this gate would bring that
ledger to 19,800 seconds. The one-trial search uses three base points and two
residual points, at most 5,000 conflicts, 32 models per identity pattern, eight
patterns, and an 1,100-second internal wall watchdog. Only exact overlap after
the import determines the actual residual-domain size. The capacity screen
for all 648 primary pilot pairs is conditional source arithmetic, not a
runtime or yield measurement.

`supervisor.py` requires a clean commit, the pinned prior audit, the two local
archive objects, a pinned offline container image and registry, and a fresh
pair-specific authorization record before it creates an output directory. It
retains stdout/stderr, Docker state, source/input/binary digests, a
conservative budget receipt, and replay status outside the checkout. Timeout,
resource cap, import failure, invalid model, and unfinished search are distinct
statuses. This gate does not collect enough relations or linear algebra to
rank complete index-calculus runtime.

The source-only controls are in `test_supervisor.py`. The retained run used
source commit `0764901992591c39566d720ae812c54cdc09ab75` and kept its
raw output at
`/Volumes/SSD990/crypto/worktrees/.n83-cross-stream-first-gate-20261013/run-0764901992`.
The independent `audit.py` replay accepted the complete raw receipt set; its
hash inventory and verdict are frozen in `audit.json`. Both objects remain in
the specified S3 factor-base prefix under `a0/objects/<BLAKE3>.jsonl.gz`.

| Measured field | First cross-stream gate |
| --- | ---: |
| Shared signed orbits | 0 |
| Base / residual signed orbits | 600 / 600 |
| Base / residual points | 99,600 / 99,600 |
| Source capacity variables / domain clauses / XOR rows | 97,080 / 269,830 / 4,980 |
| Observed build / verify / search stage wall, seconds | 194.192 / 3.717 / 275.429 |
| Observed supervisor wall, seconds | 486.233 |
| Search process wall / import, seconds | 274.522 / 3.281 |
| Search memory cgroup peak, bytes | 369,958,912 |
| Search conflicts / models / solver calls | 5,000 / 0 / 1 |
| Solver status | `UNKNOWN_search_cap` |
| Conservative additional / cumulative charge, seconds | 2,400 / 19,800 |

The first of eight possible identity patterns reached its 5,000-conflict cap.
This is a censored search observation: it supplies neither a zero-yield
estimate nor a refutation of the relation equation. The run produced no graph
rows, rank, or complete index-calculus runtime, so the boundary ledger and
selected-best field remain unchanged. No additional retained search is
authorized by this grant.
