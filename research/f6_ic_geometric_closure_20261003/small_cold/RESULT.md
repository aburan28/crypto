# Full cold n9 IC control: F6-IC recovers one public logarithm

The [precommitted protocol](../SMALL_COLD_PROTOCOL.md) fixed one previously
unseen public point `[305,466]` on registered curve
`icv1-f2m9-tm5-4a3ea183`, subgroup order 37, from fixture seed
`20261004041`. No scalar was constructed or supplied to a solver. The
workload ID is `f6d79f6dd9f2`; the exact field/curve, 21 geometric base
points, 14 distinct nonidentity subgroup-usable points, two folded
columns, three candidate IDs and input hashes are in the [freeze
receipt](freeze.txt). The three candidates differ in their PDP solver;
each uses native ordinary relation collection and dense final relation LA.

| IC arm | Inside-worker cold ms | One-target online ms | Cold F4 / arm | Online F4 / arm | Verified log |
| --- | ---: | ---: | ---: | ---: | ---: |
| Inherited F4 | 4.707 | 0.372 | 1.000× | 1.000× | 4 |
| F6-IC | **3.150** | **0.272** | **1.495×** | **1.366×** | 4 |
| Matrix F5 | 61.937 | 3.118 | 0.076× | 0.119× | 4 |

All three processes exited zero. Each built its base from scratch, issued
eight ordinary queries, produced eight distinct accepted and independently
checked relations, obtained and verified logs for both folded columns in
one dense matrix solve, then found the target relation on its first
target-dependent attempt. The source's group replay recovered `4` in all
three. A separate [native replay certificate](replay-certificate.json)
recomputed `[4]G=Q`, all 24 ordinary relation identities and six column
log identities and checked every cold and online phase sum. Its SHA-256 is
`307f02ce84e5cb463257c5bfc9a572aed8508a90f8c9b58e1c3d550589bb0d90`.
The [three keyed rows](measurements.jsonl) and [raw runs](runs/) retain the
relation matrices, query outcomes, phase ledgers, stderr and process status.

The cold interval starts inside the worker after process launch and JSON
parsing, before curve/base construction, and ends after scalar replay. Its
exclusive phases sum exactly: F6 spent 2.490 ms (79.1%) on ordinary PDP,
0.262 ms (8.3%) on target PDP, and 0.398 ms on the remaining stages,
including base construction, matrix build and final LA. Its online interval
begins after reusable preparation and ends after target replay; 96.2% of
that interval was target PDP. Thus the F6/F4 difference persists when the
ordinary relation and matrix work is charged, on this one tiny instance.

This is a **correctness and cost-attribution control**, with one target and
one process per arm on an unisolated Mac. The wall ratios are exploratory;
they give no target-to-target interval. The group has only 37 elements and
the folded relation matrix only two columns, so this result cannot be
extrapolated to larger curves. The native `ecbench` IC adapter does not yet
expose this F6 oracle. No same-target strong-rho session or independent
cross-host replay was made, leaving the primary IC/rho ratio, operation-
normalized `S`, and any attack improvement **unknown**. The current n53
`PDP4root` path does not use F6-IC.
