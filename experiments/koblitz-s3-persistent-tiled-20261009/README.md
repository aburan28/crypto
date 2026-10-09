# Persistent S3 target tile workers

One 14-worker pool preserves the fixed serial-order tile search and its first
verified target relation across all 48 N53 paired executions. Checked Sage
independently recovered all 48 scalars and verified 462 distinct relation
witnesses. The persistent implementation repeated exactly the same target
S3-call count and relation indices as the per-tile worker implementation on
each of eight public points. Its exploratory online wall time was slower:
the paired per-tile/persistent geometric ratio was 0.490×. The candidate
does not advance to the T001–T024 panel.

## Frozen pipeline

The curve `EC1N53Ckb1h738a572b91f2` uses polynomial-basis `GF(2^53)`
modulo `x^53+x^6+x^2+x+1`, the equation `y²+xy=x³+1`, a prime subgroup
of order 21044858204113, 25,864 usable factor-base points, 244 folded
columns, and 243 prepared rows. The
[parent freeze](../koblitz-s3-crosscurve-20261009/FREEZE.json) fixes the
public targets, input law, and workload IDs. The [experiment freeze](FREEZE.json)
fixes T025–T032, 14 query workers, 2,048-state tiles, two repetitions per
arm and target, source and binary hashes, and candidate IDs. The persistent
[candidate manifest](identities/candidates/IC1N53Ckb1fb25864PDP4rootRCsampleLAgaussTDdirectISO0hcb32e91f4c92.json)
retains the full algorithm identity.

| Arm | Candidate ID | Target query policy |
| --- | --- | --- |
| Immediate parallel | `IC1N53Ckb1fb25864PDP4rootRCsampleLAgaussTDdirectISO0he884cf865574` | First worker to verify a relation |
| Per-tile workers | `IC1N53Ckb1fb25864PDP4rootRCsampleLAgaussTDdirectISO0h031b01b9610a` | Spawn workers for every serial-order tile |
| Persistent workers | `IC1N53Ckb1fb25864PDP4rootRCsampleLAgaussTDdirectISO0hcb32e91f4c92` | Spawn once, synchronize every serial-order tile |

The persistent pool divides each consecutive 2,048-state tile among the 14
workers. Each worker reports its first verified relation; the coordinator
selects the smallest serial offset and proceeds to the next tile only when
none was found. All worker startup, barriers, and draining occur inside the
single target's online interval. The clock begins after reusable base, index,
and relation preparation; it includes target validation, query, S3 point
decomposition, relation check, target descent, and scalar replay. Target
fixture construction, process launch, and input loading are outside it.

## Paired result

The [canonical ledger](local/canonical_paired48.jsonl.gz) contains one row
per candidate, workload, and run ID, including phase costs, failures, and
memory. The [analysis](local/analysis.json) checks the exact identifiers,
target hashes, binaries, phase sums, relations, and root calls. Every native
run completed and verified the supplied public point. The
[Sage replay](local/sage_replay.json) independently checks group equations,
prepared relation rank, and scalar recovery using the accepted
[runtime receipt](local/sage_runtime_info.json).

| Paired measure on T025–T032 | Result |
| --- | ---: |
| Verified one-target runs | 48/48 |
| Per-tile/persistent online geometric ratio | 0.490× |
| Immediate/persistent online geometric ratio | 0.213× |
| Targets where per-tile is faster | 3/8 |
| Targets where immediate is faster | 8/8 |
| Median online time, immediate / per-tile / persistent | 5.958 / 13.859 / 27.921 ms |
| Median target S3 calls, immediate / per-tile / persistent | 14,365 / 28,511 / 28,511 |
| A/A p95 absolute log wall ratio, immediate / per-tile / persistent | 1.802 / 1.888 / 2.233 |
| Per-tile/persistent target relation and S3-call equality | 8/8 points, 32/32 runs |

The last equality compares the two tiled arms' 16 runs each. Their root-call
counts are also identical across repetitions for every point. The immediate
first-winner policy searches a different, schedule-dependent amount of work.
Wall times on this shared host are exploratory because the A/A noise gate
fails and no host-wide CPU isolation receipt is available. The controlled
speedup field remains unknown. This trial isolates the persistent pool's
barrier and scheduling cost as a poor trade for eliminating per-tile thread
creation on these inputs.

The persistent source SHA-256 is
`94603624c42c1d8642050c0c5ba185cccd3ff35cd330fb29e0a08436e82a4bb3`,
its local release binary SHA-256 is
`3743d455c9307e880411de5a93de584eb15a9913737681d7d269be618044a766`,
the plan SHA-256 is
`e2c60f0275ee61ae474e6b73c51e7ea76f29e5366ef04ea37ff8ebdadcc550c2`,
and the uncompressed canonical ledger SHA-256 is
`a8e6b9549f1d4cebf556c84a8548809d8be0567dd51ee91f486c53df789fccab`.
