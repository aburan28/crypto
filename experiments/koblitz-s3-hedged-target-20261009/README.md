# Delayed S3 target-worker launch on N53

The immediate 14-worker target search remained the fastest tested policy on
the separately frozen N53 T025–T032 public points. Starting one search at
once and delaying the other 13 by 0, 1,000 or 3,000 microseconds produced
exploratory baseline/hedged online ratios of **0.920×, 0.597× and 0.443×**.
All 64 one-target executions recovered their scalars, and the checked Sage
build independently recovered all 64 while checking 560 distinct relation
witnesses. The delayed policy was not promoted to the T001–T024 evaluation
set.

## Exact pipeline and frozen inputs

The curve is `EC1N53Ckb1h738a572b91f2`, with polynomial-basis field
`GF(2^53)` modulo `x^53+x^6+x^2+x+1`, equation `y²+xy=x³+1`, prime
subgroup order 21044858204113, 25,864 actual subgroup-usable factor-base
points, 244 folded relation columns and 243 prepared rows. The input law,
generator, public points and workload IDs are fixed in the
[parent plan](../koblitz-s3-crosscurve-20261009/FREEZE.json). The
[tuning freeze](FREEZE.json) selects T025–T032 and fixes the implementation
and three delay values. The three immutable candidate manifests are in
`identities/candidates/`.

| Arm | Candidate ID | Delayed peers |
| --- | --- | ---: |
| Immediate parallel reference | `IC1N53Ckb1fb25864PDP4rootRCsampleLAgaussTDdirectISO0he884cf865574` | 0 |
| Hedged 0 µs | `IC1N53Ckb1fb25864PDP4rootRCsampleLAgaussTDdirectISO0h0e0f31adb8e6` | 0 |
| Hedged 1,000 µs | `IC1N53Ckb1fb25864PDP4rootRCsampleLAgaussTDdirectISO0hde426773d72e` | 1,000 µs |
| Hedged 3,000 µs | `IC1N53Ckb1fb25864PDP4rootRCsampleLAgaussTDdirectISO0h349cc246bea6` | 3,000 µs |

The hedged implementation lets worker 0 search immediately. A helper waits
for the declared delay, launches 13 other rotating index starts if needed,
and is woken early when worker 0 finishes. A verified relation atomically
cancels remaining searches. Every worker drains before the online clock
stops. All target-dependent attempts, startup, cancellation and scalar
replay are charged to the one public point. Base and index preparation,
input loading and fixture construction are outside this interval.

## Paired exploratory sweep

Each of the eight points was solved twice by each arm, with order shuffled
using its frozen index. The [canonical run ledger](local/canonical_probe64.jsonl.gz)
retains all 64 reports under `(candidate_id, workload_id, run_id)`, including
failures, the stage clock, operation counts and exact point identity. Every
run completed and passed the phase-sum and native scalar checks. The
[analysis](local/canonical_analysis.json) reports target-level ratios and
A/A variation; the [Sage replay](local/sage_replay.json) verifies target
relations, prepared relation equations, matrix rank and each recovered
scalar. Its [runtime receipt](local/sage_runtime_info.json) records the
checked installed Sage build.

| Candidate delay | Exploratory online ratio, immediate / hedged | Points above 1× | S3 root-call ratio, immediate / hedged | Hedged A/A p95 absolute log ratio |
| --- | ---: | ---: | ---: | ---: |
| 0 µs | 0.920× | 4/8 | 0.945× | 1.195 |
| 1,000 µs | 0.597× | 2/8 | 0.737× | 1.604 |
| 3,000 µs | 0.443× | 1/8 | 0.909× | 1.334 |

The immediate reference A/A p95 absolute log ratio was 5.452, including
one 581.5 ms repeat where the same point took 2.49 ms on its other repeat.
The host has no isolation receipt. These timing ratios are exploratory;
the root-call results and the lack of a consistent point-level gain give
no basis to advance a delayed candidate. Concurrent first-winner search
also has schedule-dependent root-call counts, which a future fidelity
analysis must separate from host contention.

The hedged source SHA-256 is
`0d288d5f49f03d4c80eabe330c74ab959d6910a69ac585f4f409ac5a426d0e06`;
the local binary SHA-256 is
`71df02bd10284b67068c45cae3e2110e720288b517505e2f662f9b7066adc794`.
The tuning plan SHA-256 is
`44ca130e4a5748921809b2aa7ce7ef8ca0768b8d69a00b8849e582473127019d`,
and the uncompressed canonical ledger SHA-256 is
`ff0ff697de3698e289079f8a034003c975e0d02ed73d0e5f59bddc6a08a78e2b`.

The next step is a Modal fidelity measurement of the immediate parallel
candidate against group-lift, using the already prepared
[panel](../koblitz-s3-parallel-target-20261009/modal_parallel_panel.py) and
[analysis](../koblitz-s3-parallel-target-20261009/analyze_modal_parallel.py).
