# N41 C1200/C1800 same-source S3 ablation

This experiment measures how much enlarging the S3 factor base and prepared
log system contributes to single-target online performance when the target
kernel is held fixed. **Both arms run the identical source and Mach-O
executable**: the C1800 prefetch/projective implementation, with either
1,200 or 1,800 requested columns. For N41 this source uses the native
rotating-state search schedule in both configurations. It retains the same
batch-eight inversion, root lookup, scalar check, relation seed, 14 setup
workers, and one target worker. The configuration changes the actual base
size, index size, rank, and root-table capacity.

The correctness smoke ran the C1800 executable at C1200 settings on the
previously used T025 control point. It recovered the declared scalar,
reached full rank, and reproduced the earlier C1200 98,400-point factor-base
digest. Its raw report and process streams are in `local/smoke_raw.tar.zst`.
The fresh comparison uses T066–T089, 24 points from the sealed N41 workload
that were unused by either earlier N41 C1200/C1800 panel. `FREEZE.json` binds
their point and workload identities, source and binary hashes, candidate
manifests, alternating pair order, runner and analyzer, 180-second cap per
process, and 32 GiB observed peak-RSS envelope. The primary development
gate is a baseline/C1800 geometric online ratio of at least 2.0 with a
target-bootstrap 95% lower endpoint of at least 2.0. All 48 runs and Sage
scalar replays must verify before a ratio is reported.

The baseline has a new candidate ID because its target-kernel source differs
from the earlier C1200 executable. Its exact manifest is under
`identities/candidates/`; the C1800 candidate retains its prior immutable
manifest. The experiment measures the full one-target online interval through
scalar replay, with reusable setup recorded separately. Local CPU wall-time
ratios remain exploratory under the repository's host-isolation rule.

From the repository root, with the frozen macOS binary present:

```sh
python3 experiments/koblitz-s3-c1800-matchedsource-n41-20261010/run.py \
  --binaries /private/tmp/crypto-s3-release-examples-20261010 \
  --output /private/tmp/s3-n41-c1800-matchedsource24-20261010
python3 experiments/koblitz-s3-c1800-matchedsource-n41-20261010/analyze.py \
  --ledger /private/tmp/s3-n41-c1800-matchedsource24-20261010/raw.jsonl \
  --output experiments/koblitz-s3-c1800-matchedsource-n41-20261010/local/analysis.json
```

`run.py` stops and preserves the completed ledger at the first invalid run.
The analyzer audits each raw report and process stream before computing the
paired target distribution. The checked Sage launcher independently replays
each recovered scalar on the declared N41 curve.

## Result

All 48 fresh-process runs reached their declared pre-target rank, passed the
native group check, and recovered the frozen scalar. Checked Sage replayed
all 48 scalars against the 24 public points. The source and executable hashes
match exactly between arms in every row.

| N41 measure | C1200, same executable | C1800, same executable |
| --- | ---: | ---: |
| Actual usable base points | 98,400 | 147,600 |
| Prepared rank / verified preparation relations | 1,200 / 2,534 | 1,800 / 3,080 |
| Root-table slots | 134,217,728 | 268,435,456 |
| Median target S3 calls | 32 | 32 |
| Median target query, ms | 0.012541 | 0.014334 |
| Median charged target PDP, ms | 0.027542 | 0.023542 |
| Median scalar recovery check, ms | 0.003646 | 0.003771 |
| Median online wall time, ms | 0.063792 | 0.064188 |
| Median reusable setup wall time, ms | 6,257.987 | 13,860.214 |
| Maximum observed peak RSS, GiB | 8.29 | 17.43 |

The paired C1200/C1800 geometric online ratio is **0.997188**, with a
20,000-resample target-bootstrap 95% interval of **[0.886713, 1.119985]**.
The geometric target S3-call ratio is 1.218287. The frozen 2× development
and lower-endpoint gates did not pass. On these 24 N41 points, enlarging the
S3 base and index roughly doubled reusable setup and peak memory while the
same-kernel online times stayed close. This isolates the base/index effect
from the target-kernel improvement present in the earlier cross-source
candidate comparison; quantifying the kernel contribution on identical
points is the next ablation.

The five exclusive online phases, raw failures, and resource counters are in
`local/analysis.json` and the raw archive. Independent phase medians above
need not sum to the median online total; every individual report has a
checked phase sum. `local/raw_matchedsource24.tar.zst` contains all 48
reports, stdout, stderr, and the ledger; its SHA-256 is
`6824185505ace2c5c59807bdc0f1a59ce834e58c53bd5f55931e885f3ddd9dca`.
`local/sage_runtime_info.json` and `local/sage_replay.json` bind the checked
Sage execution to the same ledger. CPU timing remains exploratory without a
host-exclusive partition receipt.
