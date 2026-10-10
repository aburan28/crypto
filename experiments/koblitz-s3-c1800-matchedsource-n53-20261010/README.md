# N53 C1200/C1800 same-source S3 screen

This screen asks whether the larger N53 S3 factor base and index improve
one-target online recovery with the target implementation held fixed. Both
arms use the **same prefetch/projective source and Mach-O executable**.
C1200 uses 127,200 subgroup-usable factor-base points and 1,200 relation
columns; C1800 uses 190,800 points and 1,800 columns. The C1800 path first
searches its embedded C1200 state prefix, then the complete C1800 index.
The same relation seed, 14 setup workers, one target worker, batch-eight
inversion policy, and scalar-check implementation apply to both arms.

The T205 correctness smoke ran that executable at C1200 settings. It
reached full rank, reproduced the earlier C1200 factor-base digest, and
recovered the public-point scalar. Its report and process streams are in
`local/smoke_raw.tar.zst`. The eight public points T205–T212 were already
used in the earlier cross-source screen; this paired panel is a new
source-matched ablation on those frozen inputs. `FREEZE.json` binds their
point and workload identities, exact candidate manifests, source and
binary SHA-256, ABBA/BAAB order with two fresh-process repetitions per arm
and point, runner and analyzer, 180-second process cap, and 32 GiB observed
RSS envelope.

The development gate is a geometric C1200/C1800 online ratio of at least
2.0 and a target-plus-repeat bootstrap 95% lower endpoint of at least 2.0.
The analyzer requires all 32 native solves to verify. Checked Sage then
independently multiplies every recovered scalar by the exact N53 generator.
CPU timing remains exploratory without a host-exclusive partition receipt;
the separately frozen 34-point N53 confirmation remains the fresh-target
follow-up.

From the repository root with the exact frozen macOS executable present:

```sh
python3 experiments/koblitz-s3-c1800-matchedsource-n53-20261010/run.py \
  --binaries /private/tmp/crypto-s3-release-examples-20261010 \
  --output /private/tmp/s3-n53-c1800-matchedsource8-20261010
python3 experiments/koblitz-s3-c1800-matchedsource-n53-20261010/analyze.py \
  --ledger /private/tmp/s3-n53-c1800-matchedsource8-20261010/raw.jsonl \
  --output experiments/koblitz-s3-c1800-matchedsource-n53-20261010/local/analysis.json
```

## Result

All 32 fresh-process runs reached complete pre-target rank, passed native
group verification, and recovered the frozen scalar. Checked Sage
independently replayed all 32 scalars against the eight public points. The
same executable and target kernel were used in both arms.

| N53 measure | C1200, same executable | C1800, same executable |
| --- | ---: | ---: |
| Actual usable base points | 127,200 | 190,800 |
| Prepared rank / verified preparation relations | 1,200 / 2,338 | 1,800 / 4,970 |
| Root-table slots | 268,435,456 | 536,870,912 |
| Median target S3 calls | 584 | 192 |
| Median charged target PDP, ms | 0.1407705 | 0.0547085 |
| Median online wall time, ms | 0.174417 | 0.098250 |
| Median reusable setup wall time, ms | 7,071.720 | 18,072.961 |
| Maximum observed peak RSS, GiB | 11.19 | 23.44 |

The paired C1200/C1800 geometric online ratio is **2.056889**. A 20,000
resample bootstrap over target identities gives **[1.246035, 3.738364]**;
resampling one of the two timings per arm within each selected target gives
the repeat-aware 95% interval **[1.219081, 3.757180]**. The frozen 2×
point-estimate gate passed and the lower-endpoint gate did not. The
geometric target S3-call ratio is **2.689628**. The larger N53 index
reduces decomposition work on some targets substantially, while other
targets see similar online time; its setup and memory costs are higher.

Independent phase medians need not add to the median total; each report has
a checked five-phase sum. `local/raw_matchedsource8.tar.zst` contains all
32 reports, stdout, stderr, and the ledger. Its SHA-256 is
`0648d18a1976ced539e2ab0b2a2ce5143e70c7307b608ce13ab8e86bfbd8cef3`.
`local/analysis.json`, `local/sage_runtime_info.json`, and
`local/sage_replay.json` bind the result and checked Sage replay to the
frozen inputs. The 34 fresh N53 targets are the next distribution check;
the bounded Modal pilot remains prepared pending the reviewed transfer.
