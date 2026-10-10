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
