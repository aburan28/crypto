# Independent N53 C1800 prefix confirmation

The unchanged C1800-prefix candidate passed the independent, predeclared
target S3-call gate on T117–T124: **2.029×** C1200/C1800 geometric mean,
with a target-bootstrap 95% interval of **[1.279, 3.357]**. All 32 runs
recovered and verified their scalars. This confirms an algorithmic call
reduction; the exploratory online wall ratio on this host was **1.295×**
with interval **[0.904, 1.827]**. A larger call reduction or a qualifying
isolated host panel is needed to establish the requested 2× online wall
threshold.

| Public point | C1200 target S3 calls | C1800 prefix target S3 calls | C1200/C1800 calls | Exploratory online ratio |
| --- | ---: | ---: | ---: | ---: |
| T117 | 550 | 550 | 1.000 | 1.396 |
| T118 | 260 | 260 | 1.000 | 1.090 |
| T119 | 1414 | 672 | 2.104 | 1.254 |
| T120 | 4320 | 906 | 4.768 | 2.081 |
| T121 | 762 | 472 | 1.614 | 0.512 |
| T122 | 462 | 190 | 2.432 | 1.927 |
| T123 | 1196 | 164 | 7.293 | 2.687 |
| T124 | 462 | 462 | 1.000 | 0.752 |

[`CONFIRM_FREEZE.json`](CONFIRM_FREEZE.json) fixes a separate eight-point
panel on the previously unused T117–T124 public targets. It retains the
unchanged C1200 and C1800-prefix candidate manifests, source and release
binary hashes, one-target online interval, 14 preparation workers, one
target worker, two repetitions per arm, ABBA/BAAB order, 180-second run cap,
and 32 GiB memory envelope from the earlier T109–T116 screen. The run starts
only after this protocol is committed.

The independent confirmation gate requires all 32 scalars verified, a
target-level geometric mean C1200/C1800 target S3-call ratio of at least 2,
and the lower endpoint of a 20,000-resample, target-cluster bootstrap 95%
interval strictly above 1. The target labels and these thresholds were fixed
before any confirmation timing. Wall-time ratios remain exploratory on this
local host until an exclusive partition receipt and noise gates are available.

Run the committed protocol with:

```sh
python3 experiments/koblitz-s3-c1800-prefix-confirm-20261009/run_confirm.py \
  --binary target/release/examples/s3_partial_prepared_fullrank_c1200prefix_targetspan \
  --output /private/tmp/s3-c1800-prefix-confirm-20261009
```

The runner checks candidate and workload hashes, the exact target point,
full-rank preparation, factor-base and root-table identities, stage timing
sum, scalar replay, and search schedule. It appends all successful and
failed rows and preserves each complete raw solver report.

[`local/raw_paired32.jsonl.gz`](local/raw_paired32.jsonl.gz) retains all 32
complete run rows; [`local/raw_reports.tar.gz`](local/raw_reports.tar.gz)
retains every raw report. [`analyze.py`](analyze.py) checks their digests,
candidate fields, point pairing, phase status, and native verification before
regenerating [`analysis.json`](analysis.json). The ledger's decompressed
SHA-256 is
`3e9b7dde95055ce45d3da90b6011e614e70c3c074464ce93145579d82c5af519`;
the report archive SHA-256 is
`9addb9718b0d9cafc456c417245a3b2152e785b01d2c9f7e4c687bb44216d868`.
The host record has no exclusive partition receipt. Its A/A 95th percentiles
of absolute log online-time ratios were 0.541 for C1200 and 1.918 for
C1800, and `controlled_speedup` remains null.
