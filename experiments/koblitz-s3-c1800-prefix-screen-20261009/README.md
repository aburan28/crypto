# N53 C1800 prefix search: paired fresh-target screen

The C1200-prefix traversal over the full C1800 root index passed its frozen
advancement gate. On eight previously unused N53 public targets, the
geometric mean C1200/C1800 ratio was **2.954 for target S3 calls**. All 32
paired runs recovered and verified the expected scalar. The screen gate in
[`SCREEN_FREEZE.json`](../koblitz-s3-c1800-prefix-20261009/SCREEN_FREEZE.json)
was a call-count ratio of at least 2, fixed before these targets were run.

| Public point | C1200 target S3 calls | C1800 prefix target S3 calls | C1200/C1800 calls | Exploratory online ratio |
| --- | ---: | ---: | ---: | ---: |
| T109 | 1412 | 404 | 3.495 | 4.679 |
| T110 | 1822 | 964 | 1.890 | 1.655 |
| T111 | 924 | 924 | 1.000 | 0.816 |
| T112 | 1668 | 1182 | 1.411 | 1.282 |
| T113 | 2724 | 196 | 13.898 | 6.272 |
| T114 | 1598 | 280 | 5.707 | 4.195 |
| T115 | 3060 | 390 | 7.846 | 6.238 |
| T116 | 694 | 694 | 1.000 | 0.962 |

Each table row pairs the same public point, relation seed, resource envelope,
and one-target online interval. It combines two repetitions per arm by their
geometric mean for online time; deterministic S3-call counts agreed between
repetitions. The online interval begins after reusable factor-base, S3 index,
and full-rank log preparation, and includes all target query, decomposition,
relation check, descent, and scalar replay phases. Preparation costs and all
phase records remain in the raw reports. C1200 has 127,200 actual factor-base
points and rank 1200; C1800 has 190,800 and rank 1800. Its T025 correctness
control independently recovered the same scalar under both candidates, and
the new C1200 implementation exactly matched the prior C1200 report outside
timing and the new schedule labels.

Target-level bootstrap resampling (20,000 resamples, seed 20261009) gave a
95% interval of `[1.584, 5.628]` for the call ratio. The exploratory online
ratio was 2.446 with interval `[1.402, 4.197]`. Repetition A/A 95th
percentiles of absolute log online-time ratios were 1.882 for C1200 and 0.613
for C1800. The macOS host has no exclusive CPU partition receipt. Therefore
the online wall ratio is a local diagnostic and `controlled_speedup` is null.
The next gate is an independently frozen T117–T124 confirmation on the same
candidate and then a host-isolated paired online wall panel.

[`local/raw_paired32.jsonl.gz`](local/raw_paired32.jsonl.gz) contains every
run, failure status, phase, witness, verification flag, and source/binary
receipt. [`local/raw_reports.tar.gz`](local/raw_reports.tar.gz) retains all 32
complete reports; [`local/host.json`](local/host.json) records the host.
[`analyze.py`](analyze.py) verifies all report SHA-256 digests and candidate
identities and regenerates [`analysis.json`](analysis.json). The ledger's
decompressed SHA-256 is
`e9e9b18c2760af925ba8af3e91bb3e897af53b9bce36e69e71c132f7f2c1983b`;
the report archive SHA-256 is
`4472ee29245cc28ae1330cdc0ca8a22c5ceb5602f6e00bd1a2554105c813ea52`.
[`control_receipt.json`](control_receipt.json) records the T025 checks and
separate compressed control reports.
