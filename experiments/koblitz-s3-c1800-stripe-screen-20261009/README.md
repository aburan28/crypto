# Fresh N53 C1800 striped-search screen

The C1800 striped root search recovered and verified every scalar in the
16-run fresh screen. Across eight paired public points, the C1200/striped
geometric mean was **3.505× fewer target S3 calls** for the striped arm,
with a target-bootstrap 95% interval of **[1.821, 6.916]**. This passes the
frozen **3×** point-estimate advancement gate. The single-repetition local
online ratio was **2.809×** with interval **[1.621, 5.042]**. It is an
exploratory timing diagnostic: there is no A/A noise estimate or exclusive
host partition receipt for this screen. A new replicated target panel is
the next gate for the requested 2× online wall result.

| Public point | C1200 target S3 calls | C1800 striped target S3 calls | C1200/striped calls | Exploratory online ratio |
| --- | ---: | ---: | ---: | ---: |
| T133 | 4206 | 1134 | 3.709 | 2.922 |
| T134 | 2532 | 766 | 3.305 | 10.055 |
| T135 | 396 | 126 | 3.143 | 1.650 |
| T136 | 3738 | 172 | 21.733 | 7.945 |
| T137 | 1684 | 188 | 8.957 | 4.465 |
| T138 | 2344 | 772 | 3.036 | 2.405 |
| T139 | 12 | 12 | 1.000 | 0.966 |
| T140 | 190 | 190 | 1.000 | 0.971 |

[`SCREEN_FREEZE.json`](SCREEN_FREEZE.json) fixes eight previously unused
one-target public points T133–T140. Each point is paired between the native
C1200 state order and the C1800 root table with a 256-state exact lead and
eight interleaved C1200 state ranges. Both arms use one committed release
binary and source snapshot, 14 preparation workers, one target worker, and
the exact factor-base and search-policy manifests frozen before timing.

This one-repetition screen advances only if all 16 scalar recoveries verify
and the target-level geometric mean C1200/striped target S3-call ratio is
at least 3. AB/BA order alternates by target. A 180-second per-run timeout,
32 GiB memory envelope, and all failure rows are retained. Local online
wall-time ratios are exploratory; the separate replicated panel and host
isolation gate are needed for the requested 2× wall claim.

Run the frozen screen with:

```sh
python3 experiments/koblitz-s3-c1800-stripe-screen-20261009/run_paired.py \
  --binary target/release/examples/s3_partial_prepared_fullrank_c1800stripe_targetspan \
  --output /private/tmp/s3-c1800-stripe-screen-20261009
```

The runner verifies candidate and workload hashes, source/binary identities,
public points, full-rank factor-base logs, table size, search schedule,
online phase sum, and independently replayed scalar for every run. Its
append-only ledger retains complete reports and failures.
[`local/raw_paired16.jsonl.gz`](local/raw_paired16.jsonl.gz) preserves all
16 rows; [`local/raw_reports.tar.gz`](local/raw_reports.tar.gz) preserves the
full reports. [`analyze.py`](analyze.py) rechecks report hashes, matching
points, frozen identities, stages, and scalar verification to regenerate
[`analysis.json`](analysis.json). The ledger's decompressed SHA-256 is
`4554bf730d67d4bb5bfe6da3fa738a403444ca6246f6937a4f2439d1fcfbe1b5`;
the report archive SHA-256 is
`e8f7648e7f963ea8b00006523cab6172608e5ca3cc15f625f04543e4eb69afc4`.
The host record has no exclusive partition receipt and `controlled_speedup`
remains null.
