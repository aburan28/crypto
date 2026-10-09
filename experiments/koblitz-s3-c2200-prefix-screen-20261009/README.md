# Fresh N53 C1200/C2200 prefix call-count screen

The C2200 prefix index recovered and verified every scalar in the 16-run
fresh screen. The eight-point geometric mean C1200/C2200 ratio was **2.377×
for target S3 calls**, with a target-bootstrap 95% interval of
**[1.389, 4.359]**. This is below the frozen **3×** advancement gate, so
this C2200 search policy stops at screening. Its single-repetition local
online ratio was 1.776× with interval [0.957, 3.767]; this is exploratory
and has no A/A noise estimate.

| Public point | C1200 target S3 calls | C2200 prefix target S3 calls | C1200/C2200 calls | Exploratory online ratio |
| --- | ---: | ---: | ---: | ---: |
| T125 | 1042 | 1004 | 1.038 | 0.538 |
| T126 | 2742 | 502 | 5.462 | 3.470 |
| T127 | 444 | 150 | 2.960 | 2.433 |
| T128 | 288 | 128 | 2.250 | 1.303 |
| T129 | 82 | 82 | 1.000 | 0.809 |
| T130 | 1098 | 468 | 2.346 | 1.630 |
| T131 | 94 | 94 | 1.000 | 0.848 |
| T132 | 3824 | 332 | 11.518 | 14.936 |

[`SCREEN_FREEZE.json`](SCREEN_FREEZE.json) fixes eight previously unused
one-target workloads T125–T132 on the same N53 public curve and exact target
point in each pair. It compares C1200 with the C2200 root table and C1200
target-state prefix. Both arms use the same committed source and release
binary; candidate manifests fix the distinct factor bases, ranks, search
orders, and source hashes. Runs alternate AB/BA by point, use 14 preparation
workers and one target worker, and retain all target-dependent online phases
and scalar replay.

This is a **one-repetition screening gate** for deterministic target S3-call
counts. Advancement to a replicated fresh-target wall panel requires all 16
scalars verified and a target-level geometric mean C1200/C2200 call ratio of
at least 3. The threshold, targets, 600-second per-run cap, and 32 GiB
memory envelope were fixed before execution. A wall ratio from this local
screen is exploratory and does not establish the requested 2× online speedup.

Run the frozen screen with:

```sh
python3 experiments/koblitz-s3-c2200-prefix-screen-20261009/run_paired.py \
  --binary target/release/examples/s3_partial_prepared_fullrank_c2200prefix_targetspan \
  --output /private/tmp/s3-c2200-prefix-screen-20261009
```

The runner validates candidate/workload hashes, source and binary hashes,
public points, full-rank preparation, root-table sizes, search policies,
exclusive online phase sum, and independently replayed scalars. Its
append-only ledger retains timeouts, failures, and complete raw reports.
[`local/raw_paired16.jsonl.gz`](local/raw_paired16.jsonl.gz) preserves all
rows; [`local/raw_reports.tar.gz`](local/raw_reports.tar.gz) preserves the
complete reports; [`analyze.py`](analyze.py) verifies every report digest,
pairing, stage status, and scalar result before regenerating
[`analysis.json`](analysis.json). The ledger's decompressed SHA-256 is
`3212fd33714e4205a2ebbbfaa9a55a3b3d1e93d6823b694de2f4c84db195d88c`;
the raw-report archive SHA-256 is
`0a5cac2425f132ff3c308d32fd8c896628e9c2ea25f554f3ba904756c499c469`.
The host record has no exclusive partition receipt and `controlled_speedup`
remains null.
