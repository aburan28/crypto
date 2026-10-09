# Compact retained S3 index states

The full-rank S3 index now retains each root's polynomial `x,y` coordinates and reconstructs normal-basis coordinates during target-independent root insertion and relation preparation. The online target path still uses the same polynomial coordinates, 13-byte dense root table, eight-state inversion batch, and direct four-label scalar recovery. This changes the indexed `State` layout from 72 to 40 bytes.

The [freeze](FREEZE.json) fixes source commit `05123c9ce`, the release binary, N41 C800 and N53 C1200 T025 public-point controls, 14 preparation workers, one target worker, relation seed `20260928`, and a 180-second cap per control. Candidate identities are in [`identities/candidates/`](identities/candidates/). Each control must recover and verify its scalar, match the dense reference's factor base, root table, preparation rank, relation witnesses, S3-call count, target relation, and all other non-timing report fields. Raw failures remain evidence rows.

| Curve | Indexed states | Old state payload | Compact state payload | Exact payload reduction |
| --- | ---: | ---: | ---: | ---: |
| N41 C800 | 26,240,000 | 1,889,280,000 B | 1,049,600,000 B | 839,680,000 B |
| N53 C1200 | 76,320,000 | 5,495,040,000 B | 3,052,800,000 B | 2,442,240,000 B |

These are exact retained-vector payloads. The process peak also includes table slots, temporary root scatter vectors, allocator behavior, code, and other structures. A separately frozen RSS control measures the whole process; this structural reduction alone does not establish a matching peak-RSS reduction or online wall-time gain.

Run each frozen control in a fresh process from the repository root:

```sh
python3 experiments/koblitz-s3-compact-state-20261009/run_control.py \
  --n 41 --binary target/release/examples/s3_partial_prepared_fullrank_compactstate_targetspan \
  --receipt experiments/koblitz-s3-compact-state-20261009/local/n41_receipt.json \
  --raw experiments/koblitz-s3-compact-state-20261009/local/n41_report.json.gz
python3 experiments/koblitz-s3-compact-state-20261009/run_control.py \
  --n 53 --binary target/release/examples/s3_partial_prepared_fullrank_compactstate_targetspan \
  --receipt experiments/koblitz-s3-compact-state-20261009/local/n53_receipt.json \
  --raw experiments/koblitz-s3-compact-state-20261009/local/n53_report.json.gz
```

The T025 controls validate the representation change. Larger-base search is a separate candidate and needs its own frozen one-target inputs, resource cap, full-rank preparation, and independent Sage scalar replay before comparing it with C1200.
