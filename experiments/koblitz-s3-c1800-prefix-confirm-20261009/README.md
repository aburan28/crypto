# Independent N53 C1800 prefix confirmation

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
