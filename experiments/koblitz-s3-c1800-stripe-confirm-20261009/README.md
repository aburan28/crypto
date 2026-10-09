# Independent N53 C1800 striped confirmation

[`CONFIRM_FREEZE.json`](CONFIRM_FREEZE.json) fixes eight previously unused
public targets T141–T148, each paired on the same N53 curve and point between
the native C1200 state order and the C1800 root table with a 256-state lead
and eight striped C1200 state ranges. The candidate manifests, source and
release binary hashes are unchanged from the successful T133–T140 screen.
This protocol was committed before any confirmation timing.

Each target uses ABBA or BAAB order, two repetitions per arm, 14 preparation
workers, one target worker, a 180-second per-run cap, and a 32 GiB memory
envelope. The independent call-count confirmation gate requires all 32
scalar recoveries to verify, a target-level geometric mean C1200/striped
S3-call ratio of at least 2, and the lower endpoint of a 20,000-resample
target-bootstrap 95% interval strictly above 1. A separate exploratory
online gate requires an online ratio of at least 2 and a bootstrap lower
endpoint above 1. The local wall result remains exploratory without a
host-wide isolation receipt and passing noise gates. Preserve failures and
all target-dependent phase timings in the append-only ledger.

Run the committed protocol with:

```sh
python3 experiments/koblitz-s3-c1800-stripe-confirm-20261009/run_confirm.py \
  --binary target/release/examples/s3_partial_prepared_fullrank_c1800stripe_targetspan \
  --output /private/tmp/s3-c1800-stripe-confirm-20261009
```

The runner checks candidate and workload identities, the exact point,
full-rank preparation, factor-base and root-table identities, search
schedule, exclusive online phase sum, and independently replayed scalar.
[`preserve_confirm.py`](preserve_confirm.py) archives complete ledger rows and
raw reports. [`analyze.py`](analyze.py) checks their hashes, pairing, stage
status, and verification before regenerating the paired statistics.
