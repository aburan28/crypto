# Fresh N53 C1800 striped-search screen

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
online phase sum, and independently replayed scalar for every run.
