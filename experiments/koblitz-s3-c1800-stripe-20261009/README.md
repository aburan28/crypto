# C1800 prefix-256 striped target search

The candidate retains the first 256 point-seeded C1200 index states in their
original order. It then visits the remaining embedded C1200 states by
round-robin across eight distant contiguous ranges and finally falls back
to a complete C1800 traversal. It queries the same C1800 root table and
full-rank factor-base logs as the earlier prefix candidate. A permutation
test verifies that no embedded or fallback state is omitted or duplicated;
the first 256 states have exact old-order equality.

[`DIAGNOSTIC_FREEZE.json`](DIAGNOSTIC_FREEZE.json) fixes the previously used
T110 point as a correctness and search-order control, with one run per arm.
It cannot estimate fresh-target yield or a wall-time speedup. The new
candidate manifests in [`identities/candidates`](identities/candidates)
fix the source, root table, state order, and exact factor bases. If both
scalars verify, a separate fresh T133–T140 screen will test the search policy
against the original C1200 order on one public point per paired run.

Run the committed diagnostic with:

```sh
python3 experiments/koblitz-s3-c1800-stripe-20261009/run_diagnostic.py \
  --binary target/release/examples/s3_partial_prepared_fullrank_c1800stripe_targetspan \
  --output /private/tmp/s3-c1800-stripe-T110-diagnostic-20261009
```

The diagnostic runner validates source, binary, candidate and workload
hashes, factor-base and root-table identities, rank, scalar replay, and
exclusive one-target online phase accounting. It retains complete reports
and all statuses.
