# C1800 prefix-256 striped target search

The frozen T110 control verified both scalars. C1200 used **1,822** target
S3 calls, the earlier C1800 prefix used **964**, and the new striped C1800
path used **542**. The old C1800/striped call ratio is **1.779×** on this
previously used point. It demonstrates the new branch executes and can
shorten one search; fresh-target yield remains to be measured.

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
[`local/raw_diagnostic2.jsonl.gz`](local/raw_diagnostic2.jsonl.gz) and
[`local/raw_reports.tar.gz`](local/raw_reports.tar.gz) retain both complete
runs. [`analyze_diagnostic.py`](analyze_diagnostic.py) checks report digests,
source and candidate identities, schedule fields, scalar equality, and the
prior T110 reference before producing
[`diagnostic_analysis.json`](diagnostic_analysis.json). The decompressed
ledger SHA-256 is
`0bfd813a1cabc198606c24e3a3a19bf49cf2fa4d273fbbd95161e9993109870b`;
the report archive SHA-256 is
`5e41ea0c7063ce45c4b81b65c741b76111ef387778c4dccd3a29403fca4eeb26`.
