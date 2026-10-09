# N53 C1800 root index with a C1200 target-state prefix

The C1800 index contains the C1200 left/right/relative state set. This
candidate searches those embedded states in the exact C1200 point-seeded
order while querying the complete C1800 root table. If the prefix yields no
verified relation, it traverses all C1800 states. The fallback preserves
complete search coverage. Preparation still builds a full-rank C1800
factor-base log system before the public target enters the online clock.

[`SCREEN_FREEZE.json`](SCREEN_FREEZE.json) fixes the first fresh paired screen:
eight previously unused N53 public targets T109–T116, one target per run, two
repetitions per arm, ABBA/BAAB order, 14 preparation workers, one target
worker, a 180-second per-run cap, and the exact release binary and source
hashes. Both arms use the same binary. C1200 uses its native state order;
C1800 uses the prefix schedule. The source was committed before the freeze.
Candidate manifests in [`identities/candidates`](identities/candidates) record
the exact factor bases, search policies, and implementation hashes.

The predeclared screening gate is a geometric mean C1200/C1800 target S3-call
ratio of at least 2, computed across the eight points. If it passes and every
scalar independently replays, T117–T124 remain reserved for a separate
confirmation. Online wall time is recorded with the full target-dependent
phase sum, but local host wall-time ratios remain exploratory without an
exclusive CPU partition receipt. A controlled wall-time claim also requires
the isolated host and noise gates in
[`docs/ISOLATED_BENCHMARKS.md`](../../docs/ISOLATED_BENCHMARKS.md).

The frozen runner validates each candidate and workload hash, public point,
phase sum, full-rank status, root-table size, scalar, factor-base digest, and
search policy. It retains all failed, timed-out, and verified run rows in an
append-only ledger. Run it with:

```sh
python3 experiments/koblitz-s3-c1800-prefix-20261009/run_paired.py \
  --binary target/release/examples/s3_partial_prepared_fullrank_c1200prefix_targetspan \
  --output /private/tmp/s3-c1800-prefix-screen-20261009
```

`prepare_freeze.py` reproduces the static manifests and freeze from the
frozen source commit and release binary. Running it after a source or binary
change produces a different candidate identity and requires a new protocol.
