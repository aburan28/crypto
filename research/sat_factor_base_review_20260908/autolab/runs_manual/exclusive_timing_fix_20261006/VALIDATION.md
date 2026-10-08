# Exclusive online timing correction: validation

The direct rank fixture's `collection_ms` timer encloses final linear solving
and scalar replay. The previous online total added both again. It also left
target-specific allocation outside its interval and assigned in-collector
target generation and group verification to PDP. The revised producer charges
target-specific setup, assigns the measured nested intervals to the five
online phases, and computes the total from the inclusive outer collection
timer exactly once. The autolab claim builder now rejects overlapping legacy
rows. This is an accounting correction; this smoke run is **not** an isolated
CPU benchmark or a new speedup claim.

The source parent is `ea8892a530505a07ca1274f614cdba3b16f96c5c` on
`cursor/ic-boundary-experiments-d111`. The build used the ignored local
`Cargo.lock`, archived here as `producer-Cargo.lock`, with SHA-256
`f6e3447b98a81f515b4a73deb59bdc0fe39dae6046be0fe9a0b111d25e76e7b6`.
The revised `examples/koblitz_rank_fixture.rs` SHA-256 is
`a34aaa5177ea8b22577782979ee341a8196a5584c3d47300af1f5971f28bd32d`;
the release example binary SHA-256 is
`a9e2db10d6ac4fdb7f81caf88cbb0580e21c3de32cab7473458e502a49f2359d`.

Validation commands from the repository root:

```sh
CARGO_TARGET_DIR=/private/tmp/crypto-cold-source-20261006/target cargo test --locked --offline --release --example koblitz_rank_fixture
python3 -m unittest discover -s research/sat_factor_base_review_20260908/autolab -p test_boundary_autolab.py
python3 research/sat_factor_base_review_20260908/autolab/runs_manual/exclusive_timing_fix_20261006/check_timing_evidence.py
```

The Rust example suite passed 9/9; the Python autolab suite passed 15/15.
The exact n=13 smoke command was:

```sh
KIC_SHARED_FACTOR_LOG_PRECOMPUTATION=1 /private/tmp/crypto-cold-source-20261006/target/release/examples/koblitz_rank_fixture 13 0 1 2 123 signed_expanded independent pair_pair_16 2
```

Its raw JSONL SHA-256 is
`09caf64810a7d3d40200d611bfaf28c1e0c606e5ded27e017da45c0e096ba417`.
The checked-in `smoke-n13.jsonl.gz` is a deterministic gzip of those bytes.
The smoke yielded a verified precomputation fixture and a distinct verified
one-relation online target. The five online phases sum to its reported
0.2325 ms, and the batch charged sum matches the two fixture intervals.

`legacy-n41-online-extract.json` is the exact relevant field subset from
[cryptanalysis PR #353](https://github.com/aburan28/cryptanalysis/pull/353)
raw `n41-direct.stdout.jsonl` (full-file SHA-256
`24ee32f71b67d306b12500db18f66f257bf9382354760e776af2a58a7ec10871`).
Its old phases sum to the published 6.544 ms, but the exclusive validator
rejects the row: 1.930167 ms of final solve and replay was counted twice.
Removing that overlap gives 4.613833 ms by arithmetic. This arithmetic does
not turn the old run into a controlled speedup measurement: its setup interval
also omitted target-specific allocation, and the host was not isolated.
