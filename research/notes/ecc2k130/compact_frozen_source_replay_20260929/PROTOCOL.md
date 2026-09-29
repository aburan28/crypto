# Compact-orbit frozen-source replay

The compact S3-batch and swap-quotient evaluations (`compact_s3_batch_20260929`,
`compact_swap_quotient_20260929`) pin, in their `FROZEN.json` `source_sha256`,
the bytes of shared files that later PRs legitimately edit:
`src/cryptanalysis/koblitz_index_calculus.rs`, `Cargo.toml`, `src/lib.rs`, the
workflow files, and others. Their `run_panel.py` and `verify_panel.py` hash
the checkout they run from, so the first PR to touch any of those files broke
the evaluation, even though nothing about the frozen evaluation had changed.

Following the historical-snapshot replay used for the rotated solver-admission
receipt (PR #801) and the n53 rank-rotation replay (PR #823), the workflows
now build and run from a materialized copy of the checkout instead of
re-pinning hashes:

1. `materialize.py --freeze <FROZEN.json> --dest <dir>` copies the sparse
   checkout (without `.git` and `target`) to `<dir>`.
2. For each pinned file whose live bytes differ from the freeze, it writes the
   committed snapshot `historical_snapshots/<sha256>.gz` instead. The snapshot
   is accepted only if `SNAPSHOTS.json` lists it under the same path, its gzip
   and decompressed SHA-256 and length match, and the decompressed hash equals
   the frozen digest. A drifted file without a matching snapshot fails closed.
3. Every pinned file in `<dir>` is re-hashed against the freeze.

The frozen `run_panel.py` and `verify_panel.py` then run unmodified from
`<dir>`, so their own source checks pass against the frozen bytes. `GIT_DIR`
points at the real checkout, so `host.json` still records the commit under
test. `FROZEN.json`, `SOURCE_FROZEN.json`, fixtures, and all recorded
measurements are unchanged.

`SNAPSHOTS.json` snapshots every pinned file outside the two panel directories,
taken from main at `source_commit`, where each matched its frozen digest.
