# Current S3 content audit of all retained N83 bases

The [write-once audit receipt](audit.json) reports `PASS` for all 57 retained
factor-base objects: 54 primary pilot objects, the K=1,182 and K=2,048 size
extensions, and the dimension-12 structured-orbit object. It streamed
109,656,159 compressed bytes from their existing S3 keys. For each object,
the observed byte count and BLAKE3 digest match the committed manifest. The
receipt contains each key, its source manifest, both expected and observed
values, the four manifest SHA-256s, the auditor SHA-256, and the UTC check
time. No object was saved to disk, and no N83 solver search was launched.

The [auditor](audit.py) rejects an unexpected S3 prefix, changed inventory,
nonzero AWS exit, short or long stream, or digest mismatch. It uses at most
four concurrent transfers and a 120-second limit per transfer. The run used
AWS CLI 2.27.12 and Python `blake3` 1.0.11. Recheck against the live keys from
the repository root with:

```sh
UV_CACHE_DIR=/private/tmp/n83-uv-cache uv run --no-project --with blake3 python3 research/koblitz_n83_factor_base_sweep_20261008/verification/s3-stream-integrity-20261014/audit.py --check
```

This is a **current availability and compressed-content** check, timestamped
2026-10-14 05:21 UTC. It does not retroactively prove the original upload
time or replace the earlier construction, generic replay, and S3 round-trip
receipts. It measures no natural relation yield, useful rank, target solve, or
complete cold index-calculus runtime, so it does not select a factor base.
