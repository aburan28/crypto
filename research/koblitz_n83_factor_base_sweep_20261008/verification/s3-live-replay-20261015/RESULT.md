# Live N83 factor-base S3 archive replay

Status: `PASS_57_LIVE_ARCHIVES` on 2026-10-15. A fresh S3 listing found 72
objects under
`s3://crypto-autoresearcher/factor-bases/icv1/etc/koblitz_n83_factor_base_sweep_20261008/`:
57 compressed factor-base archives and 15 other JSON objects. All 54 archive
keys in the committed [pilot manifest](../../pilot-01/manifest.json) were
present. The three additional archives are pinned by the two v2 size-frontier
[K1182](../matched-k1182-k2048-900s-20261011/evidence-01/panel-k1182/manifest.json)
and [K2048](../matched-k1182-k2048-900s-20261011/evidence-01/panel-k2048/manifest.json)
manifests and the [v3 structured-orbit manifest](../structured-d12-orbit-20261010/manifest.json).

The 57 archives were downloaded to task-local scratch. The independent
[verifier](verify.py) read every compressed byte and full gzip stream,
checked each archive's compressed BLAKE3 against its content-addressed key,
and matched **all 57** compressed sizes, compressed BLAKE3 digests and plain
BLAKE3 digests to their four construction manifests. The total is
**109,656,159 compressed bytes** and **531,027,665 decompressed bytes**. The
retained [audit](audit.json)
contains each S3 URI, compressed SHA-256/BLAKE3, plain BLAKE3, sizes and pilot
membership, along with all four source-manifest hashes. Its saved result
passed a second byte-for-byte replay. The narrower original
[pilot-only audit](audit-pilot.json) remains byte-preserved. The
[raw listing](list.json) and 57-line [download log](download.stdout) are
retained; the download's stderr was empty. No S3 object was written or
replaced.

A separate [negative control](tamper-control.json) changed one byte in a
scratch copy of the smallest archive. The current verifier exited nonzero on
gzip's incorrect-length check and wrote no successful receipt; the
[original control](tamper-control-v1.json) also remains retained. This control
does not replace the positive replay of the original 57 archives.

This receipt establishes live storage and byte integrity of the stored
factor bases at the time of the read. It does not measure natural relation
yield, useful rank, linear algebra, individual-log work or total
index-calculus runtime. The local solver-stage ledger remains 23,100 seconds,
and the runtime-minimizing factor base remains unselected.

To reproduce from a fresh checkout with Python 3.13 and `blake3==1.0.8`,
download the archives to a new scratch directory, then check the retained
listing and receipt:

```sh
aws s3 cp s3://crypto-autoresearcher/factor-bases/icv1/etc/koblitz_n83_factor_base_sweep_20261008/ /private/tmp/n83-s3-replay-fresh/ --recursive --exclude '*' --include '*.jsonl.gz' --no-progress
python3 research/koblitz_n83_factor_base_sweep_20261008/verification/s3-live-replay-20261015/verify.py --download-root /private/tmp/n83-s3-replay-fresh/ --check
```
