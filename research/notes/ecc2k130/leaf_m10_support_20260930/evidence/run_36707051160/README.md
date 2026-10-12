# Archived leaf m10 census, Actions run 36707051160

This directory is the hosted artifact `leaf-m10-support-36707051160`
(8,223,622 bytes, digest
`sha256:cfeb16421b261a5f7c9fc36bfb3a60349cd1f5a6357c1dd49e1874d3303d3344`)
plus the second-host replay files. `result.json` and `replay.json` are the
hosted producer and verifier receipts. `second_host_replay.json` and
`HOST.json` record the independent Linux replay. `SHA256SUMS` covers every
other file in this directory.

Do not regenerate these rows. A later check replays them with

```sh
python3 research/notes/ecc2k130/leaf_m10_support_20260930/verify.py \
  --evidence research/notes/ecc2k130/leaf_m10_support_20260930/evidence/run_36707051160 \
  --out /tmp/leaf-m10-replay-receipt.json
```

The verifier refuses to overwrite an existing receipt, so the output path
must be new. That command checks the archived rows; it is not a second census.
