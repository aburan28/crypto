# Frozen rotated-subspace evidence

This archive is the sole successful measured outcome of frozen PR #762 source
head `a6be357`. It was produced by:

```sh
/opt/homebrew/bin/python3.13 research/notes/ecc2k130/rotated_subspace_support_20260925/run.py --out /private/tmp/rotated-subspace-run-20260925-114521
```

`receipt.json` records Python/platform, UTC command windows, source/input
SHA-256, every raw file hash, success status and the independent verifier
command. Its SHA-256 is
`8b35cf4ea223b112d987abeed04a0ba8bf04e2f509ec0943eead87bdcdce9eac`.
`raw.tar.gz` contains the 42 exact `raw/` files (eight toy summaries,
eight 2,003-row target censuses, eight factor lists, the toy aggregate,
five density summaries, five 256-row covariance files, five 16,384-row
density files, the density aggregate and independent `verify_report.json`).
It has 2,680,225 bytes and SHA-256
`fad4b4cc48f980349d1ef472d18327368dd0f93f3d1284fc8bbdc20f0d3086d7`.

From repository root, validate without changing the frozen evidence:

```sh
python3 research/notes/ecc2k130/rotated_subspace_support_20260925/ci_replay.py --evidence research/notes/ecc2k130/rotated_subspace_support_20260925/evidence
```

For manual inspection, extract with
`tar -xzf research/notes/ecc2k130/rotated_subspace_support_20260925/evidence/raw.tar.gz -C /tmp`.
The replay script checks the archive SHA, its full raw-file hash map, every
n=13 full-point sum and witness/negative decision, and every n=131 density
sample and covariance row against a separately implemented field/group law.
