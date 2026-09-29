# Reconstruction of potentially exposed points from seed 2026092901

The censored Actions run
[36532455386](https://github.com/aburan28/crypto/actions/runs/36532455386)
published no artifact. The runner may still have generated every fixture for
panel seed `2026092901` before measurement was canceled. Those points are
treated as exposed.

## Method

At frozen execution checkout
`765c3c5f19032bd852163805f257c56babef2040` (PR #920 head), replay the
registered restore / reference-registry / sealed-archive / controlled generic
build / `tournament.py prepare --qualification` path with seed `2026092901`,
identical cells, exclusions and resource arguments as
`run_generic_backend_qualification.py`. Stop after prepare; do not run or
verify trials.

## Digests

| Artifact | SHA-256 |
| --- | --- |
| Censored panel.json | `83c640a03b4239b918851f6f1b8450e2fe27dc710fb99f306e3481a99d8875cf` |
| Reconstructed tournament `fixtures.json` | `78055bac6b560aa67ddf89c768e3d5173bc1c5115424409ad195f7354cd7aab5` |
| Exclusion export [censored-2026092901-fixtures.json](censored-2026092901-fixtures.json) | `30f496a561aa0814adbf2822d94ec4e98c808e1a3fa9219d592b182521438004` |

Machine-readable receipt: [reconstruction-receipt.json](reconstruction-receipt.json).

## Counts

25 public points: 5 A/A, 5 smoke, 15 development. Every later campaign must
pass this file to `--exposed-fixtures` (or an equivalent sealed supplemental
exposure) before sampling a new target. Redispatch of seed `2026092901` remains
forbidden.
