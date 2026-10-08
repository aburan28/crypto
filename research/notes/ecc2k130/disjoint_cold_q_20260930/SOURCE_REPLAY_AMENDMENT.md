# Restoring the exact compact source freeze

The [disjoint-Q protocol](PROTOCOL.md) reuses the twenty-file `compact_s3_prefilter_20260930/FROZEN.json` source freeze. On the 2026-09-30 main checkout, `compact_frozen_source_replay_20260929/materialize.py` stopped at `examples/koblitz_orbit_dlp_s3_batch.rs`: its pinned SHA-256 `702a0a05709bc14bc10bafbb08edb2a2f5f86794e970d84c317da3bc65cdcf38` differed from main, and `SNAPSHOTS.json` did not contain it. Four other drifted pinned files also lacked snapshots. This was a source-replay failure, not an experiment outcome; no new Q had been generated and no arm had been timed.

PR #1095 adds the following exact bytes to the existing gzip snapshot collection. Each file was read from the named historical Git commit, checked against the **unchanged** `FROZEN.json` SHA-256, compressed with `gzip.compress(data, mtime=0)`, and recorded with its raw length and gzip SHA-256 in `SNAPSHOTS.json`.

| Pinned path | Historical commit | Raw SHA-256 |
|:--|:--|:--|
| `examples/koblitz_orbit_dlp_s3_batch.rs` | `2c207b41` | `702a0a05709bc14bc10bafbb08edb2a2f5f86794e970d84c317da3bc65cdcf38` |
| `research/notes/ecc2k130/compact_s3_prefilter_20260930/PROTOCOL.md` | `52270f7e` | `74bd16b4bb8fd21b2457e2e7690ce88f548bc2b5ea40edbeced492ff81e87bdf` |
| `research/notes/ecc2k130/compact_s3_prefilter_20260930/prepare.py` | `2c207b41` | `eda8080fa4f508c08130129479e0ad26f3235e1b628e27b5bbfc9db52ce85362` |
| `research/notes/ecc2k130/compact_s3_prefilter_20260930/run_panel.py` | `2c207b41` | `98cc9a6d20b50742a973ed7e35a1235be753b7fcf2f05390c821f0a66c073ab9` |
| `research/notes/ecc2k130/compact_s3_prefilter_20260930/verify_panel.py` | `059c74b3` | `9c6acec9cc5f0129f71392d8943aedcf37e96f7bd95d6e439fa4842f0dc355c6` |

After this addition, the materializer replayed all twenty pinned files successfully from the current sparse checkout: nine were already byte-identical and eleven came from SHA-checked snapshots. The four fail-closed materializer unit tests passed. No file in the historical source freeze was repinned or edited; the amendment only makes its existing bytes available for independent rebuilds.
