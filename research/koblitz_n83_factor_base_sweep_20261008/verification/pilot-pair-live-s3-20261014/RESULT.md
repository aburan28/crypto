# Live S3 replay of a matched K600 pilot-pair contrast

On 2026-10-14, three public primary-curve signed-Frobenius factor-base
archives were downloaded afresh from the study's S3 prefix. The fixed base is
K600 `public_x_sequential`, seed `2026100802`. The two K600 pools are
`public_x_gray_prefix` and `public_x_hash`, both seed `2026100801`. Their
exact S3 URIs, compressed-byte SHA-256 and BLAKE3 digests, manifest binding,
and independently computed orbit-set hashes are in [audit.json](audit.json).
[heads.json](heads.json) preserves the live S3 HEAD observation at
2026-10-14 16:22:26 UTC. HEAD length and ETag alone do not authenticate object
content; the downloaded compressed bytes matched each content-addressed
manifest row before decompression.

The [replay program](audit.py) checks all 99,600 point rows in each archive,
including column, phase, sign, coefficient, point range, uniqueness, and the
166-row signed-orbit structure. It canonicalizes each orbit from its *complete
point set*, checks that a shared orbit has identical points in both archives,
and independently recomputes point-set intersection. It then compares the
result with the earlier [all-pair support screen](../pilot-pair-support-screen-20261014/RESULT.md).
This is an archive and overlap audit; the earlier full-panel audit covers
group-law replay, and this audit does not make a new group-law or solver claim.

| Pool with the same K600 sequential base | Shared signed orbits | Residual signed orbits | Residual points | Literal header representatives shared |
| --- | ---: | ---: | ---: | ---: |
| Gray-prefix K600 | 503 | 97 | 16,102 | 489 |
| Hash K600 | 0 | 600 | 99,600 | 0 |

The 14-orbit gap between the Gray pool's literal representative overlap and
full-orbit overlap shows why comparing only header representatives would be
wrong. These are exact set facts for the three retained archives. The
conditional uniform-query 3+2 support ceilings remain about 0.0088301243
for Gray and 0.3378339093 for hash; the earlier specified K600-row query
floors remain 33,975 and 889. Those are counting bounds, not observed
relation probabilities or solver throughput.

The existing `run-cross-pilot` worker accepts both pinned pairs. A next
bounded comparison can hold this sequential base, the pinned public corpus,
the worker's deterministic first-trial coefficients, the 3+2 shape, conflict
cap, model cap, memory limit, and host conditions fixed while changing only
the pool. The two pool policies still differ, so a timing difference would
describe these exact configurations rather than isolate overlap as a cause.
Both arms require separately charged cold import and search receipts,
independent candidate group replay, and a later multi-query relation/rank,
linear-algebra, individual-log and matched-rho campaign before selecting a
minimum-runtime factor base. This audit launched **no N83 solver search** and
does not increase the conservative 23,100-second retained search ledger.
The saved `selected_best_total_runtime` remains null.

To replay, download the three `objects[*].s3_uri` values in `audit.json` to a
scratch directory under their `<compressed_blake3>.jsonl.gz` names, then run:

```sh
uv run --with blake3 python3 research/koblitz_n83_factor_base_sweep_20261008/verification/pilot-pair-live-s3-20261014/audit.py --archive-dir SCRATCH_DIR
uv run --with blake3 python3 -O research/koblitz_n83_factor_base_sweep_20261008/verification/pilot-pair-live-s3-20261014/audit.py --archive-dir SCRATCH_DIR
```

Both saved-receipt replays passed. In a negative control, one byte of the
scratch Gray archive was altered; optimized Python exited nonzero with
`compressed BLAKE3 mismatch`. The original S3 downloads and the altered
scratch copy remain outside the Git worktree. The exact receipt source hash
and parent input hashes are pinned in `audit.json`.
