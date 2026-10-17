# Restore the photographed-curve volcano campaign

**Status: blocked on unavailable experiment files.** This is a recovery
handoff, not a completed experiment or a publication of recovered evidence.
No solver, measurement, topology certificate, or diagram skill is included in
this change because the original files are unavailable in the current runtime.

The original workload remains **ten generators times ten paired random scalar
values per curve**, transported through the complete selected 2-isogeny
component. Rho, BSGS and actual index-calculus recovery, point/image orders,
setup/online/verification costs, failures, raw records, per-curve and per-problem
min/mean/max, and the multilevel heatmap remain required deliverables. This note
does not reduce that workload or reinterpret diagnostics as completed recovery.

## Exact source and provenance

The photographed shortlist is [crypto PR #2056](https://github.com/aburan28/crypto/pull/2056).
The selected row is `shortlist-deep-broad`, ICV1 slug
`icv1-fp57-t354478198-afcea726`, at pinned source commit
[`ef18058588077ef805c55c5dc25140b6bde2046b`](https://github.com/aburan28/crypto/blob/ef18058588077ef805c55c5dc25140b6bde2046b/research/curve_shortlist_20261016/curves.json).

| Parameter | Exact value |
| --- | --- |
| Field prime | `72057980207847001` |
| Source `a`, `b` | `18456521102804593`, `19958176276690439` |
| Curve order | `72057979853368804` |
| Prime subgroup order | `18014494963342201` |
| Cofactor | `4` |
| Trace | `354478198` |

The source row was fetched again and its parameter strings checked. That read
validates the source attribution, not the absent graph or measurements.
The historical native solver was pinned to
[`aburan28/cryptanalysis@b7352eec52676b7b82a220ba296661e8ba878915`](https://github.com/aburan28/cryptanalysis/tree/b7352eec52676b7b82a220ba296661e8ba878915).
Its historical binary and observer adapter have not been recovered here.

## Requirement status

The following values are **unverified prior checkpoint metadata**, not new
measurements and not admitted results. None can be independently replayed in
this runtime until the original archives are restored.

| Requirement | Prior reported extent | Current status |
| --- | --- | --- |
| Selected 2-isogeny component | 46 curves, 90 directed maps, five levels | reported complete; certificate unavailable |
| Matched fixtures | 100 problems/curve, ten generators and ten paired scalars | reported generated and transported; original bytes unavailable |
| Rho | 46 distinct problem/curve cases out of 4600 | partial; checkpoints unavailable |
| BSGS | 20 distinct problem/curve cases out of 4600 | partial; checkpoints unavailable |
| Full IC logarithm recovery | 0 out of 4600 | unresolved |
| IC relation-yield/regularity probes | 4600 target diagnostics reported | diagnostics only; raw records unavailable |
| Full measurements and heatmaps | parent report and paired BSGS continuation reported | not recoverable from the current filesystem |
| New solver executions in this handoff | 0 | no running solver job |

The small IC factor base was reported to yield zero random relations. That
bounded observation is not a proof of mathematical impossibility and does not
satisfy the requested full IC recovery. Restoring files also does not establish
a viable full-recovery IC configuration.

## Recover the original bytes

[recovery_manifest.json](recovery_manifest.json) records expected archive names,
byte counts, SHA-256 digests and sibling directory layout. These are prior
manifest values; the archive bytes and their durable remote location are
currently unavailable. An archive name or old local path is not a durable
storage receipt. No `store:` locator is invented here.

Place the two full archives in a recovery directory, then verify their sizes
against the manifest and run:

```sh
sha256sum -c archive_checksums.sha256
mkdir -p restored
# Inspect archive entries before extraction and use only the verified originals.
tar -tzf photo_curve_20261017_full.tar.gz
tar -tzf photo_campaign_20261017_full.tar.gz
tar -xzf photo_curve_20261017_full.tar.gz -C restored
tar -xzf photo_campaign_20261017_full.tar.gz -C restored
(cd restored/photo_curve_20261017 && sha256sum -c SHA256SUMS)
(cd restored/photo_campaign_20261017 && sha256sum -c SHA256SUMS)
```

Keep both directories as siblings. Do not regenerate random scalars and label
them as recovered fixtures, overwrite failed attempts, merge unrelated run
records, or infer measurements from missing cells. Restore and validate the
original checkpoint identity before executing its saved resume commands.
The old binaries used CPU-specific build flags: compare the host and compiler
and rebuild from preserved source when necessary, preserving the original
binary and source bindings and treating rebuild validation as a new receipt.

Once restored, put the archive in the repository's authorized immutable store
and commit its real content-addressed locator and readback verification. Do not
publish a guessed locator. The missing storage and separate downstream solver
obligations are tracked in the [problem record](../../problems/photo-volcano-checkpoints-unavailable.md).

## Validation and scope

See [VALIDATION.md](VALIDATION.md) for the exact scope of checks. Local artifact
absence is an observed infrastructure gap; the reason for that absence is
unknown. Existing canonical graphs and scoreboards are unchanged because no
new graphable measurement or recovered certificate is available. This PR does
not report timings, a speedup, an IC improvement, a Lean theorem, an external
solver ranking, or a completed experimental campaign. Historical PDFs and the
multilevel diagram skill must be restored rather than reconstructed as evidence.

## Fresh replacement authorized, 2026-10-17

A new experiment with fresh secure fixtures is now available in
[photo_volcano_restart_20261017](../photo_volcano_restart_20261017/README.md).
This authorizes continuing the full workload without restoring the old random
fixtures. The old archive gap and historical audit limitations remain; results
from the two experiments are not pooled. Public replacement inputs omit scalar
witnesses and raw randomness while publication approval is pending.
