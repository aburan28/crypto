# Reconstruction of seed 2026092902 exposures

The second registration’s measured dispatch
([Actions run 36580669479](https://github.com/aburan28/crypto/actions/runs/36580669479))
produced no auditable tournament bundle. Every public point that registration
could have generated must still be excluded from later panels.

[reconstruct_v2_exposures.py](reconstruct_v2_exposures.py) rebuilds that
25-point schedule with the pinned worker
`765c3c5f19032bd852163805f257c56babef2040`, the same sealed improvement
archives and supplemental `EXPOSED` corpora used by
`run_generic_backend_qualification_v2.py`, plus the first-run corpus
[lost-campaign-exposures.json](lost-campaign-exposures.json)
(SHA-256 `a728677b199eac02800d8338204d5306f391ec5da757c910bec1e51955fe7b41`).
It does not redispatch the campaign or invent new seeds.

The sealed export is
[lost-v2-campaign-exposures.json](lost-v2-campaign-exposures.json)
(SHA-256 `0cc792cceb8c7190a533e6f4e665e8486911153ad57609ac8d1243af282b8e54`).
Machine-readable digests are in
[v2-exposure-reconstruction-receipt.json](v2-exposure-reconstruction-receipt.json).

| Check | Result |
| --- | --- |
| Accepted distinct points | 25 (aa 5, smoke 5, development 15) |
| Overlap with first-run accepted points | 0 |
| Second local reconstruction | byte-identical to the sealed export |
| Full `tournament.py prepare` replay | all 25 accepted fixtures match; export in [independent-v2-prepare-replay-fixtures.json](independent-v2-prepare-replay-fixtures.json) |
| Seed / run id | `2026092902` / `36580669479` |

A later competitive registration must pass this file to the tournament’s
`--exposed-fixtures` path (together with the three sealed rounds, seven
supplemental corpora, and the first-run lost corpus) before sampling new
targets. Never redispatch seed `2026092902`.
