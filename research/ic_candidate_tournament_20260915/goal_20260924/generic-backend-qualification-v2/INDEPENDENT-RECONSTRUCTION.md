# Independent prepare-path reconstruction of seed 2026092901

The sealed exclusion corpus for the second registration is
[lost-campaign-exposures.json](lost-campaign-exposures.json)
(SHA-256 `a728677b199eac02800d8338204d5306f391ec5da757c910bec1e51955fe7b41`),
produced by [reconstruct_exposures.py](reconstruct_exposures.py) and merged in
PR #946.

A second reconstruction replayed the frozen campaign's full
`run_generic_backend_qualification.py` restore / build /
`tournament.py prepare --qualification --seed 2026092901` path at checkout
`765c3c5f19032bd852163805f257c56babef2040` and stopped before measurement.
All 25 accepted public points match the sealed corpus exactly (same curve
fields, target seeds and targets). The independent export is
[independent-prepare-replay-fixtures.json](independent-prepare-replay-fixtures.json);
machine-readable digests are in
[independent-reconstruction-receipt.json](independent-reconstruction-receipt.json).

This note does not change the registered v2 panel, schedule or runner. It only
corroborates that the first campaign's potentially generated points are the
ones already excluded.
