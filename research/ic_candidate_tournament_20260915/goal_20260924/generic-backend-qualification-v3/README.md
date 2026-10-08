# Generic F4/F5 subspace smoke qualification (v3)

Status: **measured, negative, closed**. The one dispatch finished with
`NEGATIVE_FAMILY_QUALIFICATION`. Scientific question, exclusions, schedule and
accounting stay frozen in [PROTOCOL.md](PROTOCOL.md); the outcome is
[RESULT.md](RESULT.md).
The intent panel is [panel.json](panel.json) (`REGISTERED_PLANNING`, SHA-256
`df92d5507785446a2a5b333bd7a04776781de4f54c63c5ea3ffa921bc99ad2e4`).
`run_generic_backend_qualification_v3.py` locks that panel, requires both
exposure corpora, and authorizes one campaign path for seed `2026093001`.
That path has now executed once. Measurement is the retained negative result,
not `not_run`.

## Dependencies

| Gate | State |
| --- | --- |
| Censored v1 exposures | Sealed (`lost-campaign-exposures.json`) |
| Censored v2 exposures | On `main` (`lost-v2-campaign-exposures.json`, SHA-256 `0cc792cceb8c7190a533e6f4e665e8486911153ad57609ac8d1243af282b8e54`) |
| Encoder-feasible F4/F5 layout | `standard_subspace` dimension 6 (d6 + recovery pilots) |
| Disclosed F5 recovery | Yes on n17a1; not a fresh-target qualification |
| SAT family | Out of scope for this panel |
| Smoke schedule support | `tournament.py --qualification-schedule smoke` |
| Campaign job | `workflow_dispatch` only; 180/35/240 minute caps |

## Closeout

The campaign job is `if: false`. Seed `2026093001` must not be dispatched again.
The ten generated points in [retained/fixtures.json](retained/fixtures.json)
join the exclusion set for any later panel.

Do not redispatch seeds `2026092901` or `2026092902`.
