# Generic F4/F5 subspace smoke qualification (v3)

Status: **planning, registration locked**. Scientific question, exclusions,
schedule and accounting are frozen in [PROTOCOL.md](PROTOCOL.md). The intent
panel is [panel.json](panel.json) (`REGISTERED_PLANNING`, SHA-256
`df92d5507785446a2a5b333bd7a04776781de4f54c63c5ea3ffa921bc99ad2e4`).
`run_generic_backend_qualification_v3.py` checks that lock and the static
encoder layout, then refuses to sample points or start a job.
`dispatch_authorized` is false. Seed `2026093001` has not been run.

## Dependencies

| Gate | State |
| --- | --- |
| Censored v1 exposures | Sealed (`lost-campaign-exposures.json`) |
| Censored v2 exposures | Sealed hash `0cc792cceb8c7190a533e6f4e665e8486911153ad57609ac8d1243af282b8e54` (`lost-v2-campaign-exposures.json`; merge via exposure-census PR) |
| Encoder-feasible F4/F5 layout | `standard_subspace` dimension 6 (d6 + recovery pilots) |
| Disclosed F5 recovery | Yes on n17a1; not a fresh-target qualification |
| SAT family | Out of scope for this panel |

## Next implementation steps

1. Land the v2 exposure corpus on `main` (`lost-v2-campaign-exposures.json`,
   SHA-256 `0cc792cceb8c7190a533e6f4e665e8486911153ad57609ac8d1243af282b8e54`).
   The checker reports this file missing as `dispatch_block` and still exits
   without measuring.
2. A later commit may add one dispatch-only workflow job with the 180/35/240
   minute measure/pack/job caps. That commit is the only place
   `dispatch_authorized` may become true, and only after the corpus above is
   on `main`.
3. Dispatch seed `2026093001` once; never retry.

Do not redispatch seeds `2026092901` or `2026092902`.
