# Generic F4/F5 subspace smoke qualification (v3)

Status: **planning**. Scientific question, exclusions, schedule and accounting
are frozen in [PROTOCOL.md](PROTOCOL.md). The intent panel is
[panel.json](panel.json) (`REGISTERED_PLANNING`). No runner, workflow or
measurement is authorized yet.

## Dependencies

| Gate | State |
| --- | --- |
| Censored v1 exposures | Sealed (`lost-campaign-exposures.json`) |
| Censored v2 exposures | Sealed hash `0cc792cceb8c7190a533e6f4e665e8486911153ad57609ac8d1243af282b8e54` (`lost-v2-campaign-exposures.json`; merge via exposure-census PR) |
| Encoder-feasible F4/F5 layout | `standard_subspace` dimension 6 (d6 + recovery pilots) |
| Disclosed F5 recovery | Yes on n17a1; not a fresh-target qualification |
| SAT family | Out of scope for this panel |

## Next implementation steps (follow-on PRs)

1. Land the v2 exposure corpus on `main` if not already merged.
2. Pin panel byte SHA-256 after any final panel edit; add
   `generic_solver_feasibility.py --require-pass` CI control.
3. Add `run_generic_backend_qualification_v3.py` and a dispatch-only workflow
   with 180/35/240 minute measure/pack/job caps and the amended packer.
4. Dispatch seed `2026093001` once from main; never retry.

Do not redispatch seeds `2026092901` or `2026092902`.
