# Generic F4/F5 subspace smoke qualification (third registration)

Status: **planning / frozen protocol only**. No measurement is authorized by
this document until a separate evidence PR lands the runner, workflow, sealed
panel digest lock, and a green feasibility preflight. Seed **2026093001** is
reserved here and must not be reused if this registration is abandoned.

## Why this registration exists

The first two complete-solve registrations are operationally censored:

| Seed | Run | Outcome |
| --- | --- | --- |
| `2026092901` | [36532455386](https://github.com/aburan28/crypto/actions/runs/36532455386) | Job-cap cancel; no bundle ([RESULT](../generic-backend-qualification/RESULT.md)) |
| `2026092902` | [36580669479](https://github.com/aburan28/crypto/actions/runs/36580669479) | 300-minute measure timeout; pack failed; no campaign artifact ([RESULT](../generic-backend-qualification-v2/RESULT.md)) |

Neither is a mathematical F4/F5/SAT failure. Ambient `subgroup_orbits` F4/F5
layouts on the v2 panel are also statically encoder-unsupported at `m=3`
([STATIC-FEASIBILITY](../generic-backend-qualification-v2/STATIC-FEASIBILITY.md)).
Disclosed-point work since then shows:

- Dimension-6 `standard_subspace` clears `MAX_VARS=64` and dispatches MatrixF4/F5
  ([d6 pilot](../generic-f4-subspace-pilot/RESULT.md)).
- Under a recovery-capable budget, disclosed n17a1 F5 completed a verified
  one-target solve; the combined F4/F5-and-SAT gate still failed
  ([recovery pilot](../generic-backend-recovery-pilot/RESULTS.md)).

This registration therefore asks a **narrower, encoder-feasible question** and
sizes the schedule to the 300-minute measure / 35-minute pack envelope.

## Question and success rule

On fresh public points that exclude every prior exposure, can at least one
source-bound generic F4/F5 arm that uses `standard_subspace` dimension 6
complete independently verified one-target IC solves on **every** A/A and smoke
point in the five-cell panel, with natural-query evidence and exclusive costs?

Success requires, for that arm, all scheduled A/A and smoke jobs verified: no
censored query audit, checked usable base and folded columns, observed PDP
dispatch into MatrixF4 or MatrixF5 with `unsupported: false`, stored relation
matrix and full rank, final relation LA, target descent, scalar replay,
complete profile/native phase accounting, and matching source/build identity.
The family gate for this registration is F4/F5-only; SAT is **out of scope**
here and remains on the separate source-complete SAT track.

The accepted incumbent and the two rho roles stay paired on the same public
point. Failures, timeouts, OOMs and unverified work retain rows with null
competitive totals. An incomplete schedule is operationally censored, not a
negative solver result. Do not widen limits, substitute points, or retry seed
`2026093001` after dispatch.

This is a development smoke qualification: no held-out confirmation, no
familywise promotion, no ECC2K-130 claim. Comparison kind is
`factor-base-policy` because the algebraic arms change the factor-base recipe
relative to the orbit-base incumbent.

## Exclusions (frozen before sampling)

Exclude, by exact curve ID and point:

1. The original improvement target history and all three sealed improvement
   archives.
2. The seven supplemental fixture corpora listed in the accepted round-three
   panel.
3. First-run censored exposures
   [lost-campaign-exposures.json](../generic-backend-qualification-v2/lost-campaign-exposures.json)
   (SHA-256 `a728677b199eac02800d8338204d5306f391ec5da757c910bec1e51955fe7b41`).
4. Second-run censored exposures
   [lost-v2-campaign-exposures.json](../generic-backend-qualification-v2/lost-v2-campaign-exposures.json)
   (SHA-256 `0cc792cceb8c7190a533e6f4e665e8486911153ad57609ac8d1243af282b8e54`),
   once that file is merged from the exposure-census PR; until then this
   planning document cites the sealed hash only.

Never redispatch seeds `2026092901` or `2026092902`. Do not reuse disclosed
recovery-pilot or d6-pilot points as fresh competitive targets.

## Frozen panel (intent)

[panel.json](panel.json) freezes the scientific intent. Before any dispatch,
recompute and lock its byte SHA-256 in the runner registration check, and run

```text
python3.12 research/ic_candidate_tournament_20260915/generic_solver_feasibility.py \
  --repo <pinned-worker-checkout> --panel panel.json --require-pass
```

against the exact worker source that will measure.

| Field | Value |
| --- | --- |
| Seed | `2026093001` |
| Cells | `n17a1,n19a0,n23a0,n23a1,n31a0` (holdout `n29a1` unused) |
| Stages | `aa`, `smoke` only (no development stage) |
| Repetitions | 1 process per point |
| Algebraic factor base | `standard_subspace`, dimension 6 |
| Algebraic arms | `generic_f4_subspace_dense`, `generic_f5_subspace_dense` |
| References | qualified pairinv incumbent, `ic_online`, cold rho, online rho |
| Continuity | one generic `pair_table` dense arm on the **same** subspace base |
| Algebraic budgets | `summands=3`, `groebner_degree=3`, `node_budget=4096`, `max_trials=256`, `batch_trials=8` |
| Child envelope | 300 s, 8 GiB, 1 CPU, `RAYON_NUM_THREADS=1` |
| Job envelope | measure ≤ 180 min; pack ≤ 35 min; job ≤ 240 min |
| Worker pin | `765c3c5f19032bd852163805f257c56babef2040` (same objects as v2) |

| Stage | Distinct points | Arms (illustrative) | Trial slots |
| --- | ---: | ---: | ---: |
| A/A | 5 | incumbent + aa_control | 10 |
| Smoke | 5 | 2 algebraic + 1 pair continuity + 2 IC refs + 2 rho | 35 |
| Total | 10 | — | 45 |

Cap scheduled slots at 60. The reduced wall budget (180 minutes measure)
reflects the disclosed F5 recovery cost (~274 s process on n17a1) and the
prior 250-slot overrun. If smoke cannot finish inside 180 minutes under these
caps, the campaign is censored; do not extend the seed.

Ambient `subgroup_orbits` F4/F5/SAT arms from v1/v2 are **not** registered
here. SAT complete-solve qualification stays on its own source-bound track.

## Cost and accounting

Primary metric: verified one-target native online wall time after reusable
preparation through scalar replay, paired against matched rho on the same
point. Report `S = total Ir / sqrt(r)` and ratios to rho only for verified
complete pairs. Missing phases stay null. Natural-query audits must distinguish
UNSAT, bounded-incomplete, timeout and absent report. Preserve zero-yield
cells. Class every change as advance / engineering / relabelling / accounting
per repository §3; a smoke pass is not an end-to-end ECC2K-130 speedup.

## Run and retention

Dispatch once from main after the runner, workflow and panel digest land.
Use the amended packer that retains `ARCHIVE_READ_UNVERIFIED` archives and
uploads the `.tar.zst`, manifest and tar stderr. Quiesce or terminate measure
orphans before tar. After the run: verify archive hash, frozen
`tournament.py verify`, natural-yield audit, F4/F5 family gate, and independent
online-pair recalculation. Publish durable artifact link and hash in an
evidence PR. Never retry this seed.
