# Fresh F4/F5/SAT complete-solve qualification

Status: **registered, not yet measured**. The scientific question, schedule,
resources, success/stop rules and accounting are frozen in [PROTOCOL.md](PROTOCOL.md).
The first registration was [operationally censored](../generic-backend-qualification/RESULT.md);
it established no solver result and must never be redispatched. This panel
uses seed `2026092902`, one process for each of 25 distinct public targets,
250 paired trial slots, and the same five toy curve cells and source-bound
algorithm arms. It does not use the sealed improvement confirmation sets.

The first run's entire possible 25-point fixture schedule is in
[lost-campaign-exposures.json](lost-campaign-exposures.json). An independent full-prepare replay matches that corpus byte-for-byte on all 25 points; see [INDEPENDENT-RECONSTRUCTION.md](INDEPENDENT-RECONSTRUCTION.md). It was reproduced
from the pinned worker and sealed historical exclusions using
[reconstruct_exposures.py](reconstruct_exposures.py); five retained prior
fixtures matched exactly and a second local reconstruction matched the checked
file's SHA-256 byte for byte. The fresh runner passes that corpus to the
tournament's exact curve/point exclusion checker before making new targets.

The workflow is
[`ic-generic-backend-qualification-v2.yml`](../../../../.github/workflows/ic-generic-backend-qualification-v2.yml).
Its PR jobs run the scientific controls and prove that a synthetic interrupted
trial can be packed and uploaded. The campaign job is dispatch-only on main.
The measured step has a 300-minute cap inside a 360-minute job; complete or
partial evidence is packed into a single `.tar.zst` file before the job cap.
Watch the live per-trial log and retain the output archive and manifest. If the
campaign is incomplete, report it as operationally censored with unknown
family qualification and competitive costs; do not retry this seed.

After the campaign, verify the archive SHA-256, extract it, inspect the frozen
`tournament/evaluator/tournament.py verify` receipt, `natural-yield.json`,
`family-gate.json`, target history, source/build record, and all failure rows.
Publish the bundle's durable artifact link, exact hash and results in a new
evidence PR. A complete, independently verified F4/F5 and SAT arm is the
qualification goal; the development panel is not a held-out or global speed
claim.
