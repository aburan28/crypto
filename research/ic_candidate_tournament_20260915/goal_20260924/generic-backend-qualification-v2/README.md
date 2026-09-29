# Generic backend qualification v2

Status: **REGISTERED_MEASUREMENT_PENDING.** Planning registration only; the
versioned runner and workflow are a follow-on PR. Do not dispatch measurement
from this directory until that implementation lands.

The first registered campaign
([generic-backend-qualification](../generic-backend-qualification/RESULT.md))
was operationally censored at the six-hour job cap with a failed artifact
upload. This directory freezes the required replacement: a smoke-first,
source-bound F4/F5/SAT complete-DLP comparison on a new seed that excludes every
reconstructed point from seed `2026092901` and every sealed prior exposure.

| Artifact | Role |
| --- | --- |
| [PROTOCOL.md](PROTOCOL.md) | Hypothesis, freeze, executed schedule, exclusions, job envelope |
| [panel.json](panel.json) | Exact arms, seed `2026092902`, 210 executed pairs, soft 300-minute measure wall; byte SHA-256 `319d8f624c7c5dcbc156ba97f85ea34850352c7d307aa06c9906f9a2d4010e91` |
| [censored-2026092901-fixtures.json](censored-2026092901-fixtures.json) | All 25 potentially exposed points from the censored prepare path; SHA-256 `30f496a561aa0814adbf2822d94ec4e98c808e1a3fa9219d592b182521438004` |
| [RECONSTRUCTION.md](RECONSTRUCTION.md) | How those points were regenerated without redispatching measurement |
| [reconstruction-receipt.json](reconstruction-receipt.json) | Compact digests and provenance |

Success is smoke family qualification (every smoke job for at least one F4/F5
arm and one SAT arm). Development fixtures are prepared and frozen as
exclusions but not executed here. Child timeout is 180 seconds. The workflow
must soft-stop measurement at 300 minutes, package one `tar.zst`, and upload
that archive inside the 360-minute job limit.

Do not redispatch seed `2026092901`. Do not treat admission controls, readiness
panels or the earlier generic/reference pair-table study as substitutes. This
registration consumes no improvement-round slot from the exhausted three-attempt
budget. Toy-panel qualification, even if positive, does not establish an
ECC2K-130 crossover or a globally fastest index-calculus method.
