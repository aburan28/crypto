# Second registered generic-backend campaign: operationally censored

The single measured [Actions run 36580669479](https://github.com/aburan28/crypto/actions/runs/36580669479),
attempt 1, executed workflow checkout
`e0fadb5ca5e99a3815c36a127fbdd1d1a5d12b43` with frozen generic worker
`765c3c5f19032bd852163805f257c56babef2040` for panel seed
`2026092902` (panel SHA-256
`d283a869b0412228d1c66260fdfd8f387d7243bd15456c7febf3c46ee5da27a8`).
It did **not** produce an auditable tournament bundle. The measured step started
at 2026-09-29 14:15:33 UTC and GitHub stopped it at 19:15:46 UTC on the
300-minute step timeout. The reserved `if: always()` pack step then tarred a
~35.4 GiB ephemeral tree while measure processes were still alive; `tar`
reported files removed and directories changing during read, `zstd` aborted,
`pack_partial_campaign.py` deleted the partial archive, and
`actions/upload-artifact@v4` failed with `if-no-files-found: error`. GitHub
listed only the PR smoke artifact
`ic-generic-pack-smoke-36580669479` (917 bytes), not the campaign archive
`ic-generic-backend-qualification-v2-36580669479-1`. The curated job log is
[workflow-job-109448846736.log](workflow-job-109448846736.log), SHA-256
`95598e7092a658c9cfd3da35d48e495041a129eadec92470da7276615695e265`.

Runner stdout shows prepare/verify steps and 46 per-job `VERIFIED` lines
through development `n17a1-000` before the timeout. Those lines are not retained
receipts: with no uploaded bundle there is no frozen
`tournament/evaluator/tournament.py verify`, no `natural-yield.json`, no
`family-gate.json`, and no durable target or failure export. The exact completed
slot count, last arm, relation matrices, and competitive costs remain
**unknown** for audit purposes even though the log suggests smoke finished and
development had barely started.

| Registered question | Required evidence | Retained result |
| --- | --- | --- |
| Full 250-slot schedule | Packed archive, frozen verifier, gate JSON | Unknown; bundle unavailable |
| F4/F5 and SAT family admission | Complete verified arm on every scheduled point | Not established |
| Natural ordinary-query yield | Independently replayed bounded reports | Unknown |
| Same-point IC/incumbent/rho comparison | Complete source-bound costs and online table | Unknown; no speedup ratio |

This is an **operationally censored campaign**, not a negative result about F4,
F5, SAT, or the optimized incumbent as mathematical algorithms. No candidate
or rho row gains a measured `S`, online time, confidence interval, or
correctness claim under repository accounting. The post-registration
[static feasibility audit](STATIC-FEASIBILITY.md) already shows the registered
ambient F4/F5 layouts are encoder-unsupported at `m=3`; that static finding is
separate from this censoring event and does not substitute for measured family
qualification.

Do not dispatch seed `2026092902` again. Treat all 25 points in
[lost-campaign-exposures.json](lost-campaign-exposures.json) as exposed together
with any targets this run may have generated while receipts were lost. A later
registration needs a new seed, a widened exposure census, a schedule that fits
the measured and pack envelopes, and packaging that quiesces or excludes live
callgrind trees before tar.
