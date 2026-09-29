# Second registered generic-backend campaign: operationally censored

The single measured [Actions run 36580669479](https://github.com/aburan28/crypto/actions/runs/36580669479),
attempt 1, executed workflow checkout
`e0fadb5ca5e99a3815c36a127fbdd1d1a5d12b43` with frozen generic worker
`765c3c5f19032bd852163805f257c56babef2040` for panel seed
`2026092902` (panel SHA-256
`d283a869b0412228d1c66260fdfd8f387d7243bd15456c7febf3c46ee5da27a8`).
It did **not** produce an auditable tournament bundle. The measured step
started at 2026-09-29 14:15:33 UTC and GitHub stopped it at 19:15:46 UTC
when the 300-minute step timeout fired. The reserved `if: always()` pack step
then tarred a ~35.4 GiB ephemeral tree while measure processes were still
alive; `tar` reported files removed and directories changing during read
(including live `callgrind.out.*` under development trees), `zstd` aborted,
`pack_partial_campaign.py` deleted the partial archive and raised
`RuntimeError('partial campaign archive failed')`, and orphan runner processes
were terminated only after packing finished. `actions/upload-artifact@v4`
failed with `if-no-files-found: error`. GitHub listed only the PR smoke
artifact `ic-generic-pack-smoke-36580669479` (917 bytes), not the campaign
archive `ic-generic-backend-qualification-v2-36580669479-1`. The retained
[complete campaign job log](workflow-job-109448846736.log) is 3,378,181 bytes,
SHA-256
`d8a7931f7b43922d2dc7968a4df2d5aa36f4804c5ea140738b9a97a28270c8d3`.

Without a checksummed archive, no frozen `tournament.py verify`,
`natural-yield.json`, `generic_backend_gate_v2.py`, or
`independent_pairs.py` cross-check can be applied. The workflow log’s
inline progress JSON is **not** a substitute for retained receipts; it
only bounds what the runner reached before interruption.

| Registered question | Required evidence | Retained result |
| --- | --- | --- |
| Full 250-pair schedule | Frozen receipts and verifier on complete bundle | Unknown; bundle unavailable |
| F4/F5 family qualification | At least one verified arm on every smoke and development job | Not established |
| SAT family qualification | Same | Not established |
| Natural ordinary-query yield | Independently replayed bounded reports | Unknown |
| Same-point IC/incumbent/rho comparison | Complete source-bound costs and verified online intervals | Unknown; no speedup ratio |

**Log-only bounds (not verified receipts).** Before step timeout the
runner emitted 82 progress lines: 46 `VERIFIED` and 36 `TIMEOUT`, covering
all 10 A/A slots, all 60 smoke slots, and only the first 12 of 180
development slots (development `n17a1-000` incomplete). No logged
`VERIFIED` row names a generic F4, F5, inherited F4, or SAT arm. Pair-table
generic arms and prepared references did verify on several cells. These
counts do not prove family failure or success; they only show the schedule
was incomplete and generic backends hit the 300-second child cap repeatedly
on this panel.

This is an **operationally censored campaign**, not a negative result about
F4, F5, SAT, or the optimized incumbent as mathematical algorithms. No
candidate or rho row gains a measured `S`, online time, confidence interval,
or correctness claim in the scoreboard from this dispatch. The
post-registration [static feasibility audit](STATIC-FEASIBILITY.md) shows the
registered ambient F4/F5 layouts are encoder-unsupported at `m=3`; that
static finding is separate from this censoring event and does not substitute
for measured family qualification.

Do not dispatch seed `2026092902` again. Treat all 25 points in
[lost-campaign-exposures.json](lost-campaign-exposures.json) as exposed together
with every public point this registration could have reached once reconstructed
(prepare replay at the pinned checkout), even though no bundle was uploaded.
A new campaign needs a new frozen panel and seed, exclusion of all three sealed
improvement archives, the 25 lost first-run points, and every point from this
second dispatch, plus packaging that quiesces or excludes live callgrind trees
before `tar` runs. Neither the sealed improvement confirmation sets nor this
interrupted run can be used as a favourable retry.
