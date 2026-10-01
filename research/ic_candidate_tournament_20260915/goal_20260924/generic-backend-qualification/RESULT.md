# First registered generic-backend campaign: operationally censored

The single measured [Actions run 36532455386](https://github.com/aburan28/crypto/actions/runs/36532455386),
attempt 1, executed the frozen PR #920 checkout
`765c3c5f19032bd852163805f257c56babef2040` for panel seed
`2026092901` (panel SHA-256
`83c640a03b4239b918851f6f1b8450e2fe27dc710fb99f306e3481a99d8875cf`).
It did **not** produce an auditable tournament result. The measured command
started at 2026-09-29 06:45:03 UTC and GitHub canceled it at 12:44:02 UTC as
the 360-minute job limit was reached. The `if: always()`
`actions/upload-artifact@v4` step then failed with `Maximum call stack size
exceeded` at 12:44:04 UTC. GitHub finalized the workflow as `cancelled` at
12:48:51 UTC, and the Actions artifact API listed zero artifacts after that
terminal state. The retained [complete job log](workflow-job-109289404217.log)
is 66,978 bytes, SHA-256
`05ffde03457696387e8acf4e46d7bdc76cf97bb8ad2d90832de228f324b5f325`.

The runner redirected each internal command's output into its ephemeral
workspace. With no uploaded bundle, the exact number of completed trial slots,
the last stage or arm, raw failures, targets, relation attempts, factor-base
inventory, and costs are **unknown**. The workflow log proves the cancellation
and upload failure; it does not prove any particular solver failed or passed.
No archived `tournament.py verify`, natural-yield audit, two-family gate or
paired incumbent/rho analysis can be applied to absent receipts.

| Registered question | Required evidence | Retained result |
| --- | --- | --- |
| Full one-target schedule | 750 native/profile receipts and frozen verifier | Unknown; bundle unavailable |
| F4/F5 and SAT family admission | At least one complete verified arm in each family | Not established |
| Natural ordinary-query yield and failed attempts | Independently replayed bounded reports | Unknown; no zero-yield inference |
| Same-point IC/incumbent/rho comparison | Complete source-bound costs and verified online intervals | Unknown; no speedup ratio |

This is an **operationally censored campaign**, not a negative result about F4,
F5 or SAT as mathematical algorithms. No candidate or rho row gains a measured
`S`, online time, confidence interval, or correctness claim. The already-merged
float-hash repair in [PR #941](https://github.com/aburan28/crypto/pull/941)
remains a future-auditor correction; this run never reached that audit.

Do not dispatch this registration again. The runner may have generated every
fixture before measurement and may have exposed more targets while the output
was lost. A new campaign must reconstruct and exclude **all potentially
generated public points** for this seed before any fresh target is sampled,
use a new frozen panel and source identity, and preserve incomplete evidence
before the job cap. Its workflow should time-limit the measured step with
enough time left to package a single archive file for upload; the new protocol
must make the expensive arm schedule fit that envelope. Neither the old
three-round confirmation sets nor the unretained 2026-09-29 run can be reused
as a favourable retry.
