# Fresh matched degree-83 IC versus optimized rho

Status: prepared and source-frozen; the main timed campaign has not started.
No new matched-panel winner is reported yet. User authorized execution
on14October2026. Worker/controller/analysis source freeze:
`c6e8daafb5a59d27e1be88a2f61f0bd6ee83ef03`.

The [prospective protocol](PROTOCOL.md) freezes the exact legacy a=1/tail149
curve,32 new point-only targets, three genuinely reused-setup IC/rho cohort
pairs, eight actual singleton cold pairs, and five rho A/A pairs. IC uses the
retained K600 and rank seed29; rho uses the previously selected Fold/reference
u64/256-lane/DP6/escape-off policy. Both are single-threaded under one native
CPU3 reservation,12GiB address-space limit and2100s full-process budget.

The primary costs include setup, solving, recovery and built-in verification.
Cold measurements start from the retained K600 file; historical factor-base
selection/construction and OS-cache flushing are outside this boundary.
Independent external arithmetic/file audit is separately metered. Producer
online boundaries differ and remain secondary. No old IC time is divided by
a newly measured rho time.

| User requirement | Current status | Evidence or remaining gate |
|---|---|---|
| Larger unseen matched panel | Fixtures independently verified | [Generation and native replay](verification/README.md) |
| Freeze both algorithms before testing | Verified before main timing | [Protocol](PROTOCOL.md), [build receipt](provenance/build_receipt.json), [frozen matrix](provenance/frozen_campaign_manifest.json) |
| Actual cold-start comparison | Not yet run | Eight predeclared singleton pairs |
| Genuine setup reuse for both methods | Implemented; batch campaign not yet run | New rho singleton matches existing worker's exported work; three32-target cohorts planned |
| Complete cost, memory and failures | Controller tests pass; campaign pending | Nine tests, including native timeout and launch/memory failure handling |
| Means, medians, slow cases and uncertainty | Native analysis14tests pass; measured input pending | Target/preparation clusters and non-finite modeled crossovers retained |
| Independent correctness and provenance | Fixtures pass; solve audit pending | Scalar and IC four-point witness replay |
| Comprehensive PDF and reviewable publication | Pending executed evidence | New addendum preserving prior70pages; extend PR1828 |

Previous reports and all frozen studies remain historical evidence. This
campaign is a local diagnostic, with calibrated S, L2/physical-host admission,
independent-host claims and global parameter optimality unset. A complete
matrix and independent verification are required for a full-workload winner.

Conductor is unavailable; the pre-edit scope attempt and exit127 are retained
under `checks/conductor.*`. Separately, current PR-CI repairs install the Lean
workflow's missing ripgrep dependency and make curve-name diff checking honor
the existing frozen-archive exemption without modifying archived bytes.
