#set page(paper: "us-letter", margin: 0.7in)
#set text(font: "Arial", size: 10pt)
#text(size: 18pt, weight: "bold")[Run record: nine-arm held-out fibers]

*Date:* 2026-10-16. *Host:* Darwin arm64 shared L0. *Revision:*
`codex/isogeny-portable-trend-20261015`. The frozen V4 spec
records exact source, binary, meter, orchestration, public and
checker-only input SHA-256 values.

*Hypothesis:* inverse-Vélu and full-image constructions select
the published destination base and have the same relation support
as the χ, direct-pair and source/quotient S3 arms. Their cold cost
may differ. See `HELDOUT_NINE_ARM_PROTOCOL_V4.md`.

*Commands:* `bash RUN_HELDOUT_NINE_ARM_V4.sh` created 108 separate
native child and meter records. `different_j_chi_heldout_check`
replayed both curves, 996 relation witnesses, row ranks, BSGS and
rho. The native summary asserted identical six-arm trial and
relation support. All checks passed.

*Preserved correction:* an initial compilation hit `json!`'s
recursion limit while building the new inverse child output.
Splitting the nested value fixed compilation before the V4 run.

*Limits:* 300-second hard wall cap per child and 512 MiB
observed Darwin RSS ceiling. No child failed or timed out.
The host has no isolation or thermal receipt. Complete portable
S and admitted wall time remain open.

See `HELDOUT_RESULTS_V4.md`, `HELDOUT_NINE_ARM_V4.summary.json`,
`HELDOUT_NINE_ARM_V4.check.json`, the raw directory and hash manifest.
