#set page(paper: "us-letter", margin: 0.7in)
#set text(font: "Arial", size: 10pt)
#text(size: 18pt, weight: "bold")[Run record: S3 and seven-arm held-out panel]

*Date:* 2026-10-16. *Host:* Darwin arm64 shared L0. *Revision:*
`codex/isogeny-portable-trend-20261015`. The frozen V2/V3 specs
record exact solver, checker, meter, script, public and private input
SHA-256 values.

*Hypothesis:* source and quotient S3 pair support matches direct
enumeration, and their complete per-target implementation cost may
differ. Any exhaustive mismatch stops the panel.

*Exhaustive commands:* run `different_j_chi_resultant exhaustive-s3`
on `FOUR_TARGET_PUBLIC_V1.json`, then independently run
`different_j_chi_resultant_check check-s3` on its 229 case rows.
Both returned exit status zero. Each curve had 136 pairs and 130
zero-pair cases.

*Panel commands:* `bash RUN_HELDOUT_SEVEN_ARM_V3.sh` created 84
separate child records and meter records. The independent
`different_j_chi_heldout_check` replayed both curves, 664 relation
witnesses, row ranks, BSGS and rho. The native summary read all
child records. All 84 children passed.

*Limits:* 300-second hard wall timeout per child; peak RSS checked
against 512 MiB after completion. Darwin rejected address-space
and data-size hard memory limits in the earlier V2 checks. This
host has no isolation or thermal receipt. Complete portable S,
L2 admitted wall time, direct inverse-Vélu and full-image relation
controls remain open.

See `HELDOUT_RESULTS_V3.md`, `HELDOUT_SEVEN_ARM_V3.summary.json`,
`HELDOUT_SEVEN_ARM_V3.check.json`, `S3_EXHAUSTIVE_V2.check.json`,
the raw directory and the hash manifest.
