#set page(paper: "us-letter", margin: 0.62in)
#set text(size: 9pt)
= Bounded isogeny control run v2
2026-10-13 16:32 UTC · shared Darwin arm64 · L0 · exit 2

The precommitted v2 setting distinguished every nonidentity class
representative after a rho step. The runner emitted all 80 intended rows:
77 verified and three six-class rho failures at F13 target indices 4, 5,
and 6. All relation and negation-rho rows verified; all F199 rows verified.
The chart shows the exact status of each source target.

#image("verification-v2.svg", width: 100%)

Raw: `raw-v2.tsv`, SHA-256 `626fb298…51c1b7415caa06e2e`.
A reconstructed v2 source reproduced both raw and log bytes under a
native rebuild; exact identity to the original v2 source and binary
bytes is unverified. Commands, parameters, and counters are in
`README.md`, `REMEDIATION_V2.md`, `validation/historical_replay.md`,
and the raw TSV. The all-arm acceptance gate was not met.
