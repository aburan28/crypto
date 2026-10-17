#set page(paper: "us-letter", flipped: true, margin: 0.5in)
#set text(size: 10pt)
= Second paired run · 2026-10-13 01:56 UTC
48/48 arm answers and both exact maps verified. Shared Darwin arm64 host, L0.
The frozen inputs and command are in `README.md`; raw record SHA-256 begins
`bdfaafb5`. Rho's phase is separated, but a reported zero for its table
lookups was only a placeholder: the library does not expose that count.
v4 removes the misleading column.

#image("../../group_ops-v2.svg", width: 100%)

The chart shows counted group-operation S, mean and observed range over eight
public targets per field. Native map, lookup, root and row work remains
separate. See `../../raw-v2.tsv` and the final `REPORT.md`.
