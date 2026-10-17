#set page(paper: "us-letter", flipped: true, margin: 0.6in)
#set text(size: 11pt)
= Third paired run · 2026-10-13 01:59 UTC
The executable verified 48/48 solver arms and both explicit maps on the
shared Darwin arm64 host. Its output schema failed: the TSV header has 47
fields and the rows have 48. The raw bytes, map log and identical same-host
replay remain at `../../raw-v3.tsv`, `../../run-v3.log`, and
`../../replay-v3.tsv`. SHA-256 of the raw TSV begins `77b7fa9f`.

#image("flow.svg", width: 100%)

No comparison is drawn from v3. The v4 runner checks its own header and
row widths, then reruns every frozen target.
