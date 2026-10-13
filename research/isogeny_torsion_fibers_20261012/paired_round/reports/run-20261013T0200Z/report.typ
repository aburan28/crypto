#set page(paper: "us-letter", flipped: true, margin: 0.5in)
#set text(size: 10pt)
= Final paired run · 2026-10-13 02:00 UTC
16 public source targets, three arms each: 48/48 verified answers. The
degree-three maps passed 361 and 44,521 ordered-pair checks. Frozen protocol
`../../PROTOCOL.md`; raw/replay SHA-256 begins `995b7bcc`. Shared Darwin
arm64 host, L0. All 47 TSV fields match each row. The replay and map logs
compare byte for byte on the same host.
An independent-code native point-law verifier passed all 48 rows and
rechecked the four curve orders and the cubic automorphism on 89,764
ordered point pairs. The later structural audit checked scaling
isomorphisms on 44,882 ordered pairs; receipt in `../../validation/`.

#image("../../group_ops-v4.svg", width: 100%)

Mean counted-group S at order 211: rho 5.938; source relations 50.186;
transported relations 41.332. Eight targets per arm; whiskers are observed
minima and maxima. Map field work, roots, lookups and row arithmetic remain
unpriced in S. The phase table is `../../phase_summary-v4.tsv`; the full
interpretation and next admission test are in `../../REPORT.md`.
