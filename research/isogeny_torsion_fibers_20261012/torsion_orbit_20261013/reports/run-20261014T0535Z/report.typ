#set page(width: 8.5in, height: 11in, margin: 0.7in)
#set text(font: "New Computer Modern", size: 10pt)
= Run: three-torsion orbit support

2026-10-14 05:35 UTC. Protocol commit `1d01a39a25`. Native Rust exact
enumeration of three arms on each of four prime-order subgroups.

The 12-point and 36-lift supports were identical: 64/103, 74/439,
75/1613 and 78/6427. The 36-distinct supports were 103/103, 351/439,
559/1613 and 642/6427. The complete source, raw data, command, hashes,
and checks are in `../../RUN.md` and `../../raw-v1.tsv`.

#image("../../graph-support.svg", width: 100%)

This run counted distinct reachable subgroup elements. It did not time a
complete solver or produce an end-to-end speed ratio. An initial compiler
type-inference failure was retained; the corrected executable completed
all twelve rows and its assertions passed.
