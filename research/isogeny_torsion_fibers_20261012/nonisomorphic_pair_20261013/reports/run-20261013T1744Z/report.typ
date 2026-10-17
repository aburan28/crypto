#set page(paper: "us-letter", margin: 0.62in)
#set text(size: 9pt)
= Non-isomorphic pair V2: clean direct execution
2026-10-13 17:44 UTC · shared Darwin arm64 · L0 · 24/24 verified arm rows

The same eight public source targets were solved by source relations,
degree-three destination relations and matched source negation rho. Both
direct executions returned zero, and their raw TSV and stderr logs compare
byte for byte. The independent Rust checker verified all recovered scalars,
the non-isomorphic map, its rational three-point kernel and eight target
fibers. Exact commands, hashes and input rules are in `README.md`.

#image("group_ops-v2.svg", width: 100%)

Across eight targets, rho used 664, source relations 3,312, and destination
relations 3,794 counted group operations. Destination/source = 1.145531 in
this unit. The map field work, roots, table probes, row work and memory are
recorded separately and lack a pinned conversion into the graph's unit.
