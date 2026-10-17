#set page(paper: "a4", margin: 12mm)
#set text(size: 8.5pt)
#set heading(numbering: none)
= Three complete ordinary isogeny classes
Additional 57/59-bit examples; run 20261017. Original A/C/D/E remain incomplete.

#table(columns: (3fr, 1fr, 1fr, .6fr), inset: 4pt,
  [Source model], [Vertices], [Directed maps], [Strata],
[icv1-fp57-t536870914-6cabb0c9], [8,192], [16,384], [14],
[icv1-fp57-tm536870941-e6d0c766], [6,561], [13,122], [9],
[icv1-fp59-t1073741954-4dd0b1ef], [10,070], [39,794], [28],)

#image("strata.svg", width: 100%)

*Verified coverage:* 24,823 vertices, 69,300 directed maps, 51 conductor strata.
All selected-degree outgoing kernels and exact map identities replayed;
all conductor populations agree. Original database hashes survived replay
on copies. Independent implementation replay: NOT_RUN.

Maximum endomorphism conductors are 8192, 6561 and 1728. Floor-model inputs
are included. Sources start at conductor 1. SQLite and CSV/JSONL exports
contain the entire vertex class and all degree-2/3 maps in the stated sets.
Every stratum has a representative with explicit outgoing rational maps.
Arbitrary-degree maps are an infinite family, outside this finite export.

*Evidence:* evidence/coverage.json, *-build.json, *-verify.json,
*-strata.json and *-stratum-edges.json. Exact integer census; uncertainty
intervals do not apply. Protocol: PROTOCOL.md, commit 3b7f4b39c; frozen
inputs/run commands: 946ab7419. Linux x86_64, 4 CPU quota, 16 GiB RAM.
Construction/replay share native arithmetic. Conductor executable missing (127).

*Still open:* full A/C/D/E generation, independent implementation audit,
graph-wide ECDLP min/max/mean and operation counts, direct EC index calculus,
factor-base yield and finite-field reductions. No solver timings or speedup
claims are made; IC scoreboards retain their prior performance values.
Full report, exact inputs, commands, hashes and provenance accompany this PDF.
