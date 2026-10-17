#set page(paper: "a4", margin: 1.8cm)
#set text(size: 9.5pt)
= A20: the folded build's streams written a line at a time (smoke only)

*Candidate and smoke, no exploration.* `LineStream` (`line-streams-on-a19.patch`, one commit on A19's candidate):
- each partition stream keeps a 64-byte line in L1;
- full lines go out with non-temporal stores;
- the order, and so the table, is unchanged.

*Tests:* 207 passed (`koblitz_index_calculus`, `koblitz_fast`).

*Smoke:* two isolated processes per arm at 2^44.3 and at 2^47.2 (`smoke.tsv`). The builds overlap; there is no visible gain, so it is not pursued.

*Probable cause, not measured:* the append is bound by first-touching 45–75 MB of fresh stream memory a build, not by line traffic.

#figure(image("build.svg", width: 90%), caption: [Each process's build time (median repetition), A19 against A20. Values from `smoke.tsv`.])
