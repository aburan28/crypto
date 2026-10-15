#set page(paper: "a4", margin: 16mm)
#set text(font: "Arial", size: 9pt)
#set heading(numbering: none)
= Post-merge validation and named-curve controls
2026-10-15 · Isogeny performance study · correctness receipts

Source for post-merge tests: #text(size: 8pt)[0191f5b924765a0f6d7882f13787b7391413c685].
The arithmetic implementation is unchanged from the final BMSS rung.
Apple M4 Pro, 14 logical cores, 48 GiB, shared macOS host; actual compiler
Homebrew rustc 1.93.1. No qualified wall timing is reported here.

== Full root and standalone validation
#table(columns: (2fr, 1fr, 1fr, 2fr), inset: 5pt,
  [Attempt], [Passed], [Failed], [Result / boundary],
  [Two test threads], [4,094], [1], [117 existing ignored; Nagao deadline expired],
  [Exact Nagao rerun], [1], [0], [One Rayon worker; all assertions retained],
  [Concurrent Rayon-one retry], [Incomplete], [Unknown], [Interrupted; no full-suite verdict],
  [Full serial Rayon-one suite], [4,095], [0], [117 existing ignored; terminal PASS],
)

The first failure explicitly reports `timed_out: true`. The exact test and
solver source are unchanged. The successful serial run retains every test,
assertion and deadline. The independent shell-exit receipt for that completed
serial tool session was no longer retrievable on resumption; its complete
terminal test verdict is retained. Full standalone release library, binary,
integration and example suites also passed. Raw failed and interrupted receipts
remain beside the successful run. Thread interference is a supported hypothesis,
without a native isolation receipt proving its precise cause.

== Dedicated-binary consumer controls
#image("CONTROLS.svg", width: 100%)

Both controls use the embedded public catalogue and independent verifier.
Consumer build source: #text(size: 8pt)[6e90a71a72b7b7c3522bb0307ad1d25299bc98ce].
P192/13 enumerates both eigenlines; P224/1117 reports one eigenline with PARTIAL
coverage. Every map passes kernel, codomain, exact rational-map, subgroup and
twenty public scalar-transport checks. These controls are excluded from the
separate 2,000-additional-map requirement. They are correctness runs, not paired
performance samples; a timing graph would have no admitted comparison.

P192 certificate SHA256:
#text(size: 7pt)[7a1809d56cb774f681e06ecf1dad9e04d97abe70a3eb2738bff9d54affc861d1]
P224 certificate SHA256:
#text(size: 7pt)[6f1812ef5c90fabc34093c8b9545053997b55f851c44a9883bfa1ea2842d05fd]

Complete commands, raw logs, maps, stderr, resource receipts and independent
replay records are retained in this directory. Read README.md and FILES.json.
The earlier frozen archives and their inventories are unchanged. Multi-CPU
constructor comparisons and the remaining algorithm families remain open.
