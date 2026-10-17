#set page(paper: "us-letter", margin: 0.58in)
#set text(size: 9pt)

= Non-isomorphic isogeny transport with rational three-point fibers
2026-10-13 · T-2405 · F617 · target subgroup 103 · 24/24 verified arm records

The explicit degree-three map sends `E: y²=x³+1` to
`E′: y²=x³+590`. Each curve has 618 rational points. The rational kernel
is `{O,(0,1),(0,616)}` and every rational image point has three
preimages. No field element has sixth power 590, so these short
Weierstrass models are not base-field isomorphic. The source generator
`(66,264)` and image generator `(483,13)` each have order 103.

Eight frozen public source targets were each solved through cold source
relations, degree-three destination relations, and matched source
negation rho. The relation factor bases consist of twelve distinct
points projected into the subgroup by `[6]`. Every row recovered a
scalar checked on the original source; destination rows also checked
their transported points. The independent Rust checker verified all
24 answers, reconstructed 16 relation solves, checked eight target fibers,
and tested 381,924 group-homomorphism pairs.

The chart uses `S=(group additions + doublings)/sqrt(103)` and shows
observed ranges over eight targets. This unit includes cold factor-base
projection, collection, recovery, descent, verification, and rho group
work. The two-class walk-step floor is `sqrt(pi/4)=0.886227`.

#image("group_ops-v2.svg", width: 100%)

#table(
  columns: (2fr, 0.7fr, 0.9fr, 1.1fr, 0.9fr),
  inset: 3pt,
  stroke: 0.3pt,
  [arm], [verified], [total ops], [mean S], [vs rho],
  [negation rho], [8/8], [664], [8.178], [1.000],
  [source relations], [8/8], [3,312], [40.793], [4.988],
  [destination relations], [8/8], [3,794], [46.729], [5.714],
)

#pagebreak()
= Explicit transport, phase accounting, and limits

#image("transport-v2.svg", width: 100%)

Source and destination collection used 231 and 290 trials across the
same targets. Their counted-total ratio is `3,794/3,312=1.145531`;
destination relations used more group operations on five targets and
fewer on three. The 482-operation difference is setup `-24`, collection
`+408`, descent `-1`, and verification `+99`; the destination
verification includes a source cross-check. The per-target paired table
and phase table are `paired_summary-v2.tsv` and `phase_summary-v2.tsv`.

V1's native child completed 24 rows, but its `/usr/bin/time -l` wrapper
exited 1 after a `sysctl kern.clockrate` denial. The frozen remediation
reran the same binary directly. V2 and its repeat each exited zero;
their raw TSV and stderr logs match byte for byte. The independent
checker receipt is `independent_v2.log`. The source, binary and raw hashes
and exact commands are in the run reports and `validation/README.md`.
The eight explicit preimage triples are in `fibers-v2.tsv`.

Lean compiled `Nonisomorphic617.ratioNotSixthPower` and
`targetCoprimeKernel` without axioms, and `kernelScalarVanishes` using
standard `propext`. These finite results support the non-isomorphism
and subgroup arguments; the concrete elliptic
group law, isomorphism classification, kernel, and fiber equalities
remain checked in Rust and written mathematics, not formalized in Lean.

The host was shared Darwin arm64, L0. Map field work (`10M+4S+2I` per
destination target), roots, table lookups, modular row work, hashes and
memory have no pinned common conversion. Map discovery and exhaustive
preflight are also outside these target rows. This is a one-size counted
group-operation diagnostic; full native cost, growing sizes and an
independent-host replay remain open. See `REPORT.md` for exact identities,
source-linked evidence, proof map, and requirement status.
