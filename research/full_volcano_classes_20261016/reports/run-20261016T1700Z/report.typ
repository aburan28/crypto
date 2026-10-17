#set page(paper: "a4", margin: 15mm)
#set text(size: 9pt)
#set heading(numbering: none)
= Complete isogeny classes: implementation and coverage
Native enumeration controls and fixed-input preflights, 16 October 2026.

The CLI is implemented and the bounded runs below passed. The full five-case
dataset request remains open. A complete vertex class and a complete graph
for selected prime degrees are separate claims. Arbitrary isogeny degrees
form an infinite family.

#image("coverage.svg", width: 100%)

== Observed correctness evidence
#table(columns: (2.5fr, 1fr, 1fr, 2fr), inset: 4pt,
  [Source model (ICV1)], [Vertices], [Directed maps], [Coverage],
  [icv1-fp7-t2-4910f2d0], [8], [16], [Complete degree-2 class; 4 strata],
  [icv1-fp8-tm1-1ed3feb7], [9], [18], [Complete degree-3 class; 3 strata],
  [icv1-fp7-t4-0c09e92f], [10], [36], [Complete degree-2/3 class; 4 strata])

All three vertex sets agree with an independent exhaustive census of every
nonsingular model in the corresponding small field. Six native release tests
also cover a nontrivial surface class group, checkpoint/resume equivalence,
corrupt-map rejection, missing-degree components, BSGS and subgroup/budget checks.
Two binary controls also check a killed writer's recovery/resume against an
uninterrupted graph and preservation of existing output files.
Raw test output: evidence/native-tests.txt. No method timing or speed comparison
is admitted by these correctness tests.

Each exact large input was expanded 1,024 times and replayed. Those outputs are
partial; unexpanded frontier vertices remain. Replay regenerated every outgoing
kernel on expanded vertices and checked every exact rational-map identity and
target scaling. It shares native arithmetic with construction. An independent
implementation audit of the new generic backend has not been performed.

== What the command does
`crypto isogeny class build --curve curve.json --out graph.sqlite --steps all`
uses certified root order, exact class populations, complete rational-kernel
enumeration and atomic SQLite vertex expansions. The equivalent standalone
commands are `isogeny class ...` and `isogeny-class ...`.
Automatic completion triggers replay and conductor-population checks. Resume
keeps committed expansions. Vertex inspection returns explicit rational maps.

Backend limits: ordinary prime fields 3 < p < 2^61, -10,000,000 <= D_K < -4,
prime degrees <=31 different from p. Other cases fail explicitly. The source
endomorphism conductor is not inferred from the Frobenius conductor alone.

== Full datasets, storage and DLP obligations
The historical degree-3 case has all 129,140,163 vertices, 18 strata and
258,280,326 directed kernel maps in its separately distributed lossless archive.
Its compressed bytes and hashes are preserved in evidence/historical-B/.
Historical timings were produced outside this repository's measurement protocol.

The five full classes contain 7,833,658,967 vertices and 515 strata. A 24-byte
layout for a,b,j alone occupies 188,007,815,208 bytes, excluding edges, indexes,
journals and benchmark evidence. This is a fixed-layout baseline, not a measured
SQLite size estimate or an information-theoretic lower bound. Compressed or
generative encodings may be smaller. The workspace has 32 GiB total.
Full remaining class generation needs an appropriate durable execution volume
or another validated storage backend; it is not declared impossible.

The solve command verifies one explicit cyclic-subgroup BSGS logarithm. It does
not complete the requested graph-wide rho/BSGS/PH and direct/finite-field IC
campaign. All-vertex end-to-end timing, total operations, relation yields and
min/max/mean remain pending. No sampled result is labeled exhaustive.

Evidence and reproduction: PROTOCOL.md, run-controls.sh, run-prefixes.sh,
tools/volcano-cli/README.md, native plans, build/replay JSON and the named storage
problem. The existing class-number and ordinary-volcano theorems are reused;
no new mathematical theorem or Lean proof is claimed.
