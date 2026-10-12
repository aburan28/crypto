# Decision: O-aware rational S3 CNF semantics passes; dense table is not a scale-up arm

The frozen one-hot DIMACS exporter in [PR #781](https://github.com/aburan28/crypto/pull/781)
implements the O-aware local table proved in [PROOF.md](PROOF.md):
`O+a→a`; finite `u+a→` every finite `S3(u,a,c)` root; and
finite `u+a→O` iff `u=a`. Every factor slot is a complete rational x
fibre. A final-state unit literal selects an exact full-point target,
including O; both signs of an affine target are realizable by global
negation. An independent verifier rebuilt every local table from full
point additions, parsed the emitted DIMACS and compared its complete
clause list, then enumerated every signed tuple and compared every CNF
model path and all 61 target assumptions. All checks passed. This is an
exhaustive semantic result for the frozen domains, not a solver run.

The pre-outcome source commit was `b27ce0fa6241dbf49086ea9d2d2ea12671927a0e`.
The [frozen manifest](FROZEN.json) SHA-256 is
`684553d2a05f7355e29acc6b3f76dc3998a8cd5f970b6aafcd508c3fe2d02cf2`.
The [protocol](PROTOCOL.md) fixed the four panels, #767 factor archive,
#774 Q+T target list, negative controls and 180 s/512 MiB child caps
before this first and only outcome. The draft PR hash-only and workflow
syntax CI passed before the run; there was no failure or retry.

| Complete panel | Factor-x tuples | Signed point tuples | CNF/model paths | Paths with O prefix | Paths ending at O | Full-point target labels | Variables / clauses |
|:--|--:|--:|--:|--:|--:|--:|:--|
| n=2, m=4, all fibres | 256 | 2,401 | 949 | 184 | 133 | 8 incl. O | 31 / 258 |
| n=3, m=5, all fibres | 32 | 243 | 70 | 49 | 20 | 4 incl. O | 22 / 66 |
| n=4, m=4, all fibres | 4,096 | 50,625 | 22,373 | 1,744 | 1,485 | 16 incl. O | 59 / 1,744 |
| n=13, m=5, rotated #767 fibres | 243 | 3,125 | 1,563 | 63 | 0 | 32 fixed Q+T + O | 1,263 / 1,195,344 |
| **Total** | **4,627** | **56,394** | **24,955** | **2,040** | **1,638** | **61** | — |

O-prefix and O-terminal path columns may overlap and are not disjoint
classes. Every one of the 24,955 parsed-CNF paths is realized by a signed
point tuple, and every complete signed tuple projects to a listed CNF
path. For every finite path, the verifier found **both** rational full
terminal signs (the x=0 singleton is its own sign). The all-point toy
panels accepted their exact O targets. The n=13 synthetic O target had
zero models and zero full-point witnesses in its complete 3,125-tuple
oracle; this is an exhaustively certified toy negative, not a result from
SAT. Every one of the 32 prior Q+T target literals was checked against
the new independent oracle.

The fixed n=13 target index 12, `Q3T0=(7256,3272)`, has three CNF paths.
One passes through O immediately after the first two factors and has
factor-x tuple `[0,0,0,6433,217]`, corresponding to #774's
exceptional-only mask `[0,0,0,2,1]`. Thus the concrete O state recovers
the witness omitted by the affine-only chain. The other two paths for
this target are affine. All four frozen negative controls were rejected:
nonliftable n=3 x=2 (`DomainError`), a single-sign n=2 x=1 fibre
(`DomainError`), an inserted clause forbidding `(0,0)→O`
(`ClauseError`), and an O-target unit replaced with the finite x=0
literal (`TargetError`). No corrupt or sign-incomplete input was
mistaken for a certified PDP negative.

| Cold child | Wall s | CPU s | Peak RSS bytes | Frozen cap |
|:--|--:|--:|--:|:--|
| Producer | 1.023 | 0.911 | 237,912,064 | 180 s / 536,870,912 bytes |
| Independent verifier | 2.071 | 1.947 | 495,173,632 | 180 s / 536,870,912 bytes |

The n=13 dense CNF alone is 17,630,779 bytes. Its 1,195,344 clauses
split into 462,381 one-hot clauses (mostly the 928-way terminal-state
pairwise exclusion) and 732,963 forbidden-transition triples. The
verifier passed with only 41,697,280 bytes below its 512-MiB cap.
Those measured costs reject **this dense table layout** as a sensible
next scale-up candidate even on the toy rung; they do not bound every
sparse or implicit encoding and say nothing about n=131 solver time.
The next representation gate should encode allowed transitions sparsely,
with auxiliary selectors or another compact scheme, and prove equivalence
to this exact CNF/full-point oracle on the same fixed labels before any
SAT/FES/Gröbner timing. A separate direct S6/S7 exporter and exact model
lifting remain necessary for a direct-versus-chain solver comparison.

[The raw evidence](evidence/README.md) includes all four solver-ready
DIMACS bases, schemas with target assumption literals, compressed complete
model paths, full producer/verifier results, stdout/stderr, cold child
receipt, and byte/SHA-256 manifest. Its first receipt SHA-256 is
`6d74ec53d6615eb1e17d38d35ebb2fadbf4023f63263b331e5ea95d35887d766`.
No end-to-end index-calculus cost S, n=131 relation yield or rank,
challenge logarithm, or matched-rho ratio has been measured here.
