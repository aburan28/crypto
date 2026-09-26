# Decision: toy full-point DAG relation passes; n131 remains not admitted

The frozen K0 point-addition predicate passed exhaustive full-point controls
at n=2 and n=3. For every canonical rational P,Q,R and every slope λ, the
DAG accepted exactly the reference group-law output, with the specified
number of witnesses. All off-curve and noncanonical point-role controls and
malformed model controls rejected. A fresh archive-only verifier replay of the
final raw rows passed using the workflow's relative evidence path. This is a
representation-correctness result, not a rotated PDP or ECDLP speed result.

| Field | #E | PQR rows | Valid λ evaluations | Invalid-role evaluations | Variables / limbs | XOR / AND nodes | Errors |
|:--|--:|--:|--:|--:|:--|:--|--:|
| GF(2²), `0x7` | 8 | 512 | 2,048 | 288 | 17 / 1 | 184 / 107 | 0 |
| GF(2³), `0xb` | 4 | 64 | 512 | 2,976 | 24 / 1 | 280 / 173 | 0 |

Each field separately exercised identity-copy, inverse, double, and generic
addition cases: the n2 pair counts were `(copy Q,copy P,inverse,double,generic)
= (8,7,7,6,36)`; n3 counts were `(4,3,3,2,4)`. The producer and independent
verifier agreed on all 576 compressed raw rows and 5,824 model evaluations
(2,560 valid-table plus 3,264 invalid-role). Ten model controls per field
passed. The synthetic 920-input expression DAG read limb 15 correctly and
rejected a missing limb and nonzero padding; it does not instantiate an n131
point circuit. The model width formula for one three-point addition edge plus
one slope is `7n+3`, hence 920 Boolean inputs and 15 64-bit model limbs at
n131. That count is a representation requirement, not a solver-cost bound.

| Final cold child | Child wall | CPU | Peak RSS | Outer wall |
|:--|--:|--:|--:|--:|
| Producer | 0.223 s | 0.210 s | 26,443,776 B | 0.365 s |
| Independent verifier | 0.306 s | 0.230 s | 24,952,832 B | 0.396 s |

The 512-MiB RSS value in the protocol was an acceptance gate checked after
each child exited, not a hard process-memory limit. The 180-second child wall
and 195-second external wall gates likewise passed. These local Python 3.13
costs are toy validation costs; no SAT run, relation collection, linear algebra,
or full-DLP comparison was performed. The attack boundary remains
matched automorphism-aware Pollard rho with all phases charged; end-to-end S
and speedup are unset.

The original frozen head `c7c8a10` produced the same raw truth-table SHA-256
`a2b760b72e477cb9ab4d8409e71102cd6a2a09db0b8a5e8f5ad170d1a1102da9`
and independently verified its relation, but its archive-only replay failed
because a relative evidence path was resolved after changing the child working
directory. That complete first attempt, original freeze, failure log, and
receipt are retained in `evidence_failure_0/`. Head `3a9eaf7` changed only
the archive path resolution and added the same-x completeness proof, then
froze again before this final cold run. The final freeze SHA-256 is
`16065b7dd0bb45ac942bc2622adf8451464b6457c35282a40f7b4c07ba0655a9`.
No first-attempt data was overwritten or relabelled.

Next admission gate: connect this DAG to a complete rotated factor-domain
constraint and an independently checked CNF/ANF exporter plus model lift on a
frozen n13/n19 panel. Compare its full export and solver cost with the
explicit O-aware onehot baseline. n131 stays `NOT_ADMITTED` until that bridge,
width, memory, useful-support, and end-to-end cost gates have evidence.
