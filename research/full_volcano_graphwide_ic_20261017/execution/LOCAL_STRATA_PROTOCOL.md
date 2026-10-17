# Conservative local strata, before applying to the C measurements

Keep the original prefix and its global NOT_ESTABLISHED stratum receipt
unchanged. Use the ordinary-volcano floor-distance correspondence already
implemented in `tools/volcano-cli/src/lib.rs::stratify` (D_K < -4, conductor
prime ell != p, positive depth v_ell(f_pi)). A completely expanded vertex
with one rational degree-ell kernel is a floor. Retain kernel multiplicity
when counting degree, including parallel maps and loops.

On the finite visible graph, the distance to a known floor is an upper bound
on distance to a floor. Distance to the union of known floors and unexpanded
frontier vertices is a lower bound: any omitted route must first leave through
the frontier. An unexpanded vertex is not certified a floor. Certify a height
only when both finite distances agree and are within the known volcano depth.
Convert all certified prime heights to the conductor with the existing
floor-distance formula. Otherwise retain missing conductor and explicit bounds.
Never infer complete class or stratum populations from these local labels.

Before using this annotation on C, calibrate against the three complete graphs
at visible id cutoffs 0, 1, 17, 1024 and all. Check every lower/upper bound and
every exact label against the existing complete-graph conductor. Require zero
disagreements. Read immutable databases and bind output to the DB SHA-256 and
existing map-replay receipt. This is a second labeling implementation using
the same standard volcano theorem and map arithmetic, not an independent
mathematical or arithmetic implementation proof.

For C use its fixed 1,869-vertex prefix and exact 1,024 expanded flags. Keep
local exact conductors, unknown conductors and all per-prime bounds separately.
Any later measurement join must match the DB hash and exact vertex ids. Report
sampled strata only; retain the original expected 165 strata and
1,934,901,248 vertices as unresolved full-class denominators. Annotation is
study preparation, not work used by the DLP solver and not included in solver
elapsed intervals. This is an implementation of existing floor-distance
reasoning, not a new algebraic identity or speed claim.


Generic-command extension, 2026-10-17: the same calculation is exposed as
`local-strata certify PREFIX_DB FRESH_CERTIFICATE`. It applies the unchanged
input/complete-kernel/depth checks to each retained original A–E prefix.
Before use, all fifteen existing complete-truth cuts pass and generic C
matches its earlier certificate byte for byte. This is a CLI extension and
certification run, not a new DLP benchmark or change to the frozen C evidence.
