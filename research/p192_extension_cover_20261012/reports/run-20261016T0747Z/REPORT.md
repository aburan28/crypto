# Three-step Frobenius trace for the P-192 Type I quotient

**Result.** Every rational Type I plane point in the inverse's open locus
maps through this explicit quotient formula to a point of
\(E(\mathbf F_{p^3})\). Its six-step \(p\)-Frobenius trace is therefore
twice its three-step trace. The existing six-trace cofactor scalar
\(k_6=h6^{-1}\bmod r\) can be replaced by
\(k_3=2k_6=h3^{-1}\bmod r\) when the anchor is also in
\(E(\mathbf F_{p^3})\). The frozen anchor passes that condition.
This removes three Frobenius conjugates and three extension-field point
additions from each trace evaluation.

The proof follows the actual inverse formulas. For rational plane
coordinates \(u,v\), its \(A_i=u+\zeta_i v+\zeta_i^2\) lie in
\(\mathbf F_{p^3}\). The inverse abscissa \(X\) lies in \(\mathbf F_p\),
and the quartic factors and norm roots lie in \(\mathbf F_{p^3}\), so
\(\lambda=A_0f_0(X)/(A_1A_2N(X))\) and
\(y_i=1/(\lambda A_i)\) lie there wherever the denominators are nonzero.
The quotient coordinates
\(c_0+d_0/(X-\rho)\) and
\(\kappa y_0/(X-\rho)^2\) also lie there because the model constants
are fixed by \(p^3\)-Frobenius. For \(P\) over that subfield,
\(F^3(P)=P\) and
\(\operatorname{Tr}_6(P)=2(P+F(P)+F^2(P))\).
The [checked Sage receipt](../../type_i_three_trace_subfield_256.json)
verifies every named model constant, the anchor, and both residual images
of 128 frozen quartic cuts. All 256 subgroup columns agree with the
six-step projector and the frozen earlier output.

| Diagnostic stage, 256 frozen images | Six-step | Three-step |
| --- | ---: | ---: |
| Trace calls, unisolated elapsed seconds | 11.213 | 7.604 |
| Frobenius conjugates and additions per trace | 6 | 3 |
| Frozen subgroup-column mismatches | 0 | 0 |

The local elapsed ratio is 1.47 for the trace stage only, with subfield
membership tests made outside both trace timers and three existing
producer processes running. It is a diagnostic for choosing the next
implementation, not a measured full-collector speedup. The running
141-pivot producer remains on its frozen six-step source; its current
receipts are unaffected.

The [Lean module](../../lean/TraceThree.lean) contains the following
proof map:

| Mathematical step | Lean declaration | Scope |
| --- | --- | --- |
| \(\operatorname{Tr}_6=2\operatorname{Tr}_3\) | `sixTrace_eq_two_threeTrace` | Additive Frobenius map; assumes \(F^3(P)=P\) |
| Cofactor scalar replacement | `threeTrace_cofactor_projection` | Assumes the existing six-trace projection identity |
| Projection from structural hypotheses | `threeTrace_cofactor_projection_from_conditions` | Uses the six-step theorem's Frobenius, order, scalar-congruence, and trace-annihilation hypotheses |
| Quotient-coordinate subfield closure | `recovered_quotient_image_cube_fixed` | Field endomorphism; fixed constituents and inverse open locus |

The [Lean verification receipt](../../lean/trace-three-verification-20261016.json)
records `lake build TraceThree`, the `#print axioms` audit, and
`lake env leanchecker --fresh TraceThree`, all successful with only
standard Lean foundations. The current Lean scope does not instantiate
the concrete P-192 elliptic curve or prove that its coordinate-wise
Frobenius action is the abstract additive map; the Sage receipt checks
those concrete coordinates. The dedicated [CI workflow](../../../../.github/workflows/p192-geometry-collector-lean.yml)
now builds, audits, and fresh-checks this module. The [run manifest](RUN.json)
freezes the proof and 256-column replay artifacts.
The first fresh-check attempt encountered a dangling local Mathlib-package
symlink. The [dependency repair log](../../lean/trace-three-lake-update-20261016.stderr)
records restoration from the pinned public revisions; the Lake manifest
hash stayed unchanged, and the repeated build, axiom audit, and fresh
kernel check passed.
