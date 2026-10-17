# P-192 Type I line support and projected point base

**Result.** An explicit inverse from the smooth genus-three plane quartic to
its Type I cover converts rational line-section points into independently
checkable elliptic-quotient classes. From the frozen 141 rational points in
the 160 fixed-origin residual divisors, all 9,870 pair-lines are distinct.
Exactly 4,892 line residual quadratics split over the P-192 prime field;
they yield 9,784 new rational points. The 141 original and 9,784 new points
are distinct, and cofactor projection gives 9,925 distinct nonzero first-
quotient (r)-torsion classes, also distinct after sign folding. The 141
original classes have 9,870 distinct unordered two-point sums.

The fully enlarged 9,925-point base has at most
\(\binom{9927}{3}=162993977475\) unordered three-sum targets. Against the
previously certified accepted nonzero scalar-set lower bound
\(L=1569275433846670190958947353987675545133315652035359680560\),
an ideal uniform accepted scalar query hits this whole support with
probability below \(2^{-152}\). This is a capacity bound for the stated
base and three-summand policy; the ordinary-query yield and complete
one-target cost are separate measurements.

The exact 141-class line-derived base supports a frozen known-scalar query
comparison. Its first 16 classes have 816 distinct three-sum targets; 206
pass the two-fiber geometry predicate. In 64 queries per law, the uniform
and geometry-conditioned laws both have zero exact three-sum hits and
incidence rank zero. The uniform arm has 11 geometry-admitted queries; the
conditioned arm uses 259 scalar trials to obtain 64 admitted queries.
The known equation is \(T_a=a(-2G)\) in each arm. The earlier
[fixed-origin cubic-section certificate](../../type_i_fixed_origin_reduction.json)
and [guided-query certificate](../../type_i_guided_negative_fixed_origin_reduction.json)
supply degree-three Jacobian representatives for the 16 original and 64
conditioned scalar classes; the 141 line-derived
point classes are checked here by their point differences and injective
odd-order quotient coordinates.

For all 9,925 selected line points, a separate [fixed-origin screen](../../line_reduced_source.py)
also constructs a unique effective degree-three representative of the
*source* class \([P-A]\), with \(A\) the exact
[quotient-base cover anchor](../../line_anchor_plane.py). It chooses rational points \(B,C\) whose line has
an irreducible quadratic residual and intersects each of \(B,C\) simply,
then fixes \(D_0=A+B+C\) and \(D_P=P+B+C\). Only \(B,C\) from the selected
set lie on that line, so no \(D_P\) is a canonical line section; for a
smooth nonhyperelliptic plane quartic, \(h^0(D_P)=1\). The receipt gives the
exact anchor and line, and the Lean proof checks
\([D_P-D_0]=[P-A]\) algebraically. The projected \(r\)-torsion class is
\(h[P-A]\), so reducing this cofactor multiple itself to degree three
remains a separate Jacobian-arithmetic step.

## Geometry and inverse map

The [native GP screen](../../line_support_screen.gp), invoked by its
[deterministic wrapper](../../produce_line_support.py), starts from the
frozen 160-class support and plane-quartic coefficients. For each pair
of the 141 selected rational points, it substitutes their line into the
quartic, divides by the two known roots, and factors the degree-two
residual. All 9,870 residuals are squarefree and of degree two: 4,978
are nonsplit and 4,892 split. Their [point list](../../type_i_line_new_points.tsv)
and [pair identities](../../type_i_line_new_point_pairs.tsv) contain
9,784 fresh points. The [independent Python replay](../../verify_line_support_receipt.py)
rebuilds all 9,870 line polynomials from the projective quartic coefficients,
checks each residual discriminant, every rational root and curve equation,
and the absence of point or line-key reuse. The source 141 points are in
[their ordered list](../../type_i_line_base_points.tsv).

The [inverse implementation](../../plane_point_inverse.py) uses the three
conjugate quartic right-hand sides \(f_i(X)\) and the canonical plane forms
\(A_i\). It finds \(X\) as the common linear factor of
\(A_0^2f_0-A_1^2f_1\) and \(A_0^2f_0-A_2^2f_2\). With \(N(X)\) the
six-factor norm, the scale
\(\lambda=A_0f_0(X)/(A_1A_2N(X))\) gives
\(y_i=1/(\lambda A_i)\). Each candidate is checked against all three
quartic equations, the norm, the descended cover model, prime-field
Frobenius, and the plane-coordinate roundtrip. The inverse has exceptional
denominators; all 9,925 selected points pass the nonexceptional checks.

The [Sage projection screen](../../line_projection_screen.py) uses the
accepted repository launcher and the frozen [runtime receipt](../../type_i_line_sage_runtime_info.json).
It first inverts 16 previously certified plane points and recovers their
exact elliptic-quotient classes. It checks all three quotient components
for the 141 original points and the first 64 new points, and checks the
first component and its (r)-torsion order for every new point. The
[independent full Sage replay](../../verify_line_projection.py) reconstructs
the inverse without importing the producer, repeats every projected key
and pair sum, and compares the ordered key digests in the
[projection receipt](../../type_i_line_projection.json). The
[receipt verifier](../../verify_line_projection_receipt.py) checks source,
input, runtime, and replay-log custody.

| Exact frozen screen | Result |
| --- | ---: |
| Original rational points / distinct pair-lines | 141 / 9,870 |
| Nonsplit / split line residuals | 4,978 / 4,892 |
| New rational line points / total selected points | 9,784 / 9,925 |
| Distinct projected base classes / pair sums | 141 / 9,870 |
| New first-quotient classes, including sign folding | 9,784 |
| New-to-base signed collisions / new-to-new signed collisions | 0 / 0 |
| Exact old-base controls / new three-component controls | 16 / 64 |

## Known-scalar geometry-guided collector

The [collector](../../line_guided_collector.py) projects each of the 141
original line-support points into the first P-192 elliptic quotient after
subtracting one frozen cover-point anchor and applying the exact cofactor.
All three quotient coordinates agree for these 141 classes, and their ordered
keys match the independently replayed projection digest. It builds the
10,011 unordered pair sums with repetition and probes every exact
three-sum for each frozen query.

Both arms use seed `0x19220261023`. The uniform arm takes the first 64
nonzero scalar draws. The conditioned arm takes the first 64 draws for which
the inverse quotient fibers of \(aG\) and \(2aG\) have nonzero square
quartic right-hand sides in \(\mathbb F_{p^3}\); the prior receipt records
all 259 trials and the 195 rejections. Both targets are
\(T_a=a(-2G)\), so each prospective match yields a known-scalar equation.
The first-16 exact support census classifies a triple target \(T\) by
\(P=-(2^{-1}\bmod r)T\), for which \(T=-2P\), and tests the same two fibers.

| Frozen screen | Uniform scalars | Geometry-conditioned scalars |
| --- | ---: | ---: |
| Queries | 64 | 64 |
| Scalar trials charged | 64 | 259 |
| Geometry-admitted queries | 11 | 64 |
| Exact three-sum matches on 141 classes | 0 | 0 |
| Quotient incidence rank | 0 | 0 |

The [independent prime-field group replay](../../verify_line_guided_collector.py)
checks the 141 archived projected points, all pair keys, both scalar lists,
every queried elliptic target, matches, and row rank without Sage or the
collector's elliptic arithmetic. It relies on the earlier independent
Sage replay for the two-fiber character tests; a further
[independent character replay](../../verify_line_guided_character.py)
recomputes the 206/816 support classification and both sample admission
counts directly from the frozen projected points. The conditional support
fraction `206/816` is specific to the first 16 line-derived classes;
the full 141-class conditioned support was not enumerated. Neither arm
produced a useful row in this bounded comparison.

The first independent projection replay reached the full 9,784-point loop
but could not write its receipt when the SSD returned `ENOSPC`; its
[failed log](../../type_i_line_projection_verification.enospc.log) is retained.
The first rebuilt Lean `LineSupport` check stopped on the same disk error;
its [failed log](../../lean/line-support-build-20261015.enospc.log) is retained.
The reruns use the same frozen sources and inputs with regenerable build
output on the system volume.

## Relation consequence and proof map

Each of the 4,892 split line sections has two new rational points. Under
the frozen point-column policy, every such line equation therefore has
private new columns. The [Lean `LineSupport` theorem](../../lean/LineSupport.lean)
proves that a family of rows with one private column per row extends every
assignment of old columns; adding these line rows alone does not determine
an old logarithm. It also checks the exact capacity and hit-bound integers.
Further function-divisor relations or projected class identifications would
change that matrix, so they are the next constructive search target.

The [Lean inverse module](../../lean/LinePointInverse.lean) proves the
six-factor quartic product identity and uniqueness of the scale from the
squared and norm equations. Both modules build under the pinned Lean 4.34.1
and Mathlib revision, replay in a fresh checker, and have audited axioms.
The concrete polynomial gcds, curve and quotient maps, and 9,925 point
instances are checked by GP, Python, and Sage; a full concrete Lean
formalization of that geometric bridge and the Riemann--Roch nonspeciality
implication remains open.

| Mathematical statement | Lean declaration | Boundary |
| --- | --- | --- |
| Six-factor norm identity and scale uniqueness | `P192ExtensionCover.three_quartics_product_is_norm_square`, `P192ExtensionCover.scale_fixed_by_square_and_norm` | Explicit algebraic hypotheses; exact P-192 inverse instances are Sage-replayed. |
| Fixed-origin point-divisor class identity | `P192ExtensionCover.fixed_origin_point_divisor_identity` | Additive divisor identity; quartic nonspeciality and RR bridge are checked or proved in the written geometry. |
| Private columns leave old-log assignments unconstrained | `P192ExtensionCover.line_private_pairs_extend` | Assumes one private column per row; exact point distinctness is independently replayed. |
| Three-sum capacity and accepted-set inequality | `P192ExtensionCover.p192_line_point_three_sum_capacity`, `P192ExtensionCover.p192_line_point_hit_bound_certificate` | Integer certificates; subgroup size and accepted-set lower bound come from earlier receipts. |

The [run manifest](RUN.json) binds sources, inputs, outputs, logs, and
proof receipts by SHA-256. The proposal keeps `candidate_id: null`: a
frozen ordinary known-scalar query experiment, projected useful rank,
linear algebra, target descent, and complete one-target accounting are
the next gates. The recorded local elapsed values are correctness-run
diagnostics under the stated host, not isolated performance comparisons.

## Reproduction

From the repository root, with the accepted checked Sage launcher:

```sh
python3 research/p192_extension_cover_20261012/produce_line_support.py
python3 research/p192_extension_cover_20261012/verify_line_support_receipt.py
/Volumes/SSD990/cryptanalysis/sage --runtime-info > research/p192_extension_cover_20261012/type_i_line_anchor_runtime_info.json
/Volumes/SSD990/cryptanalysis/sage -python research/p192_extension_cover_20261012/line_anchor_plane.py
python3 research/p192_extension_cover_20261012/line_reduced_source.py
/Volumes/SSD990/cryptanalysis/sage --runtime-info > research/p192_extension_cover_20261012/type_i_line_sage_runtime_info.json
/Volumes/SSD990/cryptanalysis/sage -python research/p192_extension_cover_20261012/line_projection_screen.py
/Volumes/SSD990/cryptanalysis/sage --runtime-info > research/p192_extension_cover_20261012/type_i_line_verify_runtime_info.json
/Volumes/SSD990/cryptanalysis/sage -python research/p192_extension_cover_20261012/verify_line_projection.py
python3 research/p192_extension_cover_20261012/verify_line_projection_receipt.py
/Volumes/SSD990/cryptanalysis/sage --runtime-info > research/p192_extension_cover_20261012/type_i_line_collector_runtime_info.json
/Volumes/SSD990/cryptanalysis/sage -python research/p192_extension_cover_20261012/line_guided_collector.py
python3 research/p192_extension_cover_20261012/verify_line_guided_collector.py
```
