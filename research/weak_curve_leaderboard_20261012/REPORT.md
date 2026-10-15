# WeakCurves: source-linked names and weakness findings

The WeakCurves snapshot combines exact model identities, published
covering cases, and selected-subgroup attack criteria in one queryable atlas.
It provides a deterministic route from the two source-locked research
databases and the repository's named-curve and transfer registries to a
static catalog. The snapshot contains 814 prior exact models, 105 named
curves (101 standards-provenance curves, two challenge curves, and two
source-pinned pairing variants), 36
published odd-characteristic parameter cases, three exact prime-field
examples and twenty-six criterion entries. The WAM1 research queue adds
38 input- and witness-specific checks, giving the site 64 searchable attack
checks in total. Each evaluated finding has a
[WID1 identifier](IDENTIFIERS.md) distinct from its curve identity.

## What the leaderboard ranks

The rank is an **evidence order**, not an estimated security-bit score. An
exact DLP criterion names a model and selected group satisfying a published
method's stated arithmetic condition. A covering-model signal names an exact
model with a positive genus-three covering label. A candidate has an exact
screen, such as low embedding degree, while its map or comparative cost gate
remains open. Full-corpus rows retain unresolved and excluded labels. A
published parameter case is not an exact model and appears in its own view.

| View | Rows | Meaning |
| --- | ---: | --- |
| Findings | 29 | 2 prime-field DLP conditions, 24 covering-model signals, 1 pairing candidate, 2 Oakley GHS invariant candidates |
| Full exact corpus | 922 | 814 source-locked models, 105 named curves, 3 prime-field examples |
| Published cases | 36 | Iijima–Momose–Chao odd-characteristic `(2,…,2)` covering cases |
| Protocol cases | 1 | SIKEp434 auxiliary-torsion key-recovery setting |
| Large-field subset | 537 | 192–252-bit models; 384 unclassified in this index, 153 outside its domain |

The large-field branch screen refines the 384 unclassified degree-six models:
247 pass the cubic necessary condition, while 137 are excluded by the proved
depth condition. All 384 fail the recorded quadratic-model necessary
conditions. The remaining 153 degree-ten or degree-fourteen models are
outside the present genus-three domain. These statuses are extracted from
each model's linked `curve_property.family_checks_json` and rechecked by the
exporter’s count assertions; passing the cubic screen does not decide a
covering or subgroup advantage.

The 105 named curves are a tracked roster, not a census of all published
curves. The 103 inherited records have ICV1 identities; the two pairing
variants use source-pinned WCP1 identities. All retain standard or challenge aliases;
an unassessed row has no WID1 finding. WID1 is assigned only to a stated
condition, with its own selected scope and evidence grade.

## Oakley Groups 3 and 4

[RFC 2409](https://www.rfc-editor.org/rfc/rfc2409.html) names the binary
elliptic curves over `GF(2^155)` and `GF(2^185)` as IKE Groups 3 and 4.
[RFC 2412 Appendix E](https://www.rfc-editor.org/rfc/rfc2412.html) gives
their field polynomials, coefficients, points and group orders. The
existing [transfer registry](../../docs/curves/JACOBIAN_TRANSFERS.md)
calculates an exact GHS magic number of five for the degree-five towers.

| Standard name | Field tower | Full order | RFC generator order | Derived selected order | WID1 | Recorded route |
| --- | --- | --- | --- | --- | --- | --- |
| Oakley Group 3 | `2^155 = (2^31)^5` | `12r` | `4r` | `r = 3805993847215893016155463826195386266397436443` | `WID1:GHS:56f762afe456` | Exact magic number 5, Type II; formula genus 16 over `GF(2^31)` |
| Oakley Group 4 | `2^185 = (2^37)^5` | `4r` | `2r` | `r = 12259964326927110866866776214413170562013096250261263279` | `WID1:GHS:4df6967b110c` | Exact magic number 5, Type II; formula genus 16 over `GF(2^37)` |

The current binary transfer audit gives nonzero relative trace for both
models, hence Type II and the genus-16 formula. The previous genus-15
selectors identified different, erroneous findings; each replacement
record names its predecessor in `corrects`. The [GHS construction](https://www.kantmath.de/publications/papers/weil_descent.pdf)
supplies the genus dichotomy, while the model-specific type is checked by
the repository's binary-field transfer computation. The revised formula
costs are about 93.424 versus 90.388 rho-reference log2 units for Group 3,
and 99.829 versus 106.586 for Group 4. The chosen destination solver
hypothesis is false for both rows, so these numbers remain projections.

The full group cofactors 12 and 4 are **not** the selected DLP order:
the RFC generators have respective orders `4r` and `2r`, and
`[4]P` or `[2]P` has exact prime order `r`. The
[checked Sage replay](validation/oakley.sage) verified irreducibility of
both field polynomials, the RFC points on the specified curves, primality of
both `r` values, and the derived subgroup orders. It did not recompute the
RFC's full curve point counts. The conditional group-order step and the
integer factorizations are proved in [Lean](lean/OakleyOrders.lean).

The GHS findings are graded `exact_invariant_candidate`: the model-specific
invariant is exact, while the cover and injective map on the selected group
have not been constructed. The transfer registry's displayed cost numbers
are formula projections whose chosen Enge–Gaudry solver hypothesis is false
for these rows. They do not enter a measured-advantage ranking.

For a comparison of attack costs, the model, selected subgroup, explicit
transfer map, target distribution, construction cost and destination DLP cost
must all be fixed. None of the source corpus's 537 large-field rows currently
has a measured selected-subgroup advantage. Covering support is a geometric
entry, rather than a completed attack comparison.

## Genus-128 binary GHS cost-stage audit

The exact-invariant `GF(4)` routes for `c2pnb176w1` and `c2pnb272w1` have
magic number eight, Type II, and projected genus 128. An [independent
integer audit](src/bin/ghs_genus128_cost_audit.rs) of the transfer registry's
conditional Enge–Gaudry smoothness model optimizes the factor-base degree
over all 84 bounds examined. Both routes select degree bound 15. The count
uses `I_d=(1/d)∑_{e|d} μ(e)4^(d/e)` monic irreducibles of degree `d` and
counts degree-128 monic polynomials with all factors of degree at most 15.
The [frozen receipt](ghs_genus128_cost_audit.json) pins the registry, cost-model, and native producer hashes,
the exact polynomial counts, all degree-bound prices, and the two subgroup
parameters.
The [run report](reports/run-20261015T0854Z/REPORT.md),
[cost graph](reports/run-20261015T0854Z/ghs_cost_sweep.svg), and
[PDF](reports/run-20261015T0854Z/report.pdf) retain this model-stage evidence.
The [later source-lock refresh](reports/run-20261015T0929Z/REPORT.md)
updates the registry hashes and dependent site exports while preserving the
same model prices.
The [cover-source alignment](reports/run-20261015T0950Z/REPORT.md) updates
the decomposition index's covers-inventory digest after the clean-runner
validator exposed a third upstream source change.

| Source curve | Nominal half-split base | Model smoothness | Model relation queries | Model total log2 bit cost | Registry cost |
| --- | ---: | ---: | ---: | ---: | ---: |
| `c2pnb176w1` | 48,921,052 places | `2^-26.328277` | `2^51.872229` | 71.336033 | 71.336 |
| `c2pnb272w1` | 48,921,052 places | `2^-26.328277` | `2^51.872229` | 71.773568 | 71.774 |

The smoothness numerator is exactly
`1374288545819365998587940548442445318355469588992205247403709040067002`
out of `4^128` monic polynomials. The model assumes that half the base
irreducibles yield usable divisor places and that polynomial smoothness
predicts divisor relation yield on the eventual cover. These are precise
construction and measurement gates: build the cover, certify a nonzero
selected-subgroup image, enumerate its actual usable places, measure
ordinary-query smoothness and independent rank, then charge final linear
algebra and target descent. The current `ghs_descent.rs` constructs an
abstract higher-magic-number tower and an unsmoothed affine model at magic
number two; it does not yet construct a smooth magic-eight cover. That route
is the next construction task. The
[high-genus index-calculus literature](https://eprint.iacr.org/2011/098.pdf)
provides the relevant algorithm family, while the table above is a frozen
cost-model replay for these two exact source curves.

## Three exact prime-field models

The Rust enumerator checks every `x,y` in the named prime field for the
short-Weierstrass equation, adds the unique point at infinity, verifies the
nonzero discriminant, and computes the displayed point order by group
addition. PARI independently returned the same curve and point orders.

| Signal | Exact model and point | Curve order | Selected point order | Result |
| --- | --- | ---: | ---: | --- |
| Anomalous | `E/F_101: y²=x³+x+32`, `P=(4,10)` | 101 | 101 | `t=1`; [Smart's trace-one reduction](https://doi.org/10.1007/s001459900052) applies to the selected group |
| Smooth selected order | `E/F_101: y²=x³+x`, `P=(3,38)` | 100=`2²·5²` | 10=`2·5` | Exhaustive maximum point order 10; `E(F_101)≅Z/10×Z/10`. [Pohlig–Hellman](https://doi.org/10.1109/TIT.1978.1055817) applies to every selected cyclic subgroup |
| Embedding-degree candidate | `E/F_103: y²=x³+x`, `P=(18,44)` | 104=`8·13` | 13 | `103≡−1 (mod 13)`, hence `ord_13(103)=2`; pairing image and target-field cost are open |

For the anomalous model, the finite-field solution count is 100 affine
points plus infinity. For the smooth model it is 99 plus infinity. For the
pairing candidate it is 103 plus infinity. Their traces follow directly from
`#E(F_p)=p+1−t`. The latter model is supersingular: the trace is zero for
`p=103>3` and it has `j=1728`; the selected subgroup and embedding-degree
checks are exact, but this snapshot does not promote a pairing reduction
without a nondegenerate map and destination cost.

The examples are deliberately small, serving as exact recognition and
construction fixtures. They do not estimate the behavior of the 192–252-bit
corpus. The three exact point counts and the embedding-degree arithmetic are
also represented by kernel-checked Lean statements; the group-law and attack
theorems remain an explicit formalization bridge.

The Rust scanner checks the order of every affine point of the smooth model;
its maximum is 10. Together with group order 100 and the classification of
finite abelian groups of rank at most two, this gives `Z/10×Z/10`. PARI's
independent `ellgroup` result is `[10,10]`.

## Mechanism map

The [WAM1 attack map](attack_map.json) extends the existing criteria with
38 source-linked checks. Its entries state the required auxiliary inputs or
implementation behavior and the exact test needed to assess a curve or
protocol. They remain `unassessed_on_indexed_models` and add no WID1 result.
The site filters these checks by family and assessment state, while the curve
leaderboard filters recorded positive atlas checks.

The graph viewer draws finite-field trace classes and certified map edges
separately. It now filters explicit edges by degree, focuses one class, and
highlights a selected node's incident maps with degree labels. The source
graph remains 1,092 nodes in 676 classes with 232 certified directed maps;
the viewer does not infer missing maps from class peers.

The [hyperelliptic decomposition guide](HYPERELLIPTIC_DECOMPOSITION.md)
separates divisor factor-base decomposition, partial-relation recombination,
high-genus smoothness, Jacobian isogeny factorization, and the elliptic
quotient/Prym complement of a cover. The cover checker already constructs
same-field `H→E` maps for supported elliptic models. The genus-two source
factor is expected; the exact factor maps, selected-subgroup projections and
costed divisor relation pipeline are the further assessment gates.

The [machine-readable registry](criteria.json) specifies domain, recognition,
required witness, construction and primary source for twenty-six criteria. The
entries make four crucial separations.

The new `ic_relation_support_rank` entry tests a fixed-base index-calculus
candidate before treating a covering or character filter as evidence of a
usable relation method. An unordered positive \(m\)-sum base of \(B\) actual
usable classes covers at most \(\binom{B+m-1}{m}\) targets. Under a uniform
accepted-query law on \(A\) targets, its hit probability is at most that
support divided by \(A\), capped at one; verified relation yield and novel
rank remain separate measurements. The [P-192 Type I continuation](../p192_extension_cover_20261012/BRANCH_CONTINUATION.md)
certifies this distinction with 64 degree-three prime-field representations
among 64 accepted guided known-scalar queries, zero matches from a 16-class
three-sum base, and rank zero. Lean checks the exact 816-image support bound
and P-192 numerical capacity threshold. The criterion records feasibility
evidence and creates no new WID1 finding for an indexed curve.
The generic finite-image and accepted-set intersection caps are also proved
in [Lean](lean/RelationSupport.lean). Their [proof map](PROOF_MAP.md)
separates the combinatorial theorem from the curve-specific summation map
and the uniform-query assumption behind the probability bound.

The `ic_scalar_anchor_overlap` entry adds the equation gate: for a frozen
signed base, it compares nonzero-unknown-row support with independently
known-scalar anchor targets and checks whether repeated target keys have
different formal signed-point vectors. The exact [P-192 Type I audit](../p192_extension_cover_20261012/reports/run-20261015T1857Z/REPORT.md)
finds zero anchor overlap and zero distinct formal equations among repeated
rank-bearing targets. This criterion is indexed by the specified base and
query law; cover geometry alone does not assign a curve-wide finding.

The `ic_adaptive_base_rank_balance` entry checks dynamic factor bases.
Defining a new class from each known-scalar query makes a relation easy to
write, but also creates a new unknown logarithm. The exact [P-192 audit](../p192_extension_cover_20261012/reports/run-20261015T1926Z/REPORT.md)
constructs 64 such classes: its 64 rows have 64 private columns, leaving
all 16 original logarithms unconstrained. The full 8,704-point pair-policy
expansion has no duplicate or old-support decomposition. Lean proves the general unique
extension for every old-log assignment. The criterion asks for rank after
eliminating new columns and for an independent second decomposition or
known-log witness before classifying a dynamic collector as useful.

First, the [MOV reduction](https://doi.org/10.1145/103418.103434) uses a
Weil pairing, while [Frey–Rück](https://doi.org/10.1090/S0025-5718-1994-1218343-6)
uses a Tate-pairing variant. Both begin with the selected order `r` and
embedding degree `k=ord_r(q)` when `gcd(r,q)=1`, but a small `k` alone does
not certify a nondegenerate pairing image or favorable finite-field DLP cost.
The F_103 record therefore remains a candidate.

Second, binary [GHS Weil descent](https://www.iacr.org/cryptodb/data/paper.php?pubkey=14115)
and the [odd-characteristic covering classification](https://eprint.iacr.org/2015/805.pdf)
have different domains. The present published-case index covers odd
characteristic and prime extension degrees 2, 3, 5 and 7. It records 36
numbered parameter cases, including five minimum-genus patterns. The
existing genus-three model index works over `F_(p^6)` with `q=p²`, using
the cubic norm-one and nonsplit quadratic branches separately. Its
`v₂(f_π)≥2` cubic rule is necessary within that family, not a universal
criterion for weak curves. The p41, p43 and p47 cubic diagnostics retain
396, 376 and 424 ordinary high-depth zero rows respectively, so the
depth screen cannot serve as a sufficient classifier even there.

Third, the atlas now distinguishes three additional routes. Direct
[extension-field decomposition](https://eprint.iacr.org/2010/157) works on
the elliptic group itself; [cover-and-decomposition](https://eprint.iacr.org/2011/020)
combines a lower-field cover, norm–pullback subgroup transfer, and a
destination relation method; and a
[cubic Jacobian isogeny](https://arxiv.org/abs/2012.07173) can transfer a
prime-order group without a literal covering map. Their published methods
are present as criteria. Exact atlas joins and complete selected-subgroup
measurements remain open; the repository's separate Joux–Vitse replay is
tracked in the coverage audit.
The [per-curve decomposition index](decomposition_properties.json) now
screens both sextic intermediate-field towers and the quadratic sieve
separately. The [published 149-bit case](published_cover_cases.json) retains
the source's exact tower, selected order, factor base and phase reports as a
published campaign outside the 922 local atlas identities.
The [coverage audit](COVERAGE_AUDIT.md) maps these methods, the 36
prime-degree cases, and same-field hyperelliptic covers to their distinct
construction and measurement gates.

Fourth, the ambient curve order, a selected point order, and an exact
transfer all answer different questions. A curve with a smooth ambient
order still needs a selected point and target membership to support a
selected-group claim. A model over an extension field only yields a
smaller-field route if the selected points descend. These distinctions
prevent a cover, field degree, or factorization from silently becoming a
measured advantage.

Fifth, the graph and filters separate finite-field isogeny-class membership
from an explicit route. Equal field size and trace put elliptic curves in the
same isogeny class by Tate's theorem. The [graph export](ISOGENY_GRAPH.md)
contains 1,092 indexed curves in 676 classes and 232 certified directed
isogeny maps from the P-192/P-224 walk records. An atlas node receives an
isogeny-transfer map tag only when it is an endpoint of a certified edge;
the 319 atlas nodes in multi-model classes are not thereby assigned a map.
The GLV filter currently finds seven exact prime-field models with a rational
`j=0` or `j=1728` automorphism screen. These are structural leads until the
selected-subgroup endomorphism action and method cost are checked. GLS has a
published recognition entry and no exact atlas-model match in this export.

Sixth, [Ethereum's exact BN_SNARK1 variant](https://eips.ethereum.org/EIPS/eip-197)
and [BLS12-381](https://datatracker.ietf.org/doc/draft-irtf-cfrg-pairing-friendly-curves/)
have separately recorded order-`r` G1 and twist G2 generators and GT output
fields. The [arithmetic replay](validate_pairing_roster.py) checks each
generator's curve equation and `[r]P=O`, the G1 cofactor relation, and
`ord_r(p)=12`. A pairing-friendly design supports its intended pairing
protocols; DDH, static DH, subgroup-input exposure, and GT DLP each have
their own witness and cost gate. [Hitt's minimal-field distinction](https://eprint.iacr.org/2006/415)
is formalized at the group-order level in `lean/MinimalPairingField.lean`.
The [SIKEp434 protocol row](protocol_cases.json) records the
[Castryck–Decru](https://eprint.iacr.org/2022/975) key-recovery result under
its auxiliary-torsion-image input, separate from the ordinary isogeny graph.

## Construction and search method

1. Fix an input law and an exact field representation. For prime fields,
   enumerate or construct nonsingular models and count points. For extension
   fields, select a published covering family or a specific norm-one or
   quadratic-torus branch.
2. Record the exact model, group order, trace, chosen point and order.
   Compute arithmetic screens: trace one, smooth selected order, embedding
   degree, supersingularity, CM depth, or proper-subfield containment.
3. For a pairing or covering, construct the auxiliary field/curve and map.
   Check genus, smoothness, kernel intersection and selected-point image.
4. Preserve failed constructions and distinct class/model labels. Run a
   destination DLP cost study only for a fixed selected subgroup and a
   specified target law. Pair it with the same-point reference under the
   repository's measurement contract.

The next import should normalize equations and finite-field realizations
for the 36 published cases, then add exact subgroup maps and pairing
certificates. The most useful large-field experiment is a frozen 192–252-bit
cohort with selected prime-order subgroups and a costed map check, rather
than a field-size-only ranking.

## Sources and reproducibility

The original model database is
[`../weak_curve_database_20261012/weak_curves.sqlite`](../weak_curve_database_20261012/weak_curves.sqlite),
and the numbered case index is
[`../weak_curve_case_index_20261012/published_cases.sqlite`](../weak_curve_case_index_20261012/published_cases.sqlite).
Both SHA-256 digests are embedded in `snapshot.json`; the exporter checks
the joined row totals on every build. The Rust source records the exact
prime-field enumeration and point addition. The [validation record](VALIDATION.md)
contains commands, independent PARI output, and the Lean audit.
