# Extension-field and Jacobian-transfer coverage audit

The atlas now names twenty-three criteria, including direct extension-field
index calculus, cover-and-decomposition, and cubic Weil-restriction Jacobian
isogenies. Its 922 exact curve rows are **not** a census of every curve to
which those methods apply. The `criteria.json` entries describe recognition
and witness requirements; they do not add WID1 findings until an exact model
and claim scope are checked.

| Route | Present evidence | Curve-level gap |
| --- | --- | --- |
| Odd-characteristic `(2,…,2)` covering for prime extension degree `d=2,3,5,7` | [36 published parameter cases](../weak_curve_case_index_20261012/REPORT.md), including five with cover genus `d` | Instantiate the equations and square-class conditions for each exact model; verify smoothness, genus, map, Jacobian and subgroup transport. The parameter cases do not classify the 814 local models. |
| Direct extension-field decomposition on `E(F_(q^d))` | [Transfer registry](../../docs/curves/JACOBIAN_TRANSFERS.md) has nine `extension_index_calculus` formula rows; [Joux–Vitse](https://eprint.iacr.org/2010/157) obtained degree-five relations and an oracle-assisted static-DH application | Import executable factor-base, relation, rank, linear-algebra, and target records. This route needs no cover. |
| Cover-and-decomposition over composite-degree odd extensions | [Joux–Vitse](https://eprint.iacr.org/2011/020) gives the method; the repository has a [separate constructed-instance replay](../jv_cover_e2e_30bit_20261009/RESULT.md) | Join exact model and subgroup identities to current atlas rows; rerun selected targets with complete cost and paired reference receipts under the current contract. |
| Cubic Weil-restriction isogeny to a genus-three Jacobian | [Tian](https://arxiv.org/abs/2012.07173) gives a prime-order-group route over `F_(q^3)` without a literal covering map | Add its exact isogeny and kernel screen, destination curve and selected-subgroup transport to the atlas. No local curve row currently carries this criterion. |
| Classical binary GHS | Oakley Groups 3 and 4 have exact tower magic-number candidates and separate WID1 markers | Construct the claimed smaller-field curve and map, then check the chosen prime-order subgroup and solver cost. |
| Generalized GHS and isogeny extension | The [transfer registry](../../docs/curves/JACOBIAN_TRANSFERS.md) has isogeny-extension formula rows; it explicitly does not enumerate Hess's generalized Artin–Schreier covers | Search exact representatives and function fields rather than treating formula counts as constructed maps. |
| Same-field hyperelliptic covers | [COVERS.md](../../docs/curves/COVERS.md) gives verified degree-two genus-two covers in odd characteristic and degree-three genus-three covers for ordinary binary models; the transfer registry contains 312 `same_field_cover` route rows | Surface these geometric facts on curve pages, with field-of-definition and cost boundaries. A same-field cover by itself does not produce a smaller-field attack or a WID1 weakness finding. |

The present source-locked exact-model corpus has 661 degree-six, 77
degree-ten, and 76 degree-fourteen odd-characteristic extension models. The
published-case index is organized by the *relative* prime degree `d` over a
chosen subfield: `F_(p^6)/F_(p^2)` has `d=3`, while degrees ten and fourteen
have relative degrees five and seven over `F_(p^2)`. Matching the number `d`
does not match a paper's coefficient, ramification, or square-class
conditions. The 103 named-curve roster includes prime-field and binary-field
standards, but no named odd-characteristic prime-degree extension standard.

The 2015 classification covers its stated odd-characteristic
`(2,…,2)`-covering family for prime degrees `2,3,5,7`; its [abstract](https://eprint.iacr.org/2015/805)
reports the large-genus issue at prime degrees at least eleven. It is not a
classification of every possible transfer. The [earlier classification
without the isogeny condition](https://eprint.iacr.org/2009/613),
[cubic construction and hyperellipticity test](https://eprint.iacr.org/2008/215),
and Tian's direct isogeny give distinct search branches. For a genus-three
destination, record whether it is hyperelliptic: a Jacobian may require a
different arithmetic and relation solver when it is not.

## Promotion procedure for one exact curve

1. Fix the full field representation, base-subfield embedding, curve, selected
   generator order `r`, and target law. Enumerate every relevant tower, rather
   than identifying a route from the absolute field degree alone.
2. Run route-specific screens: published ramification and square-class tests,
   binary magic number, rational torsion, extension decomposition support,
   or a Weil-restriction isogeny condition. Preserve rejected and unresolved
   candidates with their exact condition.
3. Construct the actual destination and map. Record equations, field of
   definition, smoothness, genus, hyperelliptic status, map degree or isogeny
   kernel, and exceptional points. For `φ:C→E` of degree `m`, the relation
   `φ_*φ^*=[m]` makes pullback injective on an order-`r` subgroup when
   `gcd(m,r)=1`; otherwise inspect the kernel intersection directly. A direct
   isogeny must pass the analogous selected-subgroup check.
4. Verify the image of the actual generator and target, divisor arithmetic,
   and scalar compatibility. Then measure factor-base construction, ordinary
   relation yield and rank, final linear algebra, target descent, and replay.
   Compare one previously unseen target with a matched source-group reference
   and keep the complete resource and failure record.
5. Issue a WID1 for a precise model or selected-subgroup condition and raise
   its evidence grade only when the matching witness is present. A universal
   same-field cover stays a geometric property; it is not automatically a
   weakness finding.

The repository's [cover and Jacobian skill](../../.agents/skills/ecdlp-surfaces/SKILL.md)
specifies the map and arithmetic checks. The [proof map](PROOF_MAP.md) lists
the Lean status of the atlas's mathematical claims. This audit adds no new
mathematical theorem; it records published route families and the exact
unproved or unmeasured bridge for each.

## Attack families and assessment layers still missing

The twenty-three-criterion registry is a method index, not a complete attack audit.
In particular, its MOV and Frey–Rück entries share an embedding-degree screen,
and several geometric entries have no constructed map on an atlas model. The
following additions need separate records for an exact curve, selected subgroup,
attack inputs, and measured or modeled cost. A structural candidate stays a
candidate until its stated gate passes.

| Proposed check | Exact recognition and evidence gate | Present gap and source |
| --- | --- | --- |
| Twist and invalid-input exposure | Count and factor the curve and its relevant twists; record input format, accepted point set, subgroup validation, cofactor handling, secret-scalar policy, and whether a multiplication result is observable. Test the implementation or protocol path against the actual small-order or invalid inputs. | `smooth_selected_order` covers only the selected subgroup on the listed model. Twist and invalid-curve attacks depend on input handling; [RFC 7748](https://www.rfc-editor.org/rfc/rfc7748) discusses low-order inputs and the all-zero check, and [Bernstein–Lange](https://eprint.iacr.org/2024/1265) analyze twist, small-subgroup, and invalid-curve exposure. |
| Cheon auxiliary-input DLP | Factor `r−1` and `r+1` for the selected prime subgroup order `r`, then record the *actual* supplied powers of one secret, such as `G,[x]G,[x^d]G` with `d | r−1`; estimate the best attack for that input set. | The factors alone never establish the attack precondition. [Cheon](https://doi.org/10.1007/11761679_1) gives the auxiliary-input algorithm; [pairing-curve analysis](https://eprint.iacr.org/2020/351) applies it to structured-reference-string power sequences. |
| Trace-zero subgroup index calculus | Specify `E/F_q`, the extension `F_(q^n)`, trace map, and selected point's membership in its kernel. Build the trace-zero factor base and record natural-query relations, rank, linear algebra, descent, and one-target cost. | `extension_index_calculus` analyzes decomposition on the whole elliptic group, while [Gorla–Massierer](https://eprint.iacr.org/2014/318) analyze a distinct algorithm on the trace-zero variety. The generator's trace-zero membership is a required gate. |
| Pairing destination cost | The new `minimal_pairing_field` and `pairing_target_field_dlp` entries record `d=ord_r(p)` and exact G1/G2/GT roles. Next, certify the pairing image, field tower, finite-field DLP algorithm, and transfer plus individual-log cost. | The two named pairing variants have checked `d=k=12` and checked G1/G2 generators. A calibrated GT discrete-log estimate remains open. [Hitt](https://eprint.iacr.org/2006/415), [Tower NFS](https://eprint.iacr.org/2015/505), and [pairing-curve cost estimates](https://eprint.iacr.org/2019/1371) distinguish the parameters and cost models. |
| Other efficient automorphisms and endomorphisms | Construct the endomorphism, verify its action and orbit size on the selected subgroup, and compare the resulting rho cost on matched inputs. Include binary Koblitz Frobenius and endomorphisms beyond the `j=0,1728` prime-field screens. | Current GLV and GLS entries are narrow structural screens. [Koblitz-curve analysis](https://eprint.iacr.org/2016/603) discusses Frobenius rho speedups; [SafeCurves' discriminant analysis](https://safecurves.cr.yp.to/disc.html) explains the small-CM-discriminant signal and its limits. |
| Prime-field summation-polynomial index calculus | Fix a prime-field factor base and exact decomposition system, then measure ordinary-query relation yield, rank, final linear algebra, target descent, and matched one-target cost. | `extension_index_calculus` is scoped to `d>1`; [Amadori–Pintore–Sala](https://eprint.iacr.org/2017/609) study the prime-field branch. A proposed decomposition formula supplies a research route, while a weakness finding needs the complete pipeline. |
| Restricted secrets, nonce leakage, and faults | Attach a separate application or implementation record: secret-key distribution, signature samples and nonce relation, observable side channel or fault model, and a reproducible key-recovery check. | These conditions concern use of a curve rather than its ECDLP subgroup. [Sun–Espitau–Tibouchi–Abe](https://eprint.iacr.org/2021/455) analyze partial ECDSA nonce leakage, and [Neves–Tibouchi](https://eprint.iacr.org/2015/1233) analyze fault-induced degenerate-curve behavior. |

### Promotion order

1. Add exact twist orders, selected subgroup orders, and input-validation
   metadata for the named-curve roster. Publish separate **curve structure**
   and **protocol exposure** badges, with a source and exact acceptance test.
2. Add Cheon and trace-zero screens only when the auxiliary-input record or
   trace-zero membership is known. Store `unassessed` for missing context.
3. Replace pairing embedding-degree rankings with reproducible finite-field
   cost estimates and an explicit pairing witness. Expand GLV/GLS screens to
   exact subgroup actions and include binary Frobenius.
4. Connect existing GHS, cover, Jacobian, and isogeny candidates to exact
   maps and complete cost records using the procedure above. The current
   graph's verified edges and the atlas's cover candidates are separate
   evidence types.

Future mathematical recognition results follow the [proof map](PROOF_MAP.md)
and the repository's Lean requirement before promotion to theorem status.

## Security goals and protocol-specific records added in snapshot 03

The new DDH criterion asks for a nondegenerate pairing on the **challenge
inputs**; the F103 supersingular example is a structural candidate pending an
explicit distortion-map pairing witness. Static DH asks for a reused-secret
oracle and therefore has no curve-only positive tag. Both tests distinguish
their security goal from elliptic DLP. The named Ethereum BN_SNARK1 and
BLS12-381 rows carry exact G1/G2/GT roles and subgroup parameters; membership
misuse is an implementation claim that needs a parser and accepted input.
The SIKEp434 case is a separate protocol record because its published
Castryck–Decru attack relies on auxiliary torsion images in SIDH public keys,
not an ordinary finite-field isogeny edge. See [Galbraith–Rotger](https://eprint.iacr.org/2004/070),
[Granger](https://eprint.iacr.org/2010/177), the [CFRG pairing draft](https://datatracker.ietf.org/doc/draft-irtf-cfrg-pairing-friendly-curves/),
[EIP-197](https://eips.ethereum.org/EIPS/eip-197), and
[Castryck–Decru](https://eprint.iacr.org/2022/975).
