# Extension-field and Jacobian-transfer coverage audit

The atlas now names thirteen mechanisms, including direct extension-field
index calculus, cover-and-decomposition, and cubic Weil-restriction Jacobian
isogenies. Its 920 exact curve rows are **not** a census of every curve to
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
