# WeakCurves: evidence-graded criteria and first leaderboard

The first WeakCurves snapshot combines exact model identities, published
covering cases, and selected-subgroup attack criteria in one queryable atlas.
It provides a deterministic route from the two source-locked research
databases to a public-facing static catalog. The snapshot contains 814
previous exact models, 36 published odd-characteristic parameter cases,
three new exact prime-field examples and ten criterion entries.

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
| Exact signals | 27 | 2 prime-field DLP conditions, 24 covering-model signals, 1 pairing candidate |
| Full exact corpus | 817 | 814 source-locked models plus 3 new prime-field examples |
| Published cases | 36 | Iijima–Momose–Chao odd-characteristic `(2,…,2)` covering cases |
| Large-field subset | 537 | 192–252-bit models; 384 unclassified in this index, 153 outside its domain |

The large-field branch screen refines the 384 unclassified degree-six models:
247 pass the cubic necessary condition, while 137 are excluded by the proved
depth condition. All 384 fail the recorded quadratic-model necessary
conditions. The remaining 153 degree-ten or degree-fourteen models are
outside the present genus-three domain. These statuses are extracted from
each model's linked `curve_property.family_checks_json` and rechecked by the
exporter's count assertions; passing the cubic screen does not decide a
covering or subgroup advantage.

For a comparison of attack costs, the model, selected subgroup, explicit
transfer map, target distribution, construction cost and destination DLP cost
must all be fixed. None of the source corpus's 537 large-field rows currently
has a measured selected-subgroup advantage. Covering support is a geometric
entry, rather than a completed attack comparison.

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

The [machine-readable registry](criteria.json) specifies domain, recognition,
required witness, construction and primary source for ten criteria. The
entries make three crucial separations.

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

Third, the ambient curve order, a selected point order, and an exact
transfer all answer different questions. A curve with a smooth ambient
order still needs a selected point and target membership to support a
selected-group claim. A model over an extension field only yields a
smaller-field route if the selected points descend. These distinctions
prevent a cover, field degree, or factorization from silently becoming a
measured advantage.

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
