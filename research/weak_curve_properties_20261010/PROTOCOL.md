# Certified structural properties for the curve catalogue

Recorded before property generation on 10 October 2026.

Continue categorizing the existing 814 exact curve models by certifiable
properties. This round adds a native structural-property extractor, an
independent PARI replay, a searchable property browser, a report with exact
scope and evidence, a diagram, an inventory graph and a reviewed PDF. Extend
the existing publication PR. Preserve the earlier sealed catalogue study and
all parent arithmetic evidence.

The sole model input is the frozen `docs/curves/weak-families/models.json`,
SHA-256 `9738f6bd6c016732f88e4e7979ac09e9aa6fd16eba96220ec8781c0017b6055e`.
Its 814 rows retain their inventory laws: 277 p7 route/conversion/control models,
25 larger-degree constructed controls, and 512 independent large-field sources.
Family membership, class-support decisions and bounded-search statuses are
copied unchanged. No inventory count becomes a uniform-class density estimate.

Certify ordinary/supersingular status, ambient-order small-prime valuations,
Frobenius discriminant and its local conductor depth at two, splitting of two
in the CM field, rational 2-torsion and full rational 4-torsion, the minimal
Frobenius degree of j, the displayed coefficients' field, and the smallest
subfield admitting a K-isomorphic model. Unknown full endomorphism orders,
selected subgroups, point containment and comparative costs remain explicit.

The conductor calculation uses the two-adic fundamental-discriminant cases.
An independent PARI method searches integral (pi-a)/2^k using trace and norm;
it covers every ordinary model without factoring its odd discriminant. Native
polynomial-basis arithmetic checks irreducibility, j, Frobenius fixed fields,
twist power classes and the root-difference criterion for full 4-torsion.
PARI independently checks those properties, returns all rational 2-roots,
and verifies full 4-torsion using the sextic division polynomial and y squares.

For an ordinary model, one rational nonidentity 2-point certifies the ordinary
2-volcano floor. Full 2- and 4-torsion give upper bounds on the endomorphism
conductor's 2-depth, and a collapsed bound certifies that local level. Full
endomorphism conductors remain unknown. The report supplies the order-basis
and separable-torsion proof, and keeps exact levels distinct from bounds.

Primitive j, displayed coefficients and K-isomorphic-model descent are separate
properties. Proper-subfield availability alone does not classify the selected
points or subgroup as contained in that subfield. Passing a necessary condition
does not promote unresolved family membership or class support.

Ship portable identity certificates and exact integer-polynomial expansions
for the load-bearing cleared identities. Use boundary controls for ordinary
versus supersingular inputs, discriminant two-adic cases, exceptional j values,
and j descent with and without isomorphic-model descent. Mutated inputs and
root witnesses must be rejected.

The quantitative graph describes this frozen inventory by cohort and property;
its counts are exact and have no sampling-error interpretation. Existing IC
ratio graphs and prior census values retain their numbers, because this round
supplies structural classification rather than a new paired timing measurement.
All arithmetic, research data generation and verification use native Rust and
compiled PARI. Browser checks use Node for presentation only. Local runtime
metadata is retained; no CPU speedup is inferred from these correctness runs.
