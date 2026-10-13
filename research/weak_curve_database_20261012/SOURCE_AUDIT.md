# Primary-source coverage and case-import backlog

The catalogue inherits the prior family registry and records four
additional primary sources in [sources.json](sources.json). This audit
checked the cited abstracts and earlier reviewed selected sections.
It records exactly which cases are already represented and which
published case tables still need equation-by-equation normalization.

| Primary source | Established contribution used here | Database action |
| --- | --- | --- |
| [Momose–Chao 2009/2013](https://eprint.iacr.org/2009/236) | Classification and density of cubic-extension elliptic curves with hyperelliptic and non-hyperelliptic weak coverings | Existing cubic and Type I/II family entries retain attribution; exact model/class support is narrower and separately evidenced |
| [Iijima–Momose–Chao 2013/2015](https://eprint.iacr.org/2013/487) | Complete odd-characteristic list under its (2,...,2) covering and isogeny condition; abstract reports ten previously known plus 17 additional classes across elliptic and genus-2/3 source curves | Import each equation as a published case after normalizing characteristic, genus, cover group, twist and equivalence; do not add 17 directly to the 21-entry count |
| [Tian 2020/2023](https://arxiv.org/abs/2012.07173) | An F_q-rational odd-prime isogeny from cubic Weil restriction to a genus-three Jacobian, including selected prime-order cases | Add a separate polarized-isogeny recipe; future records need actual kernel, polarization and map evidence |
| [Freeman–Scott–Teske 2006/2009](https://eprint.iacr.org/2006/372) | Framework for pairing-friendly curve constructions and exact embedding-degree questions | Add selected-subgroup pairing and supersingular construction recipes; embedding degree is stored as a structural condition with a destination-field gate |
| [Smith 2008/2009](https://arxiv.org/abs/0806.2995) | Genus-three hyperelliptic/non-hyperelliptic Jacobian isogeny mechanism | Retain as source for the open polarization-profile hypothesis, not as a new curve population |

The owner’s previous database includes anomalous trace-one curves,
smooth selected groups, pairing transfer, supersingular curves, binary
GHS descent, generalized coverings, quadratic and cubic torus
branches, Type I/II configurations, proper-subfield points, isogenous
representatives, CM conductor leads, endomorphism/torsion leads and
twist-input conditions. These entries overlap and are not disjoint
curve populations. Published cases absent from the model inventory
are a concrete import queue, with local status exposed in the
mechanism table.

The exact novelty question for the seven open hypotheses remains
unresolved against a full theorem-by-theorem literature comparison.
The previous bounded audit is
[here](../weak_curve_conjectures_20261012/LITERATURE.md).
Before claiming a newly discovered family, compare equations,
covering group, isogeny condition, field of definition, source genus,
twist and parameter equivalence with the complete published tables.

An external “largest database” comparison also needs a common unit.
Count distinct exact field/model IDs, distinct isogeny classes and
family predicates separately; exclude alternate representations of
one model from the first count only after an explicit isomorphism.
The current database is a large, auditable local integration. Its
coverage statement names the fixed sources and their domains.
