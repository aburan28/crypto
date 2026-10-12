# Lean statements and open proof obligations

| Statement | Lean declaration/status | Exact scope |
| --- | --- | --- |
| Two-branch inclusion-exclusion | `CensusCounting.inclusion_exclusion` | Arbitrary finite lists of Boolean branch flags; duplicated rows count as duplicated entries |
| Union plus zeros equals population | `CensusCounting.union_zero_partition` | Same lists and predicate meanings |
| Union count is bounded by population | `CensusCounting.union_bounded_by_population` | Consequence of the checked partition |
| Truth/completeness of imported family labels | Parent native/PARI certificates | Not established by these counting lemmas; native review separately checks CSV consistency and trace-domain completeness |
| `WCQ1:gap-count` | Open | Actual finite torus/CM support sets, exceptional CM units, and an asymptotic bound with quantified constants |
| `WCQ1:cm-mass` | Statistical premise open | Formalize selected order sets and occurrence law; proving count arithmetic does not validate the occupancy distribution |
| `WCQ1:nonhyper-only` | Open | An actual Type I/II curve, smooth non-hyperelliptic cover and Weil-restriction isogeny, in one exact zero class |
| `WCQ1:two-towers` | Open | Both actual coverings and their distinct base fields; ordinary isogeny-class descent exclusion, beyond j/model degree |
| `WCQ1:smaller-genus` | Open | Quotient curve over the base, actual Jacobian map onto the full restriction, dimension bound, primitive source conditions |
| Elementary scaling rejection in the smaller-genus screen | Readable algebra; Lean field lemma pending | `alpha^q=c*alpha`, c in the base, alpha nonzero gives `lambda=c` and base-defined Legendre j; this is not a proof of general quotient rigidity |
| `WCQ1:primitive-j` | Open | Actual finite-field Frobenius periods, complete anharmonic j fibre, primitive torus membership and exceptional invariants |
| `WCQ1:quartic-profile` | Open | Rational isotropic kernels, actual principal polarizations and certified geometric target types |

The compiler is pinned to Lean 4.34.1 with no external package dependency for
the new counting lemmas. Build the actual module, audit transitive axioms and
run the bundled `leanchecker --fresh --verbose CensusCounting`.

The earlier [Lean proof PR #1710](https://github.com/aburan28/crypto/pull/1710),
commit `038cad83dc146bbcf88504e82dec12b9a7fdb6f0`, supplies a separately
versioned finite-field twist/descent package. It is not assumed present on
this branch or merged. Its ordinary-curve CM/endomorphism bridges remain
open as documented there. Reuse those field statements explicitly when the
geometric conjectures become proved lemmas. Do not insert conjectures as
axioms, assume the desired degree/genus conclusion, or relabel counting
identities as proofs of weak-curve geometry.
