# Published odd-characteristic elliptic covering cases

This companion index makes the 36 numbered elliptic source-curve cases in
[Iijima–Momose–Chao, 2015, Section 5](https://eprint.iacr.org/2015/805.pdf)
queryable by extension degree `d`, binary covering rank `n`, excess genus
`e = g(C)-d`, covering genus, and allowed square class of the coefficient.
The paper's completeness domain is odd characteristic, prime extension degree
`d in {2,3,5,7}`, and a `(2,...,2)` covering over the base field. A case is a
parameter family, not one exact elliptic model or one isogeny class. The
index is separate from the exact-model database in
[`../weak_curve_database_20261012/`](../weak_curve_database_20261012/)
while that directory holds another Conductor write lease.

The five minimum-genus cases also have a queryable factor-pattern skeleton
in [`minimum_case_patterns.csv`](minimum_case_patterns.csv) and the
`minimum_case_recipe` view. These are parameter-family recipes; the paper's
full root alternatives and covering conditions still govern each realization.

The imported table is **case metadata**. Individual defining equations,
parameter equivalence, finite-field existence, model isomorphism, covering
maps, and subgroup transport are not yet normalized into this index. Its
`eta` value preserves the paper's convention: either 1 or a nonsquare in
`k_d`; `one` permits the square coefficient only. The values are transcribed
from the paper's printed pages 17–18 and retained in [`cases.csv`](cases.csv).

Build and query:

```sh
sh build.sh
sqlite3 -header -column published_cases.sqlite \
  'SELECT case_no,d,n,cover_genus,permitted_square_class FROM case_profile WHERE e=0 ORDER BY d,case_no;'
sqlite3 -header -column published_cases.sqlite \
  'SELECT * FROM minimum_case_recipe ORDER BY case_no;'
rustc --edition 2021 --test ramification_orbits.rs -o /tmp/ramification_test
/tmp/ramification_test
rustc --edition 2021 ramification_orbits.rs -o /tmp/ramification_orbits
/tmp/ramification_orbits > ramification_orbits.csv
cd lean && lake build PublishedCases && lake env lean AxiomAudit.lean
lake env leanchecker --fresh --verbose PublishedCases
```

The native screen enumerates the paper's binary ramification-mask classes:
for a listed polynomial `F | (X^d+1)` over `F_2`, set `Fhat=(X^d+1)/F`,
enumerate `a` coprime to `F`, and quotient masks `a Fhat` by cyclic rotation.
It does not decide whether a finite-field curve realizes a mask. The
[`report`](REPORT.md), [count chart](case_profile.svg), and
[construction diagram](construction_pipeline.svg) specify the remaining gates.

The [Lean proof map](PROOF_MAP.md) identifies the checked polynomial
factorizations and the structural covering statements still outside the
formalization. The [`VALIDATION.md`](VALIDATION.md) file records exact checks.
