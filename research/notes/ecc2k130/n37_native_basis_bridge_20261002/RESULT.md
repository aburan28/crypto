# Native n37 field-basis bridge and degree-73 pullback preflight

Status: **transport prerequisite passed; PDP yield and full ECDLP cost not measured**.

The frozen n37/L1024 control uses the registered source model
`icv1-f2m37-tm534059-32aad96b`, whose field modulus has low terms
`{0,1,4,6}`. The archived degree-73 descent uses a different polynomial
basis with low terms `{0,1,2,3,4,5}`. Reinterpreting a coordinate word under
the other modulus changes the point. This round tests whether the actual
frozen public workload can cross that representation boundary, reach the
archived descendant with its sign intact, and return through a rational
point pullback. It does **not** use a fixture logarithm to rebuild a target.

The inputs are the archived
[`raw.json.gz`](../../../koblitz_isogeny_descent_37_results_20260925/raw.json.gz)
(SHA-256 `eb4773d556886672e8b80735fc486c16b73765fd5069ff81133100bfc3eafa90`),
the n37/L1024 campaign
[`FROZEN.json`](../disjoint_cold_v2_20261001/FROZEN.json)
(SHA-256 `da958a3f1117dd1b88703fed2055c5c9f64255516e1e920f35d193bf53c6a88d`),
and block-00
[`points.jsonl`](../disjoint_cold_v2_20261001/fixtures/n37_L1024_b00.points.jsonl)
(SHA-256 `187ec04fe50326bbb2f17dadf76056841f04af37a8abb8ff7f502fcd531711ad`).
The point file contains 1,024 public coordinate pairs and no planted
logarithms. The separate historical fixture file has published scalars, but
the new transport tests do not read it.

`BinaryFieldBasis::discover` finds all 37 roots of the control modulus in
the archived field, chooses the least coordinate word, and checks that its
powers form an invertible binary matrix. The selected image of the control
polynomial generator is **10156182909**. `from_generator_image` reconstructs
that frozen map without rediscovery; both directions reject out-of-field
words. Discovery and construction are reusable setup and must be charged in
any cold-cost comparison.

The acceptance gate was field arithmetic preservation, a two-way point
round trip, curve and subgroup membership on both representations, signed
homomorphism through the degree-73 map, and a unique forward-checked
pullback. The native Rust tests now establish:

| Check | Observed result |
|:--|:--|
| Two distinct degree-8 field polynomials | All 256 words round-trip and preserve squaring; sampled products commute with the map. Identity and invalid-root cases are checked. |
| Frozen n37/L1024 block 00 | 256 deterministic field pairs preserve squaring, multiplication, and the basis-map round trip. All 1,024 public points map to the archived source model, map back exactly, lie on the archived curve, and retain subgroup order 230603167. The first 16 also preserve fixed-scalar multiplication and sampled additions. |
| Composed control-basis → degree-73 leaf map | All 1,024 public targets land on the certified descendant with correct subgroup order and sign; the first 16 satisfy `φ(G+Q)=φ(G)+φ(Q)` as complete points. No published fixture scalar is accessed. |
| Rational point inverse | The frozen generator, four held-out targets, the 2-torsion point, infinity, and 42 distinct sign-class points selected natively on the descendant pull back exactly. The native pullbacks also map into valid control-basis subgroup points. Every nontrivial preimage is accepted only after a full forward-map equality check. Non-coprime rational kernels and off-curve images are rejected. |

The inverse solves the degree-73 abscissa equation
`(X+x)h(x)^2 + x h(x)h'(x) + x^2 h'(x)^2 = 0` in the source field, lifts both
signs, and accepts the unique one whose complete image equals the requested
point. It applies when the isogeny degree is coprime to the rational group
order; here `gcd(73,137439487532)=1`. Root finding is a real setup cost for
each pulled-back base point, not a free correspondence.

The 42 native candidates are selected without a target by scanning descendant
abscissae in increasing order, taking the first rational lift, multiplying by
the cofactor 596, and retaining the first 42 distinct image abscissae. This
pins a reproducible candidate sign-class set; it is not yet an admitted PDP
factor base because no relation yield or rank has been measured.

Its `K=42` matches the control's *log-column count* only. A control column
can represent a signed Frobenius orbit, while an unaugmented descendant
column represents one signed pair. Matching column count alone does not
match useful relation support or action cost; both must be explicit in the
four-policy comparison.

Their signed descendant, archived-source pullback, and control-basis pullback
coordinates are frozen in [`NATIVE42.json`](NATIVE42.json), SHA-256
`bd4bd8af982bcc65234ae1fae1ad52eb0ed41a1c203bf4d2556fd8e87f788a4c`.
The Rust test recomputes and compares the entire manifest. To regenerate an
independent copy, set `N37_NATIVE_BASE_MANIFEST_OUT` to a new file outside
the repository when running the filtered composed-transport test.

Replay with `cargo test --lib cryptanalysis::binary_field_basis::tests` and
`cargo test --lib cryptanalysis::binary_velu::tests`. The previous signed-map
preflight is in [`n37_native_velu_transport_20261002`](../n37_native_velu_transport_20261002/RESULT.md).
No new relation count, rank, solver completion, online charge, or rho
comparison was measured here. The archive's endomorphism-order discriminant
change from `-7` to `-37303 = -7·73²` is structural evidence, not evidence
that a descendant-native PDP is easier.

The next scored round must freeze four target-blind bases and report both
independent log-column count and the **actual signed support under the
available, costed automorphism action**: original source, descendant-native,
source-to-descendant transported, and descendant-to-source pullback. It must
include the basis bridge and inverse costs, state any transported
automorphism action and its charge, independently check signed relations,
and measure rank and held-out-target completion. The current descendant-native
`K=42` base has an arity-five counting support ceiling of 16.10084%, so the
first arity whose counting capacity exceeds one is six; that threshold does
not predict actual yield. Only a cold full-rank logarithm recovery against
the matched strong rho control can support a speed claim. A failed yield or
rank gate remains a reported result.
