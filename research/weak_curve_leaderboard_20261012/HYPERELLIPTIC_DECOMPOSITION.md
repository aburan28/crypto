# Hyperelliptic decomposition: divisor relations and Jacobian factors

A hyperelliptic model `H/F_q` is an elliptic-DLP attack *destination* only
after the exact source group and map into `Jac(H)(F_q)` are fixed. A DLP posed
directly in `Jac(H)(F_q)` starts there. The [cover checker](../../docs/curves/COVERS.md)
already certifies same-field maps `H → E` for supported elliptic models. The
[transfer registry](../../docs/curves/JACOBIAN_TRANSFERS.md) prices several
Jacobian solvers by formula. The missing per-model records are divisor
relations, explicit Jacobian factor maps, subgroup images and complete costs.
The [WAM1 attack map](attack_map.json) gives these separate research checks.

## Per-curve detection properties

The machine-readable [property index](decomposition_properties.json) now binds
eleven hyperelliptic checks to every one of the 922 exact atlas identities.
Run `python3 build_decomposition_properties.py --check` to replay its joins,
and `python3 validate_decomposition_properties.py` for an independent
cross-source audit. The site exposes the same states through its
**Hyperelliptic property** and **Property state** filters; opening a curve
shows its complete check list and cover polynomial/map coefficients.

| Property | Recorded evidence | Result in this snapshot |
| --- | --- | --- |
| Same-field hyperelliptic cover | Exact ICV1 identity, full model hash, checker certificate or formula replay, field, genus, degree and polynomial map | 917 source-checker certificates; five formula replays |
| Odd-field cubic and quadratic genus-three branches | Exact model label and necessary screen, retaining exclusions and unclassified survivors | 14 cubic and ten quadratic positive labels; 247 large cubic screen survivors remain unclassified |
| Binary GHS lower-field route | Exact magic number when available, base/relative degrees, projected genus and transport-size screen | 49 curves with exact-invariant route records; cover construction and subgroup map open |
| Divisor factor-base and large-prime relations | Frozen base, natural-query yield, partial-relation graph, rank and target descent fields | Measurement fields present and null for every curve |
| High-genus smoothness | Size-admissible exact-invariant genus projection and smoothness-measurement field | 37 genus projections; instance smoothness open |
| Jacobian factors and elliptic quotient | Destination-field Weil polynomial, factorization, explicit projection and complement dimension | Cover geometry supplies the dimension screen; factorization and projection open |
| Subgroup transport and route cost | Selected order, numerical degree gcd where known, generator/target images, kernel, formula class and measured total | Numerical gcd is a screen; images, kernel and complete measured cost open |

The 917 certificates are imported from
[`docs/curves/covers.json`](../../docs/curves/covers.json) by exact slug,
ICV1 identity and full model SHA-256. That checker has already replayed the
polynomial identity, nonsingularity, genus and degree. The three small
prime-field examples and two pairing records are outside that source
registry. For these five, the builder instantiates the documented
`H: v²=F(u²+c)` construction, checks the resulting polynomial identity and
squarefreeness over the source prime field, and marks the row
`formula_replayed`; source primality and curve parameters remain the declared
inputs. The source cover report and transfer cost report were built at
different registry revisions, so the join checks each imported model hash
rather than asserting equality of the whole registry digests.

An exact same-field cover is common in this catalog and does not identify a
smaller-field route. The `positive_model_label` state refers only to its
specified cubic or quadratic branch. `necessary_screen_passed` does not
assert a cover, while `excluded_by_screen` excludes that recorded branch,
not every possible hyperelliptic construction. Formula costs and projected
genera stay distinct from constructed lower-field covers and measured divisor
relations. The [proof map](PROOF_MAP.md) continues to list the Lean bridge
from cover geometry to Jacobian maps as open; this property index adds no
new structural theorem.

## Decompose divisors for index calculus

For a hyperelliptic Jacobian, represent a reduced divisor in the declared
field and curve model, for example by Mumford data `(u,v)` where applicable.
Choose a factor base of prime divisors of bounded degree. A successful
decomposition expresses a sampled Jacobian element as a sum of factor-base
classes; verified relations then feed a sparse linear system modulo the
selected subgroup order. [Sarkar–Singh](https://eprint.iacr.org/2014/815)
give a decomposition method with sieving for small genus. The
[double-large-prime method](https://eprint.iacr.org/2004/153) admits partial
relations and recombines them through a collision graph. For larger genus,
[Enge](https://doi.org/10.1090/S0025-5718-01-01363-1) makes the
smoothness bound and factor-base choice central to subexponential algorithms.

The reproducible row needs `q`, genus, exact `H`, divisor encoding, selected
order `r`, factor-base construction and digest, degree bound, query law,
attempt count, full/partial/failed decomposition counts, large-prime graph,
independent relation rank, matrix size, final linear algebra, target descent,
scalar replay, and exclusive phase costs. A smoothness probability estimated
from a cost formula is a prediction; natural-query measurements provide the
instance evidence. If the original problem is on `E`, charge cover search,
map construction and subgroup transport as well.

## Decompose the Jacobian up to isogeny

A different question asks whether `Jac(H)` has lower-dimensional factors over
the **declared field**. Compute and factor its Frobenius/Weil polynomial as a
screen, then construct an actual isogeny and projection maps. The abstract
factorization and a usable map are separate evidence states;
[Ahmadi–McGuire–Rojas-León](https://arxiv.org/abs/1409.8510) study
Jacobian decomposition through L-polynomials over finite fields. For a
genus-two curve with an elliptic cover, the Jacobian is isogenous to a product
of two elliptic curves under the relevant maximal-cover hypotheses;
[the explicit genus-two study](https://arxiv.org/abs/1209.3187) gives the
corresponding degree-`n²` product isogeny. A genus-three elliptic quotient
leaves a two-dimensional complementary factor, which may need a further
split test.

Field of definition is decisive. A lower-field curve `H/F_q` can cover an
elliptic curve `E/F_(q^d)` through a map defined only after extending to
`F_(q^d)`. That map need not split `Jac(H)` over `F_q`. Factor the Weil
polynomial and construct projections over the field where the destination
DLP is supposed to run; a factor appearing only after scalar extension is a
different route with different costs.

For an exact map `φ:H→E` of degree `m`, the pullback and norm satisfy the
standard relation `φ_*φ^*=[m]`. The [coverage audit](COVERAGE_AUDIT.md)
already uses this relation as the selected-subgroup injection gate when
`gcd(m,r)=1`; a non-coprime case needs a direct kernel-intersection check.
On the recorded same-field genus-two covers, an elliptic factor related to
the source `E` is expected. A split label alone therefore gives no cheaper
route to a selected discrete logarithm. The exact generator and target must
project to a factor where the stated solver and charged route improve the
chosen task.

An isogeny-decomposition record needs the factor polynomials and field,
factor equations, isogeny degree and kernel, maps and exceptional cases,
generator/target images, exact image order, cost of each map and factor DLP,
and an independent scalar-compatibility check. A factor of the polynomial
without an explicit map stays a structural candidate. The ordinary elliptic
isogeny graph in this snapshot contains only maps between elliptic curves;
higher-dimensional Jacobian factors require a separate typed transfer edge.

## Promotion on one catalog curve

1. Fix the source field, model, chosen order-`r` generator and target law.
   Construct `H`, its declared field, genus and cover map if the source is
   elliptic. Preserve the map certificate and its cost.
2. Compute the exact Jacobian order or Weil polynomial. Screen potential
   factors; construct and verify each proposed map before assigning a split
   finding. Check the selected subgroup's kernel intersection and factor
   projection.
3. For a divisor-index-calculus route, freeze the factor base and natural
   query distribution. Preserve zero-yield, timeout and partial-relation
   cells; measure useful rank, final linear algebra and individual descent.
4. Compare a previously unseen selected target against a matched source-group
   reference with the repository's [measurement contract](../../AGENTS.md).
   An exact cover, a polynomial factor, and a cost-model row remain their
   respective evidence grades until their missing gates pass.

This method sheet indexes published algorithms and the existing cover
certificates. The [proof map](PROOF_MAP.md) tracks formalization of the
structural transfer implications; a per-model WID1 finding is promoted only
through the explicit map, selected-subgroup and cost gates above.
