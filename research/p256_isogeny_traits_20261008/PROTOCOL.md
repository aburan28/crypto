# P-256 certified isogeny trait census

Status: **preregistered; implementation and execution pending**

Date frozen: 2026-10-08

Result class: **complete under a versioned structural-trait schema, not an
exhaustive set of mathematical properties and not an ECDLP speedup**

## Question and acceptance boundary

For every curve in the already certified 1,064,000-curve P-256 union, derive
the same explicit set of class, model, path, edge, conductor, and direction
traits. Preserve exact values where the evidence determines them and typed
unknowns everywhere else. In particular, determine whether the degree-11 and
degree-13 edges are horizontal, ascending, or descending in the isogeny-volcano
sense. The grid's drawn row and spine axes are construction coordinates and
are not directional evidence.

The accepted boundary is a native streaming generator and an independent
native replay. They pass only when both bind the exact source certificates,
account for 1,064,000 distinct input and output curves, recompute every field,
agree on a chained output digest, and account for every incoming edge in
exactly one direction-status bucket. A missing calculation is a failure unless
the frozen schema explicitly permits a `null` value with a nonempty status and
reason.

## Frozen inputs

| input | curves | gzip SHA-256 | content SHA-256 | record-chain SHA-256 |
|:--|--:|:--|:--|:--|
| million grid | 1,000,000 | `87522b6de08e64804dcee17bba058eb91922d6b8dd9853cd32523f67a3d0cbdc` | `2a1a27d06f48b22b85ce3c96ac4017f5af6a0256e9d1ad277e3fb1f516475454` | `4de6e93e0d7a90ff53a26f8f1fd1cbd03ab55f701a25f8e49874cb8bead40302` |
| rows 1000 through 1063 | 64,000 | `90f22f7aadaef78cff5a72889e5539b2ae81459d57d7a84b37a9009e96a3f701` | `a24593f9d1d1c880b63f1da86d10e375981f72aaa358e1eae428e63b3d4714b1` | `7d542b93aec2f21d9237e86f9d079cee4f858b5650bc0995c2f85ca0857cf87d` |

The prior exact union receipt reports 1,064,000 input curves and 1,064,000
unique full-width j-invariants. This round consumes those immutable bytes. It
does not regenerate or adaptively replace a curve.

## Frozen schema

The artifact is deterministic gzip JSON Lines with four record kinds.

### Class record, emitted once

- exact field prime `p`, group order `N`, trace `t`, cofactor, and twist order;
- Hasse and ordinary-curve checks, group structure implied by the proved-prime
  order, and rational torsion consequences;
- Frobenius polynomial `X^2 - tX + p`, signed discriminant
  `Delta_pi = t^2 - 4p`, absolute bit length, and factorization status;
- fundamental CM-field discriminant `D_K`, Frobenius-order conductor
  `f_pi` in `Delta_pi = D_K f_pi^2`, and the exact equality check when a
  complete decomposition is available;
- factor records for `N`, the twist order, `abs(Delta_pi)`, and `f_pi`. Every
  record includes the factors found, their product, remaining cofactor,
  primality-test method, and one of `complete_proved`,
  `complete_probable_prime`, or `partial_composite_cofactor`;
- extension-field orders for degrees 2, 3, 4, and 6, plus the bounded embedding
  search already used by the native walk;
- for degrees 11 and 13: `Delta_pi mod ell`, splitting kind, exact
  `v_ell(Delta_pi)`, volcano depth `floor(v_ell(Delta_pi)/2)`, whether
  `ell` divides `f_pi`, expected rational neighbors when provable, and the
  directional theorem used.

The exact factorization budget is trial division through `2^20`, followed by
the existing deterministic-seed Pollard-Brent path for at most 64 seeds and
10,000,000 iterations per unresolved composite. The result keeps unfinished
cofactors. A fixed-base Miller-Rabin pass is labeled probable prime, never
silently promoted to a proof.

### Curve record, emitted once per curve

- input number, global `(x,y)` coordinate, full ICV1 identity, EC1 alias, full
  curve UID, and source-certificate digest;
- exact canonical `a`, `b`, and `j`, plus recomputed j consistency;
- `4a^3 + 27b^2`, standard equation discriminant
  `-16(4a^3 + 27b^2)`, `c4 = -48a`, and `c6 = -864b`, all as canonical field
  elements;
- zero/nonzero and quadratic-character values for `a`, `b`, `j`, the model
  denominator, `c4`, and `c6`;
- `j = 0`, `j = 1728`, algebraic-closure automorphism-order case, the
  canonical-model branch, `a = -3` isomorphism availability, prefix-residue
  count, and signed coefficient bit lengths;
- recorded generator validity, nonsingularity, and order status;
- symbolic path degree `11^x * 13^y`, its exact bit length, and a decimal-log
  interval. The enormous expanded integer is not duplicated per record because
  the exponent pair is an exact lossless representation.

Integer factorization of the canonical representatives of `a`, `b`, `j`, or
the equation discriminant is excluded. It changes under field-equivalent
encodings and is not an elliptic-curve invariant. Their exact field values and
quadratic characters are included instead.

### Edge record, embedded in every non-root curve record

- source and target coordinates and identities, kernel-certified prime degree,
  separability and cyclicity status, predecessor j-invariant, and construction
  axis (`row_degree_11` or `spine_degree_13`);
- splitting and conductor facts inherited from the class record;
- volcano direction (`horizontal`, `ascending`, `descending`, or `unknown`),
  evidence status, and reason;
- source and target endomorphism-ring conductors and their ratio when proved,
  otherwise explicit nulls with reasons.

An odd-prime degree is horizontal when the class evidence proves volcano depth
zero and the Frobenius polynomial splits at that prime. Coordinate orientation
alone never assigns a volcano direction.

### Summary record

Counts cover input schemas, curves, edges, construction axes, every direction
status, special-j cases, model branches, quadratic characters, existing screen
values, unresolved fields, discrepancies, and selected extrema. Large-prime
and conductor findings name the class value and evidence status rather than
being copied into one million rows.

## Independent replay and negative controls

The verifier opens both source certificates again. It validates stored and
content hashes, certificate schemas, header compatibility, summaries, record
chains, source ordering, global coordinates, exact non-overlap, and all curve
identities. It then recomputes every census record without trusting a recorded
derived field. It rejects:

1. a changed source byte or reordered source list;
2. an omitted, extra, or altered trait field;
3. a direction inconsistent with its conductor evidence;
4. a factorization whose product does not recover the declared integer;
5. a duplicated coordinate, j-invariant, ICV1 identity, EC1 alias, or UID;
6. a count or chained digest that differs from the replay; and
7. an unknown value without a reason.

Focused fixtures must include adjacent certificates that pass, overlapping
certificates that fail, a tampered source binding, a tampered trait value, a
direction downgrade/upgrade error, and a minimal valid census.

## Stop rules and evidence classification

Stop and preserve the attempt if free temporary storage falls below 512 MiB,
source validation fails, any duplicate appears, the trait artifact exceeds
512 MiB stored, or generation or replay exceeds 90 minutes. Record elapsed and
CPU time, peak RSS, compiler and source revision, input/output byte counts and
hashes, record counts, and every unresolved factorization.

The census may select a follow-up hypothesis. It cannot establish an index-
calculus or ECDLP speedup. Relation generation, linear algebra, individual-log
descent, target transport, and matched Pollard-rho cost remain unknown until a
separate preregistered one-target solver experiment measures them end to end.

## Reporting requirements

The report must contain an editable diagram distinguishing construction axes
from volcano directions, a quantitative graph of direction/model/status
counts, primary-source links for ordinary isogeny classes and volcano theory,
the complete requirements table, negative controls, resource costs, and a PDF
rendered and visually inspected page by page. No canonical speed scoreboard is
updated unless a separately qualified solver result exists.
