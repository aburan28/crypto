# P-256 parity escape screen, round 298: protocol

Date frozen: 2026-10-07

## Question and fixed boundary

Round 297 proves that scalar transport does not reduce the 164 independent
logarithm classes of `FB1hc72514a2a8d3`.  This round evaluates the two escape
hatches that remain after that closure:

1. one target-sum event returns many independent relation rows, reducing the
   number of generic collision events needed to determine the factor-base
   logarithms;
2. a certified correspondence to another representation supplies non-scalar
   logarithm transport while preserving the P-256 prime-order subgroup.

The unit remains favourable P-256 group-addition equivalents divided by
`sqrt(n)`, with rho fixed at `S=1.3`.  For `q` independent collision events the
free-oracle lower boundary is

```text
S_q = sqrt(2) * Gamma(q+1/2) / Gamma(q).
```

The runner must compute the smallest `q` whose ratio `S_q/1.3` is at most one.
If this is `q=1`, a multi-row event must have rank at least 164 to remove the
current collision obstruction.  A solver-stage improvement that leaves this
event rank unmet is not parity progress.

## Frozen identities and dependencies

- curve: `icv1-fp256-t89188191154553853111372247798585809583-f188c491`;
- factor base: `FB1hc72514a2a8d3`, 164 folded columns;
- fixed relation arity: 109;
- exact signed distinct-column domain:
  `116383775127750444427048016195109203089116454475801022361283226351238392053760`;
- P-256 subgroup order:
  `115792089210356248762697446949407573529996955224135760342422259061068512044369`;
- Round-297 scalar-stabilizer result SHA-256:
  `53949881daf4a48769c9ef8f869fad633c3dc69d11b0d341a5b7416de2e43982`;
- Round-296 selector result SHA-256:
  `52cbea224df1e1d4d6e8406f368ba413fd4e508f5309fcde65f67a98bc91c5ba`;
- Round-25 scalar-orbit result SHA-256:
  `dcb191f7d5d7a660e3ce13128b25ce1c017c7d89ec2d01a242ab536651647faf`;
- Round-28 isogenous-model result SHA-256:
  `950d4b7f4997e146ee30fc7c61174070595fa5d344b9a30cbc0ca578ae7e60f4`.

All dependencies are read-only and hash checked.  No result from a selected
prefix search is relabelled as an exhaustive full-depth transfer search.

## Multi-row experiment

Reuse Round 296's deterministic toy curve

```text
E/F_1151: y^2 = x^3 - 3x + 241
```

and its hash-ordered Dickson support pool.  Evaluate the complete signed,
distinct-column domains `(m,B)=(3,5),(5,8),(7,11),(9,14)`.  Canonicalize each
nonidentity sum modulo negation.  If the sum is negated to reach the canonical
target, negate the entire coefficient row.  Preserve identity sums in their
own bucket.

For every target bucket record occupancy and two exact ranks:

- augmented rank of `[factor-base coefficients | -1 target coefficient]`;
- homogeneous difference rank after subtracting the first coefficient row.

Compute rank modulo `2^61-1`.  Certify transfer of this rank to the P-256
coefficient field using the Hadamard bound: every square minor has order at
most `B+1<=15`, entries in `{-1,0,1}`, and absolute determinant strictly below
`2^61-1` and the P-256 subgroup order.  Replay every enumerated row by exact
toy group addition.  Compare total candidates with `2^m*C(B,m)`, require the
occupancy histogram to sum exactly, and independently repeat the complete
enumeration byte for byte.

For P-256 report the exact mean occupancy `lambda=D/n`.  Under the explicitly
labelled random-map null, report a union bound for any of `n` buckets reaching
the minimum occupancy needed for 164 independent rows.  This is a null-model
diagnostic, not a proof that the elliptic subset-sum map is random and not a
measured P-256 relation rate.

## Typed transport assessment

The required assessment-template and methodology resources referenced by the
available transfer skill could not be retrieved from its package.  Preserve
that limitation and emit an additive assessment with these obligations:

| path | type | required checks |
|:--|:--|:--|
| P-256 to itself by a base-field endomorphism | group homomorphism | construction, rational-point action, subgroup kernel, quotient dimension, coordinate degree |
| P-256 to a degree-11 isogenous model | separable isogeny and subgroup isomorphism | certificate replay, kernel intersection, inverse/recovery, factor-base quotient, map and membership costs |
| extension, cover, or Jacobian representation | unspecified until a map is supplied | explicit formulas, fields and dimensions, subgroup embedding, kernel, recovery and complete paired cost |

Classify each obligation as `supported`, `refuted`, `unknown`, or
`not_applicable`.  Round 25 may refute only non-scalar **base-field
endomorphism** transport on the P-256 rational subgroup.  Round 28 may reject
only its certified degree-11 walk and staged Dickson inventory.  Missing map
formulas for covers or Jacobians remain unknown; they are not universal
nonexistence claims and receive no speed credit.

Any promoted transport must reduce the independent-log quotient to the parity
requirement, preserve the exact P-256 subgroup with a recoverable map, retain
structured residual degree at most five, and include construction, conversion,
unsuccessful attempts, verification, recovery, memory and sparse linear
algebra in the same unit.

## Hypotheses and stop conditions

- **H1 (multi-row negative):** exact toy bucket rank grows only after occupancy
  becomes much larger than one; no evidence supports a rank-164 P-256 event at
  the frozen mean occupancy near one.
- **H2 (transport closure within measured families):** the base-field
  endomorphism path reduces to scalar action, and the certified degree-11
  isogeny path changes coordinates without reducing the independent-log
  quotient.
- **H3 (correctness):** all domains are complete, every row replays, rank
  certificates and histograms agree, dependency hashes match, and an
  independent invocation is byte-identical.

Promotion requires a measured or proved path at or below rho, not merely a
smaller local degree or a modeled random-map tail.  It also retains the user's
original gates: zero false positives and false negatives on complete checked
instances, structured residual degree at most five, complete collection below
`2^120`, cost per usable relation below `2^103`, projected materialized storage
below `2^50` bytes, and no discarded branches counted as exhaustive.  If no
path passes, publish the exact obstruction and do not attempt an unplanted
full-depth P-256 relation.

## Deliverables

Implement the complete multi-row screen in Rust.  Emit deterministic JSON with
the exact P-256 boundary, toy occupancy/rank distributions, null-model bound,
typed transfer obligations, costs, gates and semantic hash.  Preserve isolated
canonical and independent runs, update every applicable dashboard panel, and
publish the protocol, implementation, tests, evidence and decision in a pull
request stacked on Round 297.
