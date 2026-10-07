# P-256 low-width Dickson-union screen, round 295: protocol

Date frozen: 2026-10-07

## Question and boundary

Round 294 proves that an independent-log factor base cannot reach the P-256
coefficient-domain threshold with fewer than 164 columns in one fixed arity.
The smallest complete Dickson fibre actually built there has 266 columns and
an exact favourable generic-collection boundary of 17.734068 times rho.  This
round asks whether a deterministic union of complete low-depth Dickson fibres
can realize the theoretical width boundary near `B=164` without replacing the
quadratic Dickson-chain membership representation by an arbitrary point list.

For `K=B` independent logarithm classes the exact negation-folded collision
boundary remains

```text
T_K = sqrt(2*n) * Gamma(K+1/2) / Gamma(K).
```

All ratios use the frozen rho reference `S=1.3`.  They grant a perfect free
decomposition oracle, success probability one, one operation per collision
sample, and zero setup, verification, matrix, recovery, or representation
cost.  A lower ratio is a factor-base boundary advance only; it is not an
end-to-end speedup.

## Frozen identity and dependency

- curve: `icv1-fp256-t89188191154553853111372247798585809583-f188c491`;
- curve model and point encoding: the existing native
  `p256_dickson_factor_base` implementation;
- Round-294 dependency:
  `research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_variable_arity_round294_20261007/variable-arity-result.json`;
- required Round-294 SHA-256:
  `9716fc1ba7d24a023627e8701094ca870b2121a84f7a0148784f70565640815a`;
- comparison base: depth-9 `FB1h91cc12460e5a`, 266 columns, first
  covering arity 59, exact boundary 17.73406835956266 times rho;
- theoretical comparison: 164 columns, covering arities 109 and 110,
  exact boundary 13.920747397073491 times rho.

## Deterministic candidate construction

Build 256 complete depth-8 fibres and 256 complete depth-6 fibres.  Candidate
index `i` at depth `d` derives its root exponent from

```text
SHA256(curve || "/dickson-union-round295/depth-" || d || "/" || i)
```

reduced modulo `2^(96-d)` and forced odd.  Reject only a builder-declared
degenerate terminal; preserve its index and reason.  Rebuild every accepted
component and compare every emitted point row.

Form every ordered `(depth-8 index, depth-6 index)` pair.  A union column is a
unique negation-folded affine abscissa; duplicate signed points are retained
once and reported.  Component order is canonical by `(depth,index,FB1)`.  The
union uses the repository `ecbench.factor_base/v1` preimage and `FB1` hashing,
with family `dickson-torus-union`, params carrying the ordered component
identities, depths, root exponents, terminals, and point-key encoding.  The
compact artifact records its full SHA-256 and point-set SHA-256; the large
point dump is reproducible derived data and need not be committed.

## Ranking and exact arithmetic

For every distinct union width `B`, enumerate every fixed arity `2..B` with

```text
D(B,m) = 2^m * C(B,m)
```

using exact integers and verify `sum D(B,m)=3^B`.  A pair is covering when at
least one fixed arity has `D(B,m) >= n`.  Rank covering pairs by:

1. smallest `B`, equivalently the smallest exact collision boundary;
2. smallest duplicate count;
3. smallest depth-8 index;
4. smallest depth-6 index;
5. lexicographically smallest union FB1.

Independently rebuild the winner and require a byte-identical compact
identity, point digest, width, component inventory, and boundary row.

## Degree and scope

Each component remains a complete Dickson fibre whose chain equations have
generator degree two.  The union is a disjunction of two such components;
the report must keep the component-local generator degree separate from the
solving degree of a global 109/110-term union-membership system.  The latter
is unset unless directly measured.  In particular, two quadratic component
descriptions do not establish structured residual degree at most five.

No relation search, full-depth unplanted attempt, or end-to-end projection is
authorized unless all promotion gates pass.  A candidate cannot pass merely
because its factor-base dump fits in memory.

## Frozen hypotheses and gates

- **H1 (realization):** at least one exact two-fibre union has `B=164`; if
  none does, report the minimum covering width without expanding the frozen
  search.
- **H2 (factor-base boundary advance):** the selected exact boundary is below
  the depth-9 control's 17.73406835956266 times rho.
- **H3 (theoretical closure):** no valid covering union can have `B<164`, and
  any `B=164` winner exactly matches Round 294's 13.920747397073491 boundary.
- **H4 (identity):** all component and union point rows verify, the independent
  winner replay is byte-identical, and false identity/deduplication matches
  are zero.

Promotion additionally requires all of the user's original gates: structured
residual degree at most five, complete relation collection below `2^120`, cost
per usable relation below `2^103`, projected materialized storage below
`2^50` bytes, zero false positives and false negatives on complete checked
instances, and no discarded probabilistic branches counted as exhaustive.
Unset is failure.  If any gate fails, report a reproducible negative result
and do not attempt an unplanted P-256 relation.

## Deliverables

Implement and test the screen in Rust.  Emit deterministic JSON with
dependency, boundary-table, component-census, winner, replay, and semantic
hashes.  Preserve an isolation receipt and compact result report.  Update the
canonical progress data, scoreboard boundary/exponent facts, a dedicated
panel, and the generated panel index.  Publish the protocol, implementation,
tests, evidence, and decision in a pull request stacked on Round 294.
