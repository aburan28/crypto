# P-256 Dickson-union scalar stabilizer, round 297: protocol

Date frozen: 2026-10-07

## Question and boundary

Round 295 constructs the 164-column Dickson-union factor base
`FB1hc72514a2a8d3`; Round 296 finds no structured residual degree at most five
for its global Riemann--Roch selector.  This round asks a narrower exact
question: does scalar multiplication globally permute the factor base and
therefore identify factor-base logarithms without a relation search?

Let `S` be the 328 nonidentity signed points represented by the 164 folded
columns.  If `[a]S=S`, multiplication by `a` partitions `S` into free orbits
of length `ord_n(a)`, because P-256 has prime subgroup order `n`.  Consequently
`ord_n(a)` divides 328.  The exact arithmetic precheck must verify that the
only divisors of 328 that also divide `n-1` are 1, 2, 4, and 8.  Thus every
possible scalar stabilizer lies among the eight eighth roots of unity in
`F_n^*`; testing those eight scalars is exhaustive, not sampled.

For an accepted signed stabilizer subgroup `H`, negation is already folded,
so the effective column action has order `|H|/2`.  Direct orbit enumeration
must agree with the quotient dimension `K`.  The favourable generic collision
boundary is

```text
T_K = sqrt(2*n) * Gamma(K+1/2) / Gamma(K),
ratio_to_rho = sqrt(2) * Gamma(K+1/2) / Gamma(K) / 1.3.
```

This grants a perfect free decomposition oracle, success probability one, one
group-addition equivalent per sample, and zero setup, representation, solver,
verification, matrix, or recovery cost.  It is a lower boundary, not an
end-to-end speed measurement.

## Frozen identity and dependencies

- curve: `icv1-fp256-t89188191154553853111372247798585809583-f188c491`;
- factor base: `FB1hc72514a2a8d3`, 164 folded columns and 328 signed points;
- factor-base artifact:
  `research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_dickson_union_round295_20261007/selected-factor-base.json`;
- required artifact SHA-256:
  `d27516ca40a612ecf3ebabfa8ae04776084c7948110e20f221da438e1e64d8f8`;
- Round-296 result:
  `research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_rr_selector_round296_20261007/rr-selector-result.json`;
- required Round-296 SHA-256:
  `52cbea224df1e1d4d6e8406f368ba413fd4e508f5309fcde65f67a98bc91c5ba`;
- rho reference: `S=1.3`;
- unquotiented exact boundary: `K=164`, ratio `13.920747397073491`.

The runner must reject a dependency hash, curve slug, factor-base ID, point
count, duplicate column, malformed sign pair, off-curve point, or identity
point before screening a scalar.

## Deterministic exhaustive screen

Find the least integer `h >= 2` such that

```text
r = h^((n-1)/8) mod n
```

has exact order eight.  Enumerate `r^j`, `j=0..7`, in exponent order.  Verify
that these are eight distinct eighth roots and that every possible order from
the divisibility theorem occurs in the enumerated set as dictated by the
cyclic subgroup.

For every scalar and every low-y column representative `P_i`, compute
`[a]P_i` with the native constant-time P-256 scalar multiplication path and
normalize exactly.  Look up the resulting abscissa among all 164 columns and,
when present, record destination column and whether the image is the low-y or
high-y signed representative.  A scalar stabilizes the signed factor base iff
all 164 images are present and the destinations form a permutation.  Replay
every accepted mapping as an exact group equality.  Hash each complete signed
action map in canonical source-column order.

Enumerate the induced accepted action on folded columns.  Require the accepted
scalars to be a subgroup, its action maps to compose exactly, its column-orbit
partition to cover every column once, and the directly measured orbit count to
equal the reported independent-log quotient dimension `K`.

Count all scalar multiplications, bit rounds, group additions, group doublings,
normalizations, lookups, image hits, replay operations, and failures.  Emit a
canonical JSON result and require an independent invocation to reproduce it
byte for byte.

## Frozen hypotheses, success, and stop conditions

- **H1 (expected negative):** the only signed scalar stabilizers are `+1` and
  `-1`; after the already-present negation fold, `K` remains 164.
- **H2 (exhaustiveness):** the divisor theorem and eight-root enumeration cover
  every possible scalar set stabilizer; no unevaluated scalar branch remains.
- **H3 (correctness):** every accepted action replays, every rejected root has
  at least one exact missing image, all action and subgroup checks pass, and
  the independent result is byte-identical.
- **Success target:** an accepted nontrivial folded scalar action reduces the
  exact collision boundary below 13.920747 times rho.  Parity additionally
  requires a boundary at or below 1.0 times rho.

Promotion still requires the user's original gates: structured residual degree
at most five, complete relation collection below `2^120`, cost per usable
relation below `2^103`, projected materialized storage below `2^50` bytes,
zero relation false positives and false negatives on complete checked
instances, and no discarded probabilistic branches counted as exhaustive.
Unset is failure.  If H1 holds or any promotion gate remains unset, publish the
negative result and do not attempt an unplanted P-256 relation.

## Deliverables

Implement and test the screen in Rust.  Preserve the protocol commit before
the implementation or run.  Emit the canonical JSON, an isolation receipt,
the exact theorem certificate, action-map hashes, operation counts, quotient
partition, boundary table, semantic evidence hash, and compact report.  Update
every applicable dashboard panel and generated panel index.  Publish the
protocol, implementation, tests, evidence, and decision in a pull request
stacked on Round 296.
