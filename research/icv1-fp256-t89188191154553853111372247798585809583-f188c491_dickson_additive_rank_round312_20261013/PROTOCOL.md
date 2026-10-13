# P-256 Dickson-hierarchy additive-rank screen, round 312: preregistration

Date preregistered: 2026-10-13

## Objective and fixed boundary

Continue the search for an index-calculus construction at Pollard-rho parity
on `icv1-fp256-t89188191154553853111372247798585809583-f188c491`, without
changing the registered factor base `FB1h2f8621cda105`.

Round 311 showed that half-turn atomization is invertible and leaves all
65,724 atom logarithms independent.  This round tests a different possibility:
the depth-18 Dickson hierarchy may force exact low-weight group identities
among the original 131,458 columns.  Such identities would be static relation
rows, reduce the independent-log rank, and preserve Round 19's structured
residual maximum four because the factor-base point set and S17 membership
system remain unchanged.

One unit is

```text
S = P-256 group-addition-equivalent operations / sqrt(n),
```

with Pollard rho fixed at `S=1.3`.  The exact negation-folded K-th collision
boundary is

```text
T_K / sqrt(n) = sqrt(2) * Gamma(K+1/2) / Gamma(K).
```

No stage or sampled result may be described as end-to-end parity.

## Preserved requirements

| requirement | evidence required in this round |
|:--|:--|
| exact P-256 factor base | native rebuild, FB1 and point-stream hashes |
| global log-rank reduction | replayed static relations and dual-modulus rank |
| zero false positives / negatives | complete toy reference and replay of every reported P-256 identity |
| structured residual degree at most five | unchanged-base proof plus hash-pinned Round-19 maximum four |
| collection below `2^120` | complete projection after static-rank credit, with setup and sparse algebra separated |
| usable relation below `2^103` | same unit and retained S17 success probability |
| peak materialized storage below `2^50` bytes | measured logical records, peak RSS, and projected attack storage |
| no discarded branch counted exhaustive | deep samples labelled diagnostics and excluded from rank/cost promotion |
| unplanted full-depth relation | attempted only if every promotion gate passes |

The original 138,031-row collection target and duplicate/rank allowance remain
visible.  Static identities may replace dynamic matrix rows only to the extent
certified by their actual rank; sampled identities receive no such credit.

## Frozen inputs

- Round 2 factor-base result, SHA-256
  `dcaa1f66f5f89757a3670a4a8e158ca6d5c360b7ca119d384d2c4a0b2794c2c8`.
- Round 19 selector result, SHA-256
  `3096540621408e4a48cfa18963ad01da9686d3527ee26776c8b6cf6f45a71114`.
- Round 24 negation-folded collision result, SHA-256
  `d1e2fe33e1abbae8a9bf6885ca7cfef8012233363f4cd43e0d30d468a381fc9c`.
- Round 311 result, SHA-256
  `9033d28dca7c42a0464ae11e5eb06052c0501cc79533b96a13dd148e7dfd0953`.
- Curve and factor-base identities:
  `icv1-fp256-t89188191154553853111372247798585809583-f188c491`,
  `FB1h2f8621cda105`, full FB1 digest
  `2f8621cda105a51f8703dabe38dec87f3715ab855770acd7b2af0a14906f9e42`,
  point-stream digest
  `70ab21cbda76205ff43ab6bf85c68607be9f6b89df6399e6fa98e3489b888cf1`,
  131,458 columns, and depth-18 root exponent
  `0x2b6fdc73dc04e7667129`.

Reject every dependency-hash, schema, curve, FB1, point-stream, column-count,
terminal, or replay mismatch.

## Complete hierarchy domain

For a factor-base abscissa `x`, define

```text
u_0(x)=x,
u_(h+1)(x)=u_h(x)^2-2 mod p.
```

At level `h`, two columns have minimal Dickson coalescence `h` when their
`u_h` values agree and their `u_(h-1)` values differ.  Completely enumerate
every such unordered column pair for `h=1,...,7`.  For each pair `(i,j)`, form

```text
A_(i,j,+) = P_i + P_j,
A_(i,j,-) = P_i - P_j.
```

Canonicalize each nonidentity atom up to group negation, but retain the sign
needed to reconstruct its sparse row.  The complete screen covers:

1. every collision between two enumerated atoms (weight at most four);
2. every match between an atom and `[m]P_k` for
   `m in {1,2,3,4,8,16}` (weight at most three); and
3. every collision among the complete small-multiple tables themselves
   (weight at most two).

Deduplicate sparse relations by primitive integer coefficients and global
sign.  Replay every relation in the P-256 group.  Compute its rank modulo both
`2^61-1` and the P-256 subgroup order; promotion credit is the smaller exact
rank.  If the relation matrix exceeds the native exact-rank implementation's
declared materialization bound, preserve every row and report rank as unset
rather than substituting a probabilistic estimate.

## Deeper diagnostic, excluded from promotion

For each level `h=8,...,18`, generate exactly 4,096 distinct unordered column
pairs whose minimal coalescence is `h`, using the SHA-256 domain

```text
round312/deep-pair/v1/<h>/<counter>
```

to select column indices.  Screen both atom signs against each other and the
complete small-multiple tables.  These 45,056 pairs are a deterministic
diagnostic only.  Unsampled pairs are discarded explicitly; no absence claim,
rank reduction, yield, exponent, or exhaustive-search projection may be based
on them.  Any sampled hit is replayed and reported, but requires a separately
preregistered complete follow-up before promotion.

## Controls and exactness

1. Exhaust a small prime cyclic group with a complete binary hierarchy,
   comparing the optimized sorted collision census against a direct reference
   for every local pair, sign, and registered multiplier.  Include planted
   independent identities and require zero false positives and false negatives.
2. Build a deterministic scalar-labelled P-256 hierarchy with planted pair,
   atom/multiple, and multiple/multiple identities.  Detect and replay every
   planted identity.  It is a detector control and receives no factor-base
   rank credit.
3. Independently rebuild `FB1h2f8621cda105` and validate every P-256 input
   point before screening.
4. Emit counts per hierarchy level, child widths, atom records, identity
   classes, survival ratios, field/group operations, logical RAM, disk output,
   false positives, false negatives, relation rank, and deterministic stream
   digests.

## Cost and promotion accounting

Let `r` be the certified rank of complete-scope static relations and
`K=131458-r`.  Report:

- the exact free-oracle K-th collision boundary and ratio to rho;
- the unchanged Round-19 S17 relation probability and degree evidence;
- 138,031-row accounting with only certified static-rank credit;
- sparse linear algebra, verification, target recovery, RAM, and disk traffic
  separately; and
- the Round-19 unmodified baseline in the same table.

The complete local screen is preprocessing and is charged to cold one-target
cost.  A zero-rank outcome leaves the prior end-to-end projection unchanged.
Any nonzero rank is an advance only if total `S` falls; a local identity count
alone is not a speedup.

Promotion requires all of:

- zero errors on every complete checked instance and exact replay of every
  reported identity;
- certified effective log rank and structured residual degree at most five;
- projected collection below `2^120` group/field-equivalent operations;
- projected cost below `2^103` per usable dynamic row;
- peak projected materialized storage below `2^50` bytes;
- complete cold `S <= 1.3`, with every phase included; and
- no credit from the sampled deep branch.

Only if all gates pass may the runner attempt an actual unplanted full-depth
P-256 relation.  Otherwise publish the measured rank/exponent obstruction and
do not claim parity.

## Hypotheses and stop condition

- H1: Dickson-local group identities may give non-scalar static log transport
  while preserving the degree-four S17 system.
- H2: if the complete static rank is too small, the K-th collision term remains
  above rho even with a free decomposition oracle.
- H3: no collision in the complete levels-1--7 domain is an exact negative for
  that domain only, not for all 131,458-column pairs.
- H4: sampled deeper failures are diagnostics and cannot close the unexamined
  domain.

Commit this protocol, then commit the native Rust runner before executing it.
Preserve deterministic JSON, isolation receipts, exact commands, failures,
tests, a requirement ledger, the one-unit boundary table, transfer-assessment
companion, dashboard update, and stacked pull request.  The transfer skill's
main workflow was available, but its linked methodology and assessment JSON
template were unavailable; use the repository's existing typed-assessment
schema and record that limitation.
