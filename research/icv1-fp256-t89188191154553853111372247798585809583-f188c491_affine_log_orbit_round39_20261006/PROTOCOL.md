# P-256 affine-log-orbit factor-base screen, round 39: protocol

Date frozen: 2026-10-06

Round 38 left a narrow symmetry target: an implicit factor base whose pair
compatibility is close to one class without paying the full pair table.  This
round screens the maximal algebraic version of that idea,

```text
P_i = [a_i]H + [b_i]G,
```

where `G` generates the prime-order P-256 subgroup, `H` is an unknown anchor,
and `a_i,b_i` are known scalars.  The common-translate subfamily `a_i=1` has
known pair differences and is therefore the most favourable possible case for
an implicit pair selector.

This is a factor-base obstruction experiment, not a P-256 discrete-log attack.
No full-depth unplanted relation may be attempted unless every registered gate
passes.

## Frozen curve, relation model, and dependencies

- curve:
  `icv1-fp256-t89188191154553853111372247798585809583-f188c491`;
- comparison factor base: `FB1h2f8621cda105`, 131,458 columns;
- relation model: 17 distinct variable columns, signed balanced `8+9`;
- useful-row target: 138,031 independent rows including allowance;
- round-19 dependency:
  `research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_s17_multilevel_selector_round19_20261005/selector-result.json`,
  SHA-256
  `3096540621408e4a48cfa18963ad01da9686d3527ee26776c8b6cf6f45a71114`;
- round-25 dependency:
  `research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_scalar_orbit_round25_20261006/scalar-orbit-result.json`,
  SHA-256
  `dcb191f7d5d7a660e3ce13128b25ce1c017c7d89ec2d01a242ab536651647faf`.

The binary must reject a dependency hash, schema, curve, factor-base identity,
or subgroup-order mismatch.

## Exact affine-log trichotomy

For a target

```text
R = [u]G + [v]H
```

and signed relation `R = sum epsilon_i P_i`, compute modulo the prime subgroup
order

```text
Delta_a = v - sum epsilon_i a_i
Delta_b = sum epsilon_i b_i - u.
```

The experiment must prove and replay the following exhaustive cases:

1. `Delta_a=0, Delta_b=0`: the equality is an anchor-independent coefficient
   identity and supplies no discrete-log information;
2. `Delta_a=0, Delta_b!=0`: the equality is impossible;
3. `Delta_a!=0`: any true equality returns the anchor log exactly as
   `h=Delta_b/Delta_a`; finding the relation is already a direct DLP witness.

No informative affine-log relation may be counted as an ordinary index-calculus
row.  No anchor-independent identity may be counted as target information.

For the common-translate family `a_i=1`, record that the balanced signed S17
coefficient is `sum epsilon_i = +/-1`.  Exact pair differences cancel `H`, but
that cancellation may not be credited as coverage of an informative target:
the target coefficient `v=+/-1` is the identity case, while every other `v`
makes a successful relation a direct solution for `h`.

## Deterministic exact controls

Run complete toy enumerations over the prime cyclic groups of orders 19, 23,
29, and 31.  Use deterministic affine columns and enumerate every anchor,
target coefficient pair, and signed distinct-column relation of the registered
toy widths.  Compare the trichotomy with direct modular equality and report
false positives and false negatives separately.  Also enumerate the
common-translate odd balanced analogue and verify every pair-difference law.

On P-256, derive one fixed-base anchor with a recorded SHA-256-selected scalar
for controls only.  Generate 4,096 deterministic, hash-selected 17-column
affine-log instances and their 4,096 one-unit negative controls.  Each positive
must be replayed by native group arithmetic, recover the exact anchor scalar
in the informative branch, and classify the identity branch correctly.  Each
negative must fail exact group replay.  Generate another 4,096 common-translate
pair controls and replay `P_i-P_j=[b_i-b_j]G` exactly.

Record candidate counts, case counts, modular inversions and multiplications,
P-256 additions/doublings, logical RAM, materialised bytes, false positives,
false negatives, and a deterministic semantic digest.  Re-running the binary
must reproduce byte-identical canonical JSON.

## Degree and end-to-end accounting

Keep degree claims scoped by factor base:

- import round 19's completed structured residual maximum of 4 only for
  `FB1h2f8621cda105`;
- leave the unsplit S17 degree of regularity unknown;
- do not transfer the degree-4 result to affine-log or scalar-orbit bases;
- for a materialised common-translate base of `B` distinct points, record the
  degree-`B` univariate support-polynomial boundary, not a residual degree;
- the affine-log family's structured residual degree remains unknown unless a
  completed measurement is produced in this round.

An informative relation is algebraically equivalent to solving the original
one-dimensional P-256 DLP.  Report that reduction exactly; do not turn the
repository's `1.3*sqrt(n)` rho reference into a proved lower bound, and do not
claim a measured runtime ratio from the algebraic identity.  Project relation
collection only if ordinary independent rows survive the trichotomy.

## Promotion and stop conditions

Promotion requires all of:

- zero false positives and false negatives on complete checked instances;
- exact group replay of every reported relation and pair identity;
- structured residual degree of regularity no greater than 5 for this factor
  base, not an imported result for a different base;
- projected relation collection below `2^120` field/group-equivalent
  operations and below `2^103` per usable relation;
- peak projected materialised storage below `2^50` bytes;
- an informative selector demonstrably outside generic rho/claw search; and
- no discarded branch counted as exhaustive coverage.

If every informative relation directly exposes the anchor DLP while all cheap
relations are identities, publish a reproducible negative result.  Do not
implement a full selector, benchmark an attack, or attempt an unplanted
full-depth P-256 relation after that stop condition fires.

