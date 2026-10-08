# P-256 affine-invariant factor-base screen, round 27: protocol

Date frozen: 2026-10-06

Rounds 25 and 26 tested one scalar orbit at the closest attainable width.
This round broadens the screen to every proper signed factor base preserved by
an affine known-log action

```text
f(P) = [a]P + [b]G
```

and to unions of scalar orbits, while keeping the registered P-256 curve and
signed distinct-column S17 relation model.  The purpose is to determine
whether additive progressions, affine recurrences, or a multi-orbit base can
simultaneously provide a cheap memoryless action and enough columns.  A
known-log representation or a generic group walk must not be relabelled as an
index-calculus advance.

This is a bounded factor-base obstruction screen, not a P-256 discrete-log
attack.

## Frozen curve and dependencies

- curve:
  `icv1-fp256-t89188191154553853111372247798585809583-f188c491`;
- relation model: signed, distinct-column S17 with balanced 8+9 colours;
- comparison base: `FB1h2f8621cda105`, 131,458 columns;
- maximum comparison width: 139,592 columns, the Round 25 closest one-orbit
  width;
- useful-row target: 138,031;
- Round 25 result SHA-256:
  `dcb191f7d5d7a660e3ce13128b25ce1c017c7d89ec2d01a242ab536651647faf`;
- Round 26 result SHA-256:
  `463c100dd951e2a7a7a419118f84bde90fc419358d8279a430a5a0e6191db029`.

The binary must reject either dependency if its hash or frozen curve identity
differs.

## Affine classification

Work on logarithm coefficients in the prime field `F_n`.  A signed factor
base is a nonempty set `S = -S`, with zero excluded.  Prove and replay the
following dichotomy.

1. If `f(x)=a*x+b` permutes a proper signed set and `b != 0`, then conjugating
   negation by `f` supplies the reflection `x -> -x+2b`; composing it with
   negation supplies the nonzero translation `x -> x+2b`.  Since `n` is an
   odd prime, that translation is transitive on `F_n`, contradicting
   properness.  Thus `b=0`.
2. With `b=0`, every proper invariant signed factor base is a union of
   multiplication-by-`a` orbits on `F_n^*/{+-1}`.  The column-orbit width is
   `ord(a)/2` for even order and `ord(a)` for odd order.
3. `a=1` is the identity and `a=-1` is the identity after sign folding; neither
   mixes columns.  Every nontrivial candidate therefore has order at least
   three and folded width no larger than the factor-base width.

Also construct the comparison arithmetic progression
`{[u+jv]G : 0 <= j < B}`.  Replay every internal translation edge and the
wrap defect.  The wrap is a constant group translation only if `B*v=0 mod n`,
which is impossible for `0 < B < n` and nonzero `v`.  Do not count a
non-wrapping prefix as a cyclic or exhaustive walk.

## Complete scalar-action census

Reuse and independently verify Round 25's complete factorisation of `n-1`
and primitive root.  Enumerate all 10,240 divisors.  Retain every order whose
negation-folded orbit width is at most 139,592, excluding orders one and two.
For every retained order, enumerate every exact-order element of its unique
subgroup.  Fold `a` and `-a` only when they induce the same column action, and
compute both:

- binary double-and-add cost (`bits + Hamming weight - 2`); and
- the rigorous doubling-only lower bound (`bits - 1`).

Record exact-order elements visited, duplicate folded actions removed, the
minimum-cost action for every order, and a SHA-256 digest of the ordered
census.  No sampled generator search may be labelled complete.

For each retained width, choose every integer orbit count `K` for which
`131458 <= K*width <= 139592` and enough distinct signed orbits exist.  Emit
the nondominated frontier over total columns, independent unknown anchors
`K`, binary action cost, and doubling lower bound.  Separately identify the
minimum-cost action in the width band and the minimum-anchor row.

## Exact controls and native identity

Materialise the selected minimum-cost width-band candidate from deterministic
SHA-256-labelled anchors, rejecting and advancing an anchor if its complete
signed orbit intersects an earlier orbit.  Sort canonical columns by the
repository 33-byte point key and form the native `ecbench.factor_base/v1`
FB1 identity.  An optional dump, if emitted, must use
`ecbench.factor_base_dump/v1-wide`.

Replay every selected orbit edge, wrap edge, curve-membership check, and
uniqueness-up-to-sign check in the P-256 group.  Independently replay at least
4,096 hash-selected affine actions on exact points.  Report all group and
field operations, bytes, false positives, and false negatives.

If materialising the complete minimum-cost candidate would exceed a measured
resource cap, emit the exact algebraic census plus a preregistered prefix and
classify the native base as incomplete; do not promote it.

## End-to-end accounting and gates

Use `1.3*sqrt(n)` as the repository rho reference.  A coalescing global action
must fit below Round 26's maximum `1.0369824472221065` group operations per
sample even before relation solving.  Count the target-affine correction on
the right colour.  Bases with known anchor logs are generic representation or
claw searches.  Bases with unknown anchors retain at least `K` independent
log classes before relations; local degree or cheap transport alone is not a
complete attack.

Promotion requires all of:

- zero false positives and false negatives on complete checked instances;
- exact replay of every materialised transport edge and control;
- a complete, coalescing action costing at most 1.036983 group operations per
  sample including target correction;
- structured residual degree of regularity at most five;
- projected relation collection below `2^120` operations and `2^103` per
  usable row;
- peak projected materialised storage below `2^50` bytes;
- complete one-target cost below matched rho; and
- a selector and anchor treatment demonstrably outside generic rho/claw
  search.

Attempt no full-depth unplanted P-256 relation unless every gate passes.  If
no candidate passes, publish the complete action count, Pareto frontier,
native replay evidence where feasible, and the dominant obstruction.  Emit
deterministic JSON with hashes, update the canonical dashboard, and prepare a
stacked PR.
