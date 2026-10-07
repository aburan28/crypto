# P-256 transfer-anchor and structured-multiplicity screen, round 299: protocol

Date frozen: 2026-10-07

## Question and boundary

Round 298 leaves two sharply defined possibilities: a correspondence into a
cover/Jacobian representation might change more than coordinates, or a
deliberately structured factor base might place enough independent rows in one
known-target bucket to remove the 164-row collection boundary.  This round
tests both without attempting a P-256 discrete logarithm.

The production identity remains
`icv1-fp256-t89188191154553853111372247798585809583-f188c491` with factor base
`FB1hc72514a2a8d3`, 164 folded columns and arity 109.  The corrected Round-298
result is frozen at SHA-256
`8718bbe26751ff0191164b0665618226a97556c64c6e2feeef72468a38206965`.

## Prime-order transfer theorem to certify

Let `H=<G>` have prime order `n`, and let `phi:H->A` be the group homomorphism
induced by a proposed lift, cover, isogeny, Jacobian map or other typed
correspondence on the relevant subgroup.

The runner/report must certify the following cases, keeping curve morphisms
separate from their induced group maps.

1. If `phi(G)=0`, then `phi` kills all of `H`; it is not recoverable transport.
2. Otherwise `ker(phi) intersect H` is trivial, `phi(H)` has order `n`, and
   `log_{phi(G)}(phi([x]G))=x mod n`.  Transport preserves the one-dimensional
   unknown scalar and cannot by itself create a factor-base logarithm quotient.
3. Every endomorphism of the cyclic image `phi(H)` acts as a scalar modulo
   `n`.  A known scalar orbit is generic log transport; an unknown anchor
   remains one independent class.
4. For same-target coefficient rows, homogeneous differences have the actual
   factor-base log vector in their kernel.  Their rank over the subgroup field
   is therefore at most `B-1`.  If the target logarithm is known and nonzero,
   `B` distinct coefficient rows can have inhomogeneous rank `B`; if the target
   logarithm is unknown, augmented rows determine at most a projective class
   and do not recover an absolute discrete logarithm.

The scope is group-homomorphic restriction to the P-256 order-`n` subgroup.
It does not refute a nonlinear decomposition algorithm or a correspondence
with separately proved, cheaply usable destination relations.  Such a claim
must supply explicit models, fields, dimensions, formulas, kernel and
exceptional-point audit, recovery, relation pullback, and complete paired cost.

## Prime-subgroup multiplicity experiment

Use the existing deterministic toy curve

```text
E/F_1151: y^2 = x^3 - 3x + 241,  #E(F_1151)=1192=8*149.
```

Find the first nonidentity cofactor-cleared point `H=[8]P`, verify `[149]H=O`,
and enumerate its complete order-149 subgroup.  Fold its 148 nonidentity
points by negation into 74 candidate columns.  Scalar labels are retained only
for audit; geometry-defined families may not use them for selection.

The primary cell is `(B,m)=(6,3)`.  Every admissible base has raw signed domain
`2^3*C(6,3)=160`, exactly 80 rows after global-negation normalization, and 75
canonical subgroup targets, for fixed distinct mean occupancy `80/75`.  For
every base, enumerate the entire domain, replay every raw row by exact elliptic
curve addition, and compute over `F_149`:

- maximum distinct bucket occupancy;
- maximum rank of coefficient rows for a known target;
- maximum homogeneous difference rank;
- maximum augmented `[row|-target]` rank;
- counts and histograms of buckets attaining ranks 5 and 6;
- the residual anchor dimension and all violations of the certified ceilings.

Compute the same coefficient ranks modulo `2^61-1` as a deliberately labelled
formal-rank control.  A formal rank that does not respect the order-149 log
kernel is not transferable evidence for a P-256 relation bucket.

Screen these deterministic factor-base families:

| family | selection rule | scope |
|:--|:--|:--|
| affine progression | every `(a,d)` in `{1,...,148}^2`, six folded terms `a+jd` | complete within family; scalar-defined generic control |
| geometric progression | every `(a,r)` in `{1,...,148}^2`, six folded terms `a*r^j` | complete within family; scalar-defined generic control |
| pair-sum closure | every `1<=a<b<c<=74`, folded `{a,b,c,a+b,a+c,b+c}` | complete within family; scalar-defined generic control |
| hash scalar | 4,096 SHA-256-selected six-subsets | deterministic sampled generic control |
| low-coordinate geometry | every six-subset of the 12 columns with least encoded `(x,y)` | 924 complete candidates; no scalar labels in selection |
| hash-coordinate geometry | every six-subset of the first 12 columns under `SHA-256("round299/geometry" || x || y)` | 924 complete candidates; no scalar labels in selection |

Deduplicate bases within and across families while preserving family
membership.  Report family proposal, admissible and unique counts rather than
silently treating the finite families as an exhaustive search of all
`C(74,6)` bases.

## Exhaustive reference and controls

Exhaust all `C(74,4)=1,150,626` four-column bases at arity two.  Each has 24 raw
and 12 normalized rows.  Compare the scalar-sum engine, coefficient
normalization, ranks and target labels against exact group replay for every
row.  The primary six-column screen must also replay every raw row; no sampled
replay is accepted.

For every result row labelled “relation,” require exact equality of the group
sum and canonical target.  Record false positives and false negatives,
candidate and row counts, scalar/field operations, group additions, peak RSS,
wall/user/system time, deterministic enumeration hashes and independent
byte-identical replay.

## Hypotheses and promotion boundary

- **H1 (anchor conservation):** every nonzero homomorphic transfer of the
  P-256 prime-order subgroup is injective and preserves discrete-log scalars;
  transport alone does not reduce the 164 factor-base classes.
- **H2 (structured multiplicity):** scalar-defined families can produce sparse
  high-multiplicity/full-inhomogeneous-rank buckets, but those bases arrive
  with generic scalar transport and leave one unknown anchor if the generator
  is unknown.
- **H3 (geometry control):** coordinate-defined families do not show a
  reproducible full-rank advantage over matched hash controls at the fixed
  mean near one.  A contrary complete toy result is retained as a candidate,
  not projected to P-256 without a construction and scaling experiment.
- **H4 (correctness):** all registered families and the all-base `B=4`
  reference are complete within their stated boundaries; every raw row
  replays; subgroup-rank ceilings hold; independent output is byte-identical.

No candidate is promoted unless it is geometry-defined without discrete-log
knowledge, reaches known-target coefficient rank `B` at distinct mean at most
`1.1`, has zero replay/rank-control errors on complete cells, supplies a
P-256 construction with structured residual degree at most five, and gives a
complete projected cost at or below rho.  The original gates also remain:
cost per usable relation below `2^103`, collection of 138,031 independent rows
below `2^120`, storage below `2^50` bytes, and no discarded branch counted as
exhaustive.  A toy success alone triggers a next obligation, not a P-256
relation attempt.

## Deliverables

Implement the screen in Rust, emit deterministic result and transfer-
assessment JSON, preserve the failed and successful bounded outcomes, update
the comparison dashboard, and prepare a pull request stacked on corrected
Round 298.  The transfer skill's main workflow is followed.  Its referenced
`references/methodology.md` and `assets/assessment-template.json` resources
were unavailable from the package at registration time, so the assessment
uses a repository-native schema and records that limitation.
