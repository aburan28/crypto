# P-256 affine-log-orbit factor-base screen, round 39: result

Date run: 2026-10-06

The strongest remaining exact-log-transport family does not produce an index
calculus improvement.  For every factor-base column

```text
P_i = [a_i]H + [b_i]G
```

and target `R=[u]G+[v]H`, a true signed relation has exactly two outcomes: it
is a coefficient identity that contains no information about the unknown
anchor, or it returns the anchor discrete log directly.  The latter is the
original one-dimensional DLP, not an independent relation row.

The especially favourable common-translate family `P_i=H+[b_i]G` does make
every pair difference known for free.  That apparent `K=1` pair compatibility
exists only because the anchor cancels.  In a balanced signed S17 relation the
anchor coefficient is `8-9=-1`: matching a target with coefficient `v=-1` is
an anchor-independent offset identity, while any other successful target
coefficient directly solves for the anchor log.

The screen is negative.  No selector was implemented or benchmarked and no
full-depth unplanted relation was attempted.

## Exact reduction

Define, modulo the prime P-256 subgroup order,

```text
Delta_a = v - sum(epsilon_i a_i)
Delta_b = sum(epsilon_i b_i) - u.
```

The group equality reduces exactly to `Delta_a h = Delta_b`, where `H=[h]G`:

| case | equality | information |
|:--|:--|:--|
| `Delta_a=0, Delta_b=0` | always true | predetermined coefficient identity; no information about `h` |
| `Delta_a=0, Delta_b!=0` | impossible | no relation |
| `Delta_a!=0` | true only for `h=Delta_b/Delta_a` | the relation is already a direct DLP witness |

This closes the entire affine-log family, not merely one choice of coefficients.
It does not assert a new generic-group lower bound or a measured rho ratio.
Rather, an informative relation and the original DLP solution are the same
algebraic witness.  Therefore zero ordinary independent rows survive for a
138,031-row relation matrix, and relation-collection projections are not
defined for this family.

## Complete toy controls

The binary exhaustively enumerated every anchor, target coefficient pair and
signed `2+3` distinct-column relation for both general affine-log and
common-translate columns in four prime cyclic groups.

| order | family | candidates | identities | impossible | informative true | informative false | pair checks |
|--:|:--|--:|--:|--:|--:|--:|--:|
| 19 | general affine | 68,590 | 190 | 3,420 | 3,420 | 61,560 | 0 |
| 19 | common translate | 68,590 | 190 | 3,420 | 3,420 | 61,560 | 380 |
| 23 | general affine | 121,670 | 230 | 5,060 | 5,060 | 111,320 | 0 |
| 23 | common translate | 121,670 | 230 | 5,060 | 5,060 | 111,320 | 460 |
| 29 | general affine | 243,890 | 290 | 8,120 | 8,120 | 227,360 | 0 |
| 29 | common translate | 243,890 | 290 | 8,120 | 8,120 | 227,360 | 580 |
| 31 | general affine | 297,910 | 310 | 9,300 | 9,300 | 279,000 | 0 |
| 31 | common translate | 297,910 | 310 | 9,300 | 9,300 | 279,000 | 620 |

Totals:

```text
candidates checked             1,464,120
informative inversions             51,800
pair-difference checks              2,040
anchor-recovery failures                0
pair-difference failures                0
false positives                         0
false negatives                         0
```

## Native P-256 replay

The native controls use 4,096 deterministic hash-selected S17 coefficient
sets with 17 distinct columns each.  Half are anchor-independent identities;
half are informative planted equalities.  Every factor-base point is
constructed as `[a_i]H+[b_i]G`, all 17 signed terms are accumulated directly,
and the result is independently compared with the aggregated coefficient
form.  A one-unit target perturbation supplies a negative control for every
instance.

```text
positive term-by-term relation replays        4,096 / 4,096
negative one-unit controls rejected            4,096 / 4,096
informative anchor scalars recovered            2,048 / 2,048
termwise/aggregate mismatches                       0
distinct-column failures                           0
common-translate pair identities replayed       4,096 / 4,096
false positives / false negatives                  0 / 0

native scalar multiplications                   184,321
native doublings                             44,920,099
native additions                             22,734,751
```

These operation counts are verification work, not a projected attack cost.
The deterministic instance-stream digest is
`03121fdbebf0a9c8079a0fec6d92f76b850ab3b7552c39a0c7e579546e8cdfbc`.

## Degree boundary

The audit resolves the apparent degree-gate inconsistency without transferring
evidence across factor bases:

- round 19 completed a structured residual maximum of **4** for
  `FB1h2f8621cda105`; that factor base passes the registered local structured
  residual `<=5` gate;
- the unsplit S17 degree of regularity remains unknown;
- round 25 and the later near-rho selector line use different scalar/low-delta
  factor-base families, so they correctly retain an unknown structured
  residual degree;
- the affine-log/common-translate family screened here also has no completed
  structured residual degree measurement.  A materialised 131,458-point
  coordinate support polynomial has degree 131,458.

Thus the affine-log degree gate fails.  The old degree-4 result remains valid
and scoped to the original Dickson base; it was never an end-to-end degree for
every subsequent factor base.

## Resource and gate status

Materialising 33-byte point keys plus two 32-byte coefficients per column
would take 12,751,426 bytes, comfortably below `2^50`.  That is the only
resource gate this family passes.

| requirement | status | evidence or gap |
|:--|:--|:--|
| zero false positives and false negatives | passed | 1,464,120 exhaustive toy candidates and 12,288 native P-256 relation/pair controls |
| exact group replay | passed | all 17 terms replayed for every positive relation; independent aggregate comparison |
| structured residual degree at most 5 | failed | unknown for this family; comparison-base degree 4 is not transferable |
| collection below `2^120` | failed | no ordinary independent relation rows survive |
| per usable row below `2^103` | failed | no usable ordinary row exists |
| materialised storage below `2^50` | passed | 12,751,426 bytes for points and affine coefficients |
| non-generic informative selector | failed | every informative equality is the original DLP witness |
| discarded branches counted as exhaustive | no | the algebraic trichotomy is exhaustive; no probabilistic coverage claim |

## Decision

Reject affine-log-orbit and common-translate factor bases as a route to rho
parity.  They achieve exact log transport and, in the common-translate case,
free pair differences, but the same symmetry makes every cheap relation
uninformative.  An informative relation solves the DLP directly and supplies
no relation-collection advantage.

The next viable factor-base search must preserve the original Dickson base's
measured structured residual degree at most 4 while adding reusable global
selector structure that does **not** make all factor-base logs affine functions
of one hidden scalar.  Exact batched pair generation over independently chosen
starts remains open; affine log transport is now closed.

## Reproduction

```bash
cargo test --bin p256_affine_log_orbit
cargo run --release --bin p256_affine_log_orbit -- \
  --round19 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_s17_multilevel_selector_round19_20261005/selector-result.json \
  --round25 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_scalar_orbit_round25_20261006/scalar-orbit-result.json \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_affine_log_orbit_round39_20261006/affine-log-orbit-result.json
```

An independent release replay reproduced the complete artifact byte for byte.
The semantic evidence digest is
`b2f35e27f91794594dc6850a96b1ab345ee4c9300ee5e41a6951390062e0b0a3`.

Canonical result artifact: `affine-log-orbit-result.json`, 9,211 bytes,
SHA-256
`189e48a45ab2919dd3e279bf058940230c4a0461c745fba42dffb26f785ba079`.

