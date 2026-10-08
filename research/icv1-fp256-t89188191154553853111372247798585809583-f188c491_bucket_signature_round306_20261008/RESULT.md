# P-256 relation-bucket signature rank, round 306: result

## Verdict

Many decompositions in one target or walk bucket do **not** provide many
independent post-static equations.  After quotienting a maximal homogeneous
relation space, every exact relation

```text
sum_i c_i P_i - [b]Q = [a]G
```

has two-variable signature `(c dot ell,-b)`.  Exact group replay in a fixed
`(a,b)` bucket forces `c dot ell=a+b*d`, so every decomposition in that bucket
has the same signature `(a+b*d,-b)`.  Differences between bucket rows are
static.  Arbitrary occupancy therefore adds at most one quotient dimension.

The complete order-seven census checked 19,200 full buckets containing
6,585,600 retained rows and all 451,200 unordered signature pairs.  It replayed
7,488,000 cyclic-group equations and found zero false positives, false
negatives, replay failures, or rank discrepancies.  Every 343-row bucket had
rank four after a rank-three static space.  A pair reached rank five exactly
when its two signatures were nonproportional.

The P-256 scalar-labelled control confirms the same obstruction.  A 327-row
same-signature system and a second proportional signature both remain rank 164
in 165 unknowns.  An independent generator signature `(1,0)` reaches rank 165.
All 984 P-256 equations replayed exactly, and all ranks agree over the subgroup
order and `2^61-1` in both row orders.

No non-scalar-labelled event emits the required two nonproportional signatures.
No relation was reported on `FB1hc72514a2a8d3` or the registered 17-term base
`FB1h2f8621cda105`.  Structured degree, relation yield, sparse linear algebra,
and complete recovery cost remain unset.  The optimistic two-event oracle is
already **1.446505 times rho**, and the executable boundary remains
**13.920747 times rho**.  The promotion gates fail, so no unplanted full-depth
P-256 relation was attempted.

## Requirement-to-evidence status

| requested direction or gate | status | Round-306 evidence or gap |
|:--|:--:|:--|
| many independent rows per event | rejected for one fixed bucket | all rows share one quotient signature, so occupancy does not increase post-static rank |
| complete finite reference | verified | 19,200 complete buckets and 451,200 signature pairs over all 400 projective profiles |
| exhaustive group replay | verified | 6,585,600 bucket rows and 902,400 pair representatives replayed |
| determinant versus full rank | verified | 48,000 dependent and 403,200 independent pairs, with zero discrepancies |
| P-256 same-signature multiplicity | verified for labelled control | 327 rows retain rank 164 over both fields and row orders |
| P-256 proportional second signature | verified for labelled control | `(d,-1)` and `(2d,-2)` retain rank 164 |
| P-256 independent second signature | verified for labelled control | `(d,-1)` and `(1,0)` reach rank 165 |
| zero false positives / negatives | verified | zero toy or P-256 discrepancies |
| non-labelled two-signature event | not established | no construction emits target coupling and an independent generator anchor |
| complete relation, decomposition, and recovery | not established | no new relation mechanism was constructed |
| structured residual degree `<=5` | not established | no residual relation system exists for this audit |
| collection below `2^120` | failed / unset | relation yield, solver, linear algebra, and recovery are absent |
| usable relation below `2^103` | failed / unset | no relation on either actual factor base was produced |
| storage below `2^50` | failed / unset | control matrices are not an attack projection |
| complete cost at or below rho | failed | even the free two-signature oracle is 1.446505 times rho |
| actual unplanted P-256 relation | correctly not attempted | promotion gates fail |

## Boundary table

The ratios are imported from the hash-pinned Round-304 and corrected Round-305
evidence.  Oracle rows omit the cost of constructing the required signatures.

| variant | ratio to rho | status |
|:--|--:|:--|
| Pollard rho | **1.000000** | complete reference |
| unquotiented geometry-defined `FB1hc72514a2a8d3` | **13.920747** | executable free-oracle boundary unchanged |
| free maximal static kernel plus two independent signature events | **1.446505** | optimistic oracle; construction and recovery omitted |
| registered 17-term `FB1h2f8621cda105` | **394.425280** | unchanged free-perfect-oracle comparison |

Class: **negative**.  No executable boundary moved.

## Complete toy signature census

For each of all 400 normalized nonzero log profiles in `(Z/7Z)^4`, the runner
enumerated all 48 nonzero `(a,b)` buckets and all `7^4=2,401` coefficient
vectors per bucket.  Exactly 343 coefficient rows satisfy each bucket
equation.

| quantity | exact count |
|:--|--:|
| projective log profiles | 400 |
| nonzero signatures per profile | 48 |
| coefficient vectors tested | 46,099,200 |
| complete buckets checked | 19,200 |
| retained bucket rows / equation replays | 6,585,600 |
| unordered signature pairs checked | 451,200 |
| pair-representative equation replays | 902,400 |
| dependent / independent pairs | 48,000 / 403,200 |
| false positives / false negatives | 0 / 0 |
| replay failures / rank discrepancies | 0 / 0 |

Every bucket had rank four in five unknowns.  Every dependent pair had rank
four, and every independent pair had rank five.  The 940,800 rank computations
materialized 17,798,400 rows and 88,992,000 entries, performing 4,569,600
modular inversions and 192,443,868 modular multiplications.  These are exact
algorithm counters, not a projected P-256 attack cost.

## P-256 control

The deterministic control reconstructs the corrected Round-305 164-column log
witness and target.  It is scalar-labelled and receives no construction credit.

| system | rows | columns | rank over both fields and orders | nullity | P-256 group replays |
|:--|--:|--:|--:|--:|--:|
| static homogeneous | 163 | 164 | 163 | 1 | 163 |
| fixed-target bucket | 164 | 165 | 164 | 1 | 164 |
| enlarged fixed-target bucket | 327 | 165 | 164 | 1 | 327 |
| target plus proportional signature | 165 | 165 | 164 | 1 | 165 |
| target plus generator signature | 165 | 165 | 165 | 0 | 165 |

Across 20 rank computations, the runner materialized 3,936 rows and 648,788
entries.  Exact group replay used 2,948 scalar multiplications, 752,444 scalar
bits, 399,536 set bits, and 2,782 point additions.  All 165 constructed points
are on the curve; all 984 equations replayed with zero failures.

## Scope and next obligation

The result is exact for fixed-signature buckets after a maximal homogeneous
quotient in a prime-order cyclic subgroup.  It does not bound an event that
emits multiple nonproportional signatures.

The weakest surviving target is therefore narrower than “many rows”: construct
and replay one non-scalar-labelled P-256 event that emits both target coupling
and an independent generator signature.  Only then is there a residual system
whose structured degree, relation yield, sparse solve, and recovery cost can be
priced.  Raw bucket occupancy must not be counted as independent rank.

## Isolation and artifact integrity

Canonical and independent runs were isolated, uncontended, and emitted
byte-identical result and assessment JSON.

| run | exit | wall | user | system | peak RSS | status |
|:--|--:|--:|--:|--:|--:|:--|
| canonical | 0 | 6.963504 s | 6.952598 s | 0.007926 s | 12,160 KiB | final, uncontended |
| independent | 0 | 7.152032 s | 7.144299 s | 0.007120 s | 12,160 KiB | byte-identical, uncontended |

| artifact | bytes | SHA-256 |
|:--|--:|:--|
| `bucket-signature-result.json` | 12,292 | `08e11b9c0c829210100730f01710d2b68d80702bcd5d1768b3e594ca88713a31` |
| `transfer-assessment.json` | 4,077 | `cf518acd5d1434907196d23b3d5214edf3c0909dec6b5fc213c47b22cfbd70fd` |
| `isolation.jsonl` | 2,497 | `31f7835385293d4cbae0fd1e1f83e1eb34f92758b942d046f6ff45cd2bee64cf` |
| `PROTOCOL.md` | 7,562 | `fc528d174f5b3e5f60e90d7b7d7025fe52b5b7479a32e2a141f485d8d0059106` |

The result semantic-evidence SHA-256 is
`74ba0bc692c471dda28ae8a3c81409011aa5dffbe303423bfe57964f9783ccfd`.
The transfer-assessment semantic hash is
`1c02cadfc0040e2699e95a1c4a900460ac1cb77cd15056d46b5e1eb8e228b17f`.

The transfer workflow's referenced `references/methodology.md` and
`assets/assessment-template.json` resources remained unavailable.  The typed
correspondence graph, obligations, controls, scope, and open obligation are
emitted directly, and the limitation is recorded in the assessment JSON.

## Reproduction

```bash
cargo test --bin p256_bucket_signature
cargo clippy --bin p256_bucket_signature -- -D warnings
cargo build --release --bin p256_bucket_signature
python3 tools/isolated_bench.py run --wait --cpus 4 \
  --label round306-canonical \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_bucket_signature_round306_20261008/isolation.jsonl -- \
  target/release/p256_bucket_signature \
  --round298 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_parity_escapes_round298_20261007/parity-escape-result.json \
  --round304 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_algebraic_quotient_rigidity_round304_20261008/algebraic-quotient-rigidity-result.json \
  --round305 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_homogeneous_anchor_round305_20261008/homogeneous-anchor-result.json \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_bucket_signature_round306_20261008/bucket-signature-result.json \
  --assessment research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_bucket_signature_round306_20261008/transfer-assessment.json
```
