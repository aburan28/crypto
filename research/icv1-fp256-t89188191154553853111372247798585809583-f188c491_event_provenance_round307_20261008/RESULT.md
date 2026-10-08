# P-256 one-event relation-provenance rank, round 307: result

## Verdict

Replayable post-processing cannot turn one primitive relation event into two
independent quotient equations.

Write an exact equation as the homogeneous row

```text
r = (c_0,...,c_(K-1),-b,-a)
```

for `sum c_i P_i-[b]Q=[a]G`, and let `W0` be the span of all independently
verified equations known before the event.  Scalar multiplication, group
addition/subtraction, column transport, and affine replay using identities
already charged to `W0` produce rows only in

```text
W0 + span(r).
```

Their syntax and right-hand sides may differ, but their primitive quotient
rank increment is at most one.  A transport identity outside `W0` or a second
group equality is new primitive information and must be counted as such.

The complete order-seven census enumerated and replayed all 6,722,800 exact
augmented equations across 400 projective log profiles.  Every one of the
3,200 complete 2,058-row projective-line closures had rank four over the
rank-three static space.  Every one of the 11,200 distinct-line pairs had rank
five.  There were zero partition, replay, rank, false-positive, or
false-negative discrepancies.

The P-256 scalar-labelled control generated 1,024 distinct rows `u*r+w` with
nonzero, varying right-hand sides.  The 1,187-row combined system retained
rank 164 in 166 homogeneous columns over both coefficient fields and both row
orders.  A proportional second primitive also retained rank 164.  An
independent generator equation reached rank 165.  All 1,844 P-256 equations
replayed exactly.

No non-scalar-labelled mechanism supplies two independent primitive equations.
No relation was reported on `FB1hc72514a2a8d3` or the registered 17-term base
`FB1h2f8621cda105`.  Structured degree, relation yield, sparse linear algebra,
and complete recovery cost remain unset.  The optimistic two-primitive oracle
is already **1.446505 times rho**, and the executable boundary remains
**13.920747 times rho**.  The promotion gates fail, so no unplanted full-depth
P-256 relation was attempted.

## Requirement-to-evidence status

| requested direction or gate | status | Round-307 evidence or gap |
|:--|:--:|:--|
| many rows from one event | bounded exactly | all post-processing using `W0` and one primitive remains in `W0+span(r)` |
| complete finite reference | verified | every exact augmented equation over all 400 projective order-seven profiles |
| full one-line closure | verified | 3,200 closures, each containing all 2,058 nonstatic rows on that quotient line |
| two-line independence | verified | all 11,200 distinct pairs add exactly two dimensions |
| exact P-256 group replay | verified for labelled control | 1,844/1,844 equations replay over five stages |
| zero false positives / negatives | verified | zero partition, replay, or rank discrepancies |
| P-256 derived-row multiplicity | verified for labelled control | 1,024 distinct derived rows retain one quotient dimension |
| non-labelled two-primitive event | not established | no construction independently certifies target coupling and a generator anchor |
| complete relation, decomposition, and recovery | not established | no new relation mechanism was constructed |
| structured residual degree `<=5` | not established | no residual relation system exists for this audit |
| collection below `2^120` | failed / unset | relation yield, solver, linear algebra, and recovery are absent |
| usable relation below `2^103` | failed / unset | no relation on either actual factor base was produced |
| storage below `2^50` | failed / unset | control matrices are not an attack projection |
| complete cost at or below rho | failed | the free two-primitive oracle is 1.446505 times rho |
| actual unplanted P-256 relation | correctly not attempted | promotion gates fail |

## Boundary table

The ratios are imported from hash-pinned Round-304 evidence.  Oracle rows omit
construction and complete recovery unless a direct control is named.

| variant | ratio to rho | status |
|:--|--:|:--|
| Pollard rho | **1.000000** | complete reference |
| unquotiented geometry-defined `FB1hc72514a2a8d3` | **13.920747** | executable free-oracle boundary unchanged |
| one unknown anchor plus two independent primitive events | **1.446505** | optimistic oracle; post-processing cannot reduce primitive count |
| known anchor plus one primitive event | **0.964336** | ideal oracle below rho; direct labelled control is 7.714692 times rho |
| registered 17-term `FB1h2f8621cda105` | **394.425280** | unchanged free-perfect-oracle comparison |

Class: **negative**.  No executable boundary moved.

## Complete toy provenance census

For every normalized nonzero `ell in (Z/7Z)^4`, the runner chose all four base
coefficients and `b`, derived the unique exact `a=c dot ell-b*d`, and replayed
the resulting six-coordinate homogeneous row against `(ell,d,1)`.

| quantity | exact count |
|:--|--:|
| projective log profiles | 400 |
| exact equations enumerated / replayed | 6,722,800 |
| static rows | 137,200 |
| nonstatic rows | 6,585,600 |
| complete projective-line closures | 3,200 |
| rows per line closure | 2,058 |
| distinct-line pairs | 11,200 |
| false positives / false negatives | 0 / 0 |
| replay / partition / rank discrepancies | 0 / 0 / 0 |

Every line closure had rank four; every distinct-line pair had rank five.  The
28,800 forward/reverse rank computations materialized 115,248,000 rows and
691,488,000 entries, performing 137,600 modular inversions and 1,860,539,016
modular multiplications.  These exact counters certify closure completeness;
they are not a projected P-256 attack cost.

## P-256 provenance control

The deterministic control reconstructs the Round-305/306 164-column witness
and target.  It is scalar-labelled and receives no construction credit.

| system | rows | homogeneous columns | rank over both fields and orders | quotient increment | P-256 group replays |
|:--|--:|--:|--:|--:|--:|
| static pre-event space `W0` | 163 | 166 | 163 | 0 | 163 |
| one primitive target event | 164 | 166 | 164 | 1 | 164 |
| 1,024-row derived closure | 1,187 | 166 | 164 | 1 | 1,187 |
| proportional second primitive | 165 | 166 | 164 | 1 | 165 |
| independent second primitive | 165 | 166 | 165 | 2 | 165 |

The 1,024 derived rows are distinct; their digest is
`5b11496f8fae9b004cda33803cabc3cd9e2deaf5b69c526aa00662e8263a48c7`.
Across 20 rank computations, the runner materialized 7,376 rows and 1,224,416
entries.  Exact group replay used 5,905 scalar multiplications, 1,504,906
scalar bits, 748,294 set bits, and 4,711 point additions.  All 165 constructed
points are on the curve; every equation replayed with zero failures.

## Scope and next obligation

The result is exact for the linear closure of one primitive group equation
modulo independently replayed pre-event knowledge.  It covers scalar
multiplication, group-law combinations, and affine/column transports whose
identities are already charged to `W0`.

It does **not** rule out a nonlinear search that genuinely discovers a second
independent equality.  Such an equality is a second primitive, even if both
are emitted by one software invocation or share intermediate computation.
The weakest surviving target is therefore a joint P-256 event that
independently certifies two quotient lines without secret scalar labels.  Its
shared work, residual degree, yield, sparse solve, and recovery must then be
priced below rho.

## Isolation and artifact integrity

Canonical and independent runs were isolated, uncontended, and emitted
byte-identical result and assessment JSON.

| run | exit | wall | user | system | peak RSS | status |
|:--|--:|--:|--:|--:|--:|:--|
| canonical | 0 | 21.998703 s | 21.868033 s | 0.116309 s | 26,592 KiB | final, uncontended |
| independent | 0 | 22.110082 s | 21.937659 s | 0.169746 s | 26,528 KiB | byte-identical, uncontended |

| artifact | bytes | SHA-256 |
|:--|--:|:--|
| `event-provenance-result.json` | 12,817 | `3751468a6708d1851598e366bfb99d0795d2da251377a3a1bbf7918cbb36034d` |
| `transfer-assessment.json` | 3,874 | `1d9b0044e50bc710ccac33d3cc6d0c3b369c50fbeff9a94e8766a484fa4fa067` |
| `isolation.jsonl` | 2,375 | `a9d3838df6f3fbf429dd0f46280c043088b92cdfe731596999b0683ca6a6e427` |
| `PROTOCOL.md` | 8,270 | `72e24549c0976753ff0feadf9385fc60d7d9e3e6e816cd6fefeeb4a9a96504f1` |

The result semantic-evidence SHA-256 is
`ce17b74cf28adb4294e00a92984290d7e76b6b405a1c1697e897642ced164d7b`.
The transfer-assessment semantic hash is
`9c079dff2fc5e8b9617ca49b9ed759b6678fa578f7d43b9e9791f7b9247bc0b4`.

The transfer workflow's referenced `references/methodology.md` and
`assets/assessment-template.json` resources remained unavailable.  The typed
correspondence graph, obligations, controls, scope, and open obligation are
emitted directly, and the limitation is recorded in the assessment JSON.

## Reproduction

```bash
cargo test --bin p256_event_provenance
cargo clippy --bin p256_event_provenance -- -D warnings
cargo build --release --bin p256_event_provenance
python3 tools/isolated_bench.py run --wait --cpus 4 \
  --label round307-canonical \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_event_provenance_round307_20261008/isolation.jsonl -- \
  target/release/p256_event_provenance \
  --round39 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_affine_log_orbit_round39_20261006/affine-log-orbit-result.json \
  --round304 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_algebraic_quotient_rigidity_round304_20261008/algebraic-quotient-rigidity-result.json \
  --round306 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_bucket_signature_round306_20261008/bucket-signature-result.json \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_event_provenance_round307_20261008/event-provenance-result.json \
  --assessment research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_event_provenance_round307_20261008/transfer-assessment.json
```
