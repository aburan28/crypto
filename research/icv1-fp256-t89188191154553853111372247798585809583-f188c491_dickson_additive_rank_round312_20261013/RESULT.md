# P-256 Dickson-additive static-rank screen, round 312: result

## Verdict

Round 312 finds no non-scalar logarithm transport in its preregistered
complete domain.  On the actual `FB1h2f8621cda105` factor base, it exhausts
all 4,185,775 unordered column pairs whose minimal Dickson coalescence is at
levels 1--7, materializes both `P_i+P_j` and `P_i-P_j`, and compares every
output against every other output and against the complete tables
`[1,2,3,4,8,16]P_k`.  The resulting 8,371,550 complete-scope atom records
produce zero identities and therefore a certified static relation rank of
zero.

The deeper levels 8--18 contribute exactly 4,096 deterministic hash-selected
pairs per level, or 90,112 signed atom records.  They also find zero
identities, but they are samples and receive no rank, cost, or exhaustive
coverage credit.  The complete order-101 control has 91 reference relations
with zero false positives or negatives, and the scalar-labelled P-256 control
replays all 70 detected relations.  There are no reported actual-base
relations, so exact replay of that empty set is vacuous and explicitly
labelled as such.

The registered base consequently retains all 131,458 independent log
classes within the screened scope.  Its imported structured residual degree
remains at most 4, which passes the local degree gate but is not a new degree
reduction.  The free-relation-oracle boundary remains `394.425280` times rho;
the executable S17 projection remains `2^148.080980` field-multiplication
equivalents per usable relation and `2^165.155613` for 138,031 rows.  Imported
attack storage remains `2^148.926255` bytes.  Those values fail the three
cost/storage gates by wide margins, so no unplanted full-depth P-256 relation
was attempted.

This is a reproducible negative result.  Dickson coordinate coalescence does
not, by itself, transport elliptic logarithms.  Further local tuning of the
same additive atom family cannot claim a boundary improvement without a
new, complete source of nonzero rank.

## Requirement-to-evidence status

| requested direction / gate | status | Round-312 evidence or gap |
|:--|:--:|:--|
| actual ICV1 curve and repository factor base | verified | `icv1-fp256-t89188191154553853111372247798585809583-f188c491`, `FB1h2f8621cda105`, 131,458 columns; native rebuild and dependency hashes verified |
| all columns remain eligible | complete in stated domain | every cross-child pair at minimal Dickson levels 1--7 is enumerated; no fixed atom columns |
| additive and subtractive atom channels | complete in stated domain | both `P_i+P_j` and `P_i-P_j` for 4,185,775 unordered pairs |
| small-multiple transport | complete | `[1,2,3,4,8,16]P_k` for every column, including all atom/multiple and multiple/multiple collisions |
| deeper Dickson hierarchy | sampled only | 4,096 hash-selected pairs at each level 8--18; no projection credit |
| zero false positives / negatives | verified | complete order-101 control: 91 optimized and reference relations, 0 / 0 errors |
| exact P-256 replay | verified controls; vacuous actual set | 70 / 70 labelled control relations replay; zero actual-base relations reported |
| exact static rank | verified | rank 0 in both row orders modulo `2^61-1` and the P-256 subgroup order |
| structured residual degree `<=5` | passes, unchanged | imported Round-19 local join degree 2 and residual maximum 4; unsplit S17 regularity remains unset |
| collection below `2^120` | failed | `2^165.155613` FME for 138,031 usable rows |
| per usable relation below `2^103` | failed | `2^148.080980` FME |
| materialized attack storage below `2^50` | failed | `2^148.926255` bytes imported from the unchanged S17 attack projection |
| no discarded branch called exhaustive | verified | levels 8--18 are labelled samples and excluded from every improvement projection |
| rho parity | failed | even a free relation oracle remains `394.425280` times rho; complete cold ratio remains unset |
| unplanted full-depth P-256 relation | correctly not attempted | promotion gates fail |

## Exact screen

The factor base is grouped by its verified Dickson chain.  A pair belongs to
level `d` exactly when its two columns first meet at depth `d`.  This makes
levels 1--7 disjoint and complete rather than overlapping neighbourhood
samples.

| level | coverage | unordered pairs | signed atoms | maximum left / right width |
|--:|:--|--:|--:|--:|
| 1 | complete | 32,862 | 65,724 | 1 / 1 |
| 2 | complete | 65,722 | 131,444 | 2 / 2 |
| 3 | complete | 131,847 | 263,694 | 4 / 4 |
| 4 | complete | 263,518 | 527,036 | 8 / 8 |
| 5 | complete | 527,123 | 1,054,246 | 15 / 15 |
| 6 | complete | 1,055,373 | 2,110,746 | 26 / 25 |
| 7 | complete | 2,109,330 | 4,218,660 | 44 / 48 |
| 8--18 | 4,096 pairs per level | 45,056 | 90,112 | sampled; not applicable |
| **total** | mixed as labelled | **4,230,831** | **8,461,662** | -- |

No atom-output bucket has multiplicity, no atom equals a screened multiple,
and no screened-multiple bucket has multiplicity.  The complete static row
set is therefore empty before and after deduplication.

## Correctness controls

The order-101 reference enumerates every representation pair directly.  Its
optimized collision detector returns the exact same 91 relations among 56
atom and 48 multiple records; their independent rank is 7.  False positives,
false negatives, and replay failures are all zero.

The separate scalar-labelled P-256 control uses eight columns.  It detects
70 relations from 56 atom and 48 multiple records, replays all 70 in the
actual group, and obtains rank 7 over both audit moduli.  This control proves
that the P-256 detector can find and replay the identity family that the
actual base lacks.

## Work and storage ledger

| measured quantity | count |
|:--|--:|
| Dickson squarings | 2,366,244 |
| hierarchy pairs examined | 4,230,831 |
| projective additions | 8,593,120 |
| projective doublings | 525,832 |
| batch inversions | 158 |
| normalization multiplications | 46,252,050 |
| atom records materialized | 8,461,662 |
| small-multiple records materialized | 788,748 |
| deterministic deep-pair hash draws | 16,607,132 |
| actual relation-replay additions / doublings | 0 / 0 |

The screen's conservative logical peak is 499,564,472 bytes
(`2^28.896096`): 372,313,128 atom bytes, 34,704,912 multiple bytes,
79,926,464 hierarchy-key bytes, and 12,619,968 input-point bytes.  Disk
traffic is limited to the JSON artifacts and isolation receipts.  This
measured diagnostic footprint is distinct from the failed
`2^148.926255`-byte end-to-end materialization projection.

## One-unit boundary

`S` means group-addition-equivalent operations divided by `sqrt(n)`, with
Pollard rho fixed at `S=1.3`.  Round 312 executes only the finite static-rank
screen.  The attack costs below are imported projections and are not
reclassified as measurements.

| variant | S or cost | ratio / gate | evidence class |
|:--|--:|--:|:--|
| Pollard rho | `S=1.3` | `1.000000` times rho | fixed repository reference |
| free relation oracle, unchanged base | `S=512.752864` | `394.425280` times rho | optimistic lower projection |
| executable S17 per usable relation | `2^148.080980` FME | exceeds `2^103` | imported projection |
| 138,031-row collection | `2^165.155613` FME | exceeds `2^120` | imported projection |
| materialized attack storage | `2^148.926255` bytes | exceeds `2^50` | imported projection |
| complete cold DLP | unset | no rho claim | no qualifying end-to-end execution |

Sparse linear algebra cannot rescue the candidate because the static screen
reduces neither the number of columns nor the required dynamic row count.
The same five-percent rank/duplicate allowance therefore remains 138,031
rows.

## Isolation and artifact integrity

The protocol was committed as `9ae69dff6` and the runner as `eed682890`
before execution.  Both runs were pinned to CPU 4 of an Intel Xeon Platinum
8573C host and were uncontended.  The independent result and transfer
assessment are byte-identical to the canonical files.

| run | exit | wall | user | peak RSS | other CPU |
|:--|--:|--:|--:|--:|--:|
| canonical-v1 | 0 | 88.310934 s | 86.133510 s | 584,576 KiB | 0.11 s |
| independent-v1 | 0 | 82.970713 s | 82.086678 s | 584,704 KiB | 0.11 s |

| artifact | bytes | SHA-256 |
|:--|--:|:--|
| `PROTOCOL.md` | 8,537 | `14bfec43f51afb190f5464050a812e992822d6d084f52cd089511b50355069bb` |
| `p256_dickson_additive_rank.rs` | 69,619 | `ce5d7dbb991ccbb1f5b22500934ad819b51e67af19cbdee96758eaddb3b74ec7` |
| `dickson-additive-rank-result.json` | 13,550 | `0cdec91f0097b498d1642c80c358aadac279c8a0bae951ed433381e1bf53db8e` |
| `transfer-assessment.json` | 3,618 | `53ba2358cd799ad511c7eae3232960868d3e3c2432b9d338996eea0ec13d00de` |
| `isolation.jsonl` | 4,966 | `f77895b7a42359a301a2bb15d1aa4fc88e25e742dce3e23ff2e73d8169df0167` |

Result semantic evidence SHA-256:
`f3a2fa9dd245ddea9bdde0224cee258c59feda60e0958cf04318e3bcb67c743a`.

The companion transfer assessment records semantic evidence SHA-256
`26aba4beefd29ba382f98e8f3b439ef9a31bbc0dc651d5c63ce841617e7c5cc2`.
The transfer skill's linked methodology and template resources were not
available in this environment; the assessment records that limitation and
does not borrow external claims.

## Reproduction

```bash
cargo test --bin p256_dickson_additive_rank
cargo clippy --bin p256_dickson_additive_rank -- \
  -D warnings -A clippy::nonminimal-bool -A unknown-lints
cargo build --release --bin p256_dickson_additive_rank --bin isolated_bench
RAYON_NUM_THREADS=1 target/release/isolated_bench run --wait --cpus 4 \
  --label p256-dickson-additive-rank-round312-canonical-v1 \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_dickson_additive_rank_round312_20261013/isolation.jsonl -- \
  target/release/p256_dickson_additive_rank \
  --round2 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_dickson_terminal_degree_round2_20261004/factor-base-result.json \
  --round19 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_s17_multilevel_selector_round19_20261005/selector-result.json \
  --round24 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_negation_cross_colour_round24_20261006/negation-result.json \
  --round311 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_involution_atom_round311_20261013/involution-atom-result.json \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_dickson_additive_rank_round312_20261013/dickson-additive-rank-result.json \
  --assessment research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_dickson_additive_rank_round312_20261013/transfer-assessment.json
```

The two lint allowances are scoped to this binary.  They accommodate one
pre-existing library `nonminimal_bool` diagnostic and the repository's stale
lint name under Rust 1.94.1; the new binary itself passes the resulting
`-D warnings` invocation.

## Next obligation

Levels 1--7 now have a complete negative certificate for this atom family.
Blindly raising that cutoff has quadratic growth and no observed signal.  A
credible next parity experiment must first supply either (a) a mathematically
defined nonlocal identity family with a completeness argument, or (b) a
collector that yields more than one independent row per charged group event.
It must demonstrate nonzero rank before importing any reduction in relation
count, and it must preserve exact target recovery and the one-unit rho
accounting.
