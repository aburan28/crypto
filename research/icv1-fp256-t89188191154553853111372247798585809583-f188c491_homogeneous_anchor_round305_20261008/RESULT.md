# P-256 homogeneous factor-base anchor conservation, round 305: result

## Verdict

Free static factor-base relations cannot remove the last logarithm scale on
P-256, even when every replayable homogeneous identity is granted without
discovery cost.

For `P_i=[ell_i]G`, every homogeneous relation satisfies

```text
sum_i h_i P_i = O  =>  h dot ell = 0 mod n.
```

The nonzero vector `ell` is therefore in the matrix kernel, so a 164-column
static relation matrix has rank at most 163.  This conclusion does not depend
on whether the identities were found algebraically, combinatorially, or by a
non-homomorphic construction.

For `Q=[d]G`, every unanchored decomposition similarly obeys

```text
sum_i c_i P_i = Q  =>  (c,-1) dot (ell,d) = 0 mod n.
```

Even all such decompositions determine only the projective class of
`(ell,d)`: augmented rank is at most 164 in 165 unknowns.  A known nonzero
generator right-hand side is still needed to fix the global scale.

The exact controls establish two full-rank routes after a maximal static
kernel:

1. one target-coupling row plus one independent generator anchor; or
2. two independent generator-anchored target rows.

The complete order-seven census checked 400 projective log profiles, 960,400
coefficient vectors, and 412,400 group equations.  The P-256 control ranked
six systems over two fields in both row orders and replayed 984 group
equations.  There were zero false positives, false negatives, replay
failures, or rank discrepancies.

These are scalar-labelled controls, not a construction on
`FB1hc72514a2a8d3`.  No event supplying the required target coupling and
independent scale anchor was found.  The optimistic one-anchor boundary is
already **1.446505 times rho**, so the executable boundary remains
**13.920747 times rho**.  Structured degree and complete recovery cost remain
unset, and no unplanted full-depth relation was attempted.

## Correction retained in the record

Canonical v2 omitted the combined `fixed target + generator anchor` system
and consequently stated the requirement too strongly as two
generator-anchored rows in every case.  Post-run review caught the error
before publication.  The protocol amendment, withdrawn v2 result and
assessment, and all isolation receipts are preserved.

The corrected statement is that two independent **post-static equations**
are required.  One may be an unanchored fixed-target coupling, in which case
only one additional generator anchor is needed.  Canonical and independent
v3 runs include that control and are the sole basis of this report.

## Requirement-to-evidence status

| requested direction or gate | status | Round-305 evidence or gap |
|:--|:--:|:--|
| arbitrary static factor-base relations | bounded exactly | every homogeneous row annihilates `ell`; rank is at most `K-1` |
| non-homomorphic discovery mechanisms | included in bound | the proof depends only on exact group replay, not on how a row was found |
| fixed-target multiplicity | bounded exactly | every augmented row annihilates `(ell,d)`; one global scale remains |
| complete finite controls | verified | 400/400 projective profiles have ranks 3, 4, 4, 5, and 5 at the declared stages |
| exact P-256 group replay | verified for labelled control | 984/984 materialized equations replay; ranks 163, 164, and 165 match over both fields and row orders |
| zero false positives / negatives | verified | zero toy or P-256 discrepancies |
| non-labelled target-plus-anchor event | not established | the full-rank P-256 rows were derived from witness scalars |
| complete relation, decomposition, and recovery | not established | no new relation mechanism was constructed |
| structured residual degree `<=5` | not established | no residual relation system exists for this audit |
| collection below `2^120` | failed / unset | relation yield, solver, linear algebra, and recovery are absent |
| usable relation below `2^103` | failed / unset | no relation on the actual factor base was produced |
| storage below `2^50` | failed / unset | control matrices are not an attack projection |
| complete cost at or below rho | failed | the free one-anchor oracle is 1.446505 times rho |
| actual unplanted P-256 relation | correctly not attempted | promotion gates fail |

## Boundary table

The ratios are imported from the hash-pinned Round-304 evidence.  Oracle rows
omit construction unless a direct control is shown.

| variant | free-oracle ratio to rho | direct ratio to rho | status |
|:--|--:|--:|:--|
| Pollard rho | **1.000000** | **1.000000** | complete reference |
| unquotiented geometry-defined `FB1hc72514a2a8d3` | **13.920747** | unset | 164 independent log classes |
| free maximal homogeneous kernel | **1.446505** | unset | target coupling plus scale anchor still required |
| Round-25 unknown-anchor scalar-labelled control | **1.446505** | **11.572038** | direct control above rho |
| known-anchor scalar-labelled control | **0.964336** | **7.714692** | ideal oracle below rho; direct control above rho |
| registered 17-term `FB1h2f8621cda105` | **394.425280** | unset | unchanged comparison |

Class: **negative**.  No executable boundary moved.

## Complete toy census

The runner normalized every nonzero vector in `(Z/7Z)^4` by its first
nonzero coordinate, producing all 400 projective classes.  For every class it
enumerated all 2,401 coefficient vectors.

| quantity | exact count |
|:--|--:|
| raw nonzero log vectors | 2,400 |
| projective profiles | 400 |
| coefficient vectors tested | 960,400 |
| static rows retained | 137,200 |
| fixed-target rows retained | 137,200 |
| fixed-target difference rows | 136,800 |
| anchored rows replayed | 1,200 |
| cyclic-group equations replayed | 412,400 |
| false positives / false negatives | 0 / 0 |
| replay failures / rank discrepancies | 0 / 0 |

Every profile produced the same exact ranks:

| system | rows per profile | columns | rank | nullity |
|:--|--:|--:|--:|--:|
| complete static kernel | 343 | 4 | 3 | 1 |
| fixed-target augmented space | 343 | 5 | 4 | 1 |
| fixed-target differences | 342 | 4 | 3 | 1 |
| static + one anchored row | 344 | 5 | 4 | 1 |
| static + two independent anchored rows | 345 | 5 | 5 | 0 |
| fixed-target space + one generator anchor | 344 | 5 | 5 | 0 |

The 4,800 forward/reverse rank computations materialized 1,648,800 rows and
7,696,000 entries.  The operation ledger is frozen in the result JSON.

## P-256 control

The P-256 control hash-expanded 164 nonzero witness logs and one target log,
constructed all 165 points with the repository's exact group law, and checked
that every point is finite and on the curve.  It is explicitly
scalar-labelled and receives no construction credit.

| system | rows | columns | rank over both fields and orders | nullity | P-256 group replays |
|:--|--:|--:|--:|--:|--:|
| static homogeneous | 163 | 164 | 163 | 1 | 163 |
| fixed-target augmented | 164 | 165 | 164 | 1 | 164 |
| fixed-target differences | 163 | 164 | 163 | 1 | 163 |
| static + one anchored | 164 | 165 | 164 | 1 | 164 |
| static + two anchored | 165 | 165 | 165 | 0 | 165 |
| fixed target + one generator anchor | 165 | 165 | 165 | 0 | 165 |

Across the 24 rank computations the runner materialized 3,936 rows and
648,136 entries; the exact counters are retained in JSON.  Group replay used
2,462 scalar multiplications, 626,549 scalar bits, 318,166 set bits, and 2,293
point additions.  All 984 equations replayed with zero failures.

## Scope and next obligation

The result is exact for homogeneous static identities and unanchored
fixed-target decompositions in the prime-order P-256 subgroup.  It does not
rule out a mechanism that produces target coupling and an independent known
generator right-hand side together.

That is now the weakest surviving target: construct and replay one P-256
event, without secret scalar labels, that supplies both equations and then
price its relation yield, structured solver, sparse linear algebra, and
target recovery below rho.  Static symmetry or many unanchored decompositions
alone cannot pass.

## Isolation and artifact integrity

The failed initial dependency check and both withdrawn v2 receipts are
preserved.  Canonical v3 and independent v3 were isolated, uncontended, and
emitted byte-identical final result and assessment JSON.

| run | exit | wall | user | system | peak RSS | status |
|:--|--:|--:|--:|--:|--:|:--|
| initial | 1 | 0.006878 s | 0.006652 s | 0 s | 4,484 KiB | preserved one-ULP parser-check failure |
| canonical v2 | 0 | 3.205042 s | 3.196356 s | 0.008356 s | 9,872 KiB | withdrawn: missing combined control |
| independent v2 | 0 | 3.172061 s | 3.167780 s | 0.004074 s | 9,928 KiB | withdrawn: missing combined control |
| canonical v3 | 0 | 4.150124 s | 4.135883 s | 0.011281 s | 11,208 KiB | final, uncontended |
| independent v3 | 0 | 4.035391 s | 4.011924 s | 0.019690 s | 11,208 KiB | final, uncontended |

| artifact | bytes | SHA-256 |
|:--|--:|:--|
| `homogeneous-anchor-result.json` | 15,262 | `ecc749ad8edf3e90981f9fb816382d4f4a0913539e9a05e454b6ef8d0a6f4cd7` |
| `transfer-assessment.json` | 4,790 | `e27605d643f065270c5cad298e98690f6a9cea9383409aedd9f99682450f7fe5` |
| `isolation.jsonl` | 13,402 | `924ea2866db73dd0e4ac595e25e5412bb2989ce1c0c1d88fb09d8aed6f61a2b8` |
| `PROTOCOL.md` | 9,900 | `32e2e4dbfc0ba25dce7381ebca2232e3a3d5e8491edce229425c0ba34d6b1d7e` |
| withdrawn v2 result | 13,918 | `f43375c2ba1a43c2ea7cf32230cb83868bb91b7079c3efcec9a2c337ac1c7b7a` |
| withdrawn v2 assessment | 4,537 | `e65813bee8fea05907156f4e69cf6775534f0e4599ce5a2f1b97537e893616f9` |

The final result semantic-evidence SHA-256 is
`2dc50589d8ccf8dd09de36a3187f04babbac6eb4dce252e890014ffe577f1628`.
The transfer-assessment semantic hash is
`69bd838358bd2a5751de2b5fac6c56a7fa6969136448600b6b09cf1588bdf7db`.

The transfer workflow's referenced `references/methodology.md` and
`assets/assessment-template.json` resources remained unavailable.  The typed
correspondence graph, obligations, controls, scope, and open obligation are
emitted directly, and the limitation is recorded in the assessment JSON.

## Reproduction

```bash
cargo test --bin p256_homogeneous_anchor
cargo clippy --bin p256_homogeneous_anchor -- -D warnings
cargo build --release --bin p256_homogeneous_anchor --bin isolated_bench
target/release/isolated_bench run --wait --cpus 4 \
  --label round305-canonical-v3 \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_homogeneous_anchor_round305_20261008/isolation.jsonl -- \
  target/release/p256_homogeneous_anchor \
  --round295-result research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_dickson_union_round295_20261007/dickson-union-result.json \
  --round295-factor-base research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_dickson_union_round295_20261007/selected-factor-base.json \
  --round298 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_parity_escapes_round298_20261007/parity-escape-result.json \
  --round299 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_transfer_anchor_round299_20261007/transfer-anchor-result.json \
  --round304 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_algebraic_quotient_rigidity_round304_20261008/algebraic-quotient-rigidity-result.json \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_homogeneous_anchor_round305_20261008/homogeneous-anchor-result.json \
  --assessment research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_homogeneous_anchor_round305_20261008/transfer-assessment.json
```
