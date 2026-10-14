# P-256 Dickson--Walsh global transport, round 313: result

## Verdict

Round 313 is a complete negative result for the preregistered global
Walsh-character family.  On the actual `FB1h2f8621cda105` factor base, the
runner reconstructs all `2^18 = 262,144` canonical Dickson-fibre slots, places
the 131,458 low-y factor-base representatives in their original slots, places
the group identity in all 130,686 non-lifting slots, and exhausts every one of
the 262,144 binary characters.

Every Walsh output is nonidentity.  No two outputs agree up to elliptic
negation, and no output agrees up to negation with `[m]P_j` for any of the
131,458 columns and any `m in {1,2,3,4,8,16}`.  Consequently the exact row set
is empty and the certified rank is zero in both row orders modulo both
`2^61-1` and the P-256 subgroup order.  The active torus indices have full
affine dimension 18, so this is not an accidental character-row alias.

The independent final run is byte-identical to the canonical result.  The
complete cyclic control has 18 optimized and 18 directly constructed
relations with zero false positives, false negatives, transform mismatches,
or replay failures.  The scalar-labelled P-256 control finds both planted
relation types and replays both in the actual group.

This family therefore removes no logarithm classes.  The most favourable
free-relation-oracle lower projection remains `S=512.752864`, or
`394.425280` times the fixed `S=1.3` rho reference.  The imported executable
S17 projection remains `2^148.080980` FME per usable relation and
`2^165.155613` FME for 138,031 rows; attack materialization remains
`2^148.926255` bytes.  The inherited structured residual degree remains 4,
but that local degree does not reduce the end-to-end cost.

No `ecbench` session was fabricated.  The complete family has rank zero and
does not define a logarithm-recovering candidate, while the current native
prime path accepts only `p < 2^62`.  Under the preregistered fail-closed gate,
a stage-only timing or an analytic P-256 record would not be a rho comparison.
No unplanted full-depth relation was attempted.

## Requirement-to-evidence status

| requested direction / gate | status | Round-313 evidence or gap |
|:--|:--:|:--|
| registered ICV1 curve and FB1 storage | verified | exact curve slug, `FB1h2f8621cda105`, 131,458 columns, 262,916 signed rows; native dump and dependency hashes verified |
| global variable-column family | complete | every factor-base column remains active in every applicable character row; no fixed 16-column atom set |
| all Dickson characters | complete | 262,144 / 262,144 characters evaluated; zero discarded or probabilistic branches |
| identity, collision, and small-multiple matching | complete | every output and all 788,748 `[1,2,3,4,8,16]P_j` values covered |
| zero false positives / false negatives | verified | cyclic control: 18 optimized = 18 independent direct-reference rows; 0 / 0 errors |
| exact P-256 replay | verified controls; vacuous actual set | 2 / 2 planted control rows replay; zero actual-base rows reported |
| exact static rank | verified | rank 0 in forward and reverse orders over both audit moduli |
| structured residual degree `<=5` | passes, unchanged | imported Round-19 residual maximum 4; no new degree claim |
| per usable relation below `2^103` | failed | `2^148.080980` FME |
| 138,031-row collection below `2^120` | failed | `2^165.155613` FME |
| attack storage below `2^50` | failed | `2^148.926255` bytes |
| measured cold rho parity | failed / unset | free-oracle floor is 394.425280 times rho; no eligible end-to-end candidate |
| `ecbench` comparison | correctly ineligible | rank-zero screen is not a DLP solver; native P-256 execution is unsupported |
| unplanted full-depth P-256 relation | correctly not attempted | promotion gates fail |

## Complete family and result

For every `s in F_2^18`, the runner computes

```text
W_s = sum_i (-1)^<s,i> P_i
```

with a complete group-valued Walsh--Hadamard transform.  Negation
canonicalization retains the orientation needed to reconstruct a row.  The
screen separately covers identity outputs, equal/opposite output pairs, and
equal/opposite small-multiple matches.

| quantity | exact count |
|:--|--:|
| fibre slots / active factor-base slots | 262,144 / 131,458 |
| non-lifting identity slots | 130,686 |
| active-index affine span dimension | 18 |
| characters exhausted / discarded | 262,144 / 0 |
| nonidentity / identity outputs | 262,144 / 0 |
| output collision candidates | 0 |
| output--small-multiple candidates | 0 |
| exact unique rows / certified rank | 0 / 0 |

The absence of matches is exhaustive only for this registered linear
character family.  It is not evidence that every nonlinear or target-coupled
nonlocal identity family is empty.

## Work, memory, and traffic

| measured or exact quantity | count |
|:--|--:|
| Walsh butterflies | 2,359,296 |
| Walsh group additions | 4,718,592 |
| small-multiple additions / doublings | 131,458 / 525,832 |
| batch inversions | 2 |
| normalization multiplications | 6,305,350 |
| canonical-key comparisons | 39,557,274 |
| actual relation replay additions / doublings | 0 / 0 |
| conservative logical materialization | 236,811,486 bytes (`2^27.819164`) |
| measured final canonical peak RSS | 335,468 KiB |
| logical dependency reads | 40,316 bytes |

The native isolator does not expose block-I/O counters, so observed physical
disk bytes are explicitly unset.  The logical read count covers the three
frozen JSON dependencies.  Result and isolation writes are visible as the
artifact sizes below and are not presented as block-device measurements.

## One-unit boundary and `ecbench` decision

The comparison unit is charged P-256 group-addition equivalents divided by
`sqrt(n)`.  Only the finite Walsh screen is measured here; attack values are
unchanged imported projections and remain labelled as such.

| variant | S or cost | ratio / gate | evidence class |
|:--|--:|--:|:--|
| Pollard rho | `S=1.3` | `1.000000` times rho | fixed repository reference |
| free relation oracle | `S=512.752864` | `394.425280` times rho | optimistic lower projection |
| executable S17 per usable relation | `2^148.080980` FME | exceeds `2^103` | imported projection |
| 138,031-row collection | `2^165.155613` FME | exceeds `2^120` | imported projection |
| attack materialization | `2^148.926255` bytes | exceeds `2^50` | imported projection |
| complete cold DLP | unset | no parity claim | no qualifying implementation or solve |

`ecbench` is deliberately not invoked for this candidate.  It can measure a
replayable word-size `ic.pipeline` method, but zero certified rank supplies no
such method, and the repository's prime executor cannot directly instantiate
P-256.  Timing the 14--18 second screen against rho would compare unlike
outputs and violate the one-target end-to-end boundary.

## Isolation and artifact integrity

The original protocol was committed as `210b03538`, its fail-closed gate
clarification as `783816ffd`, the runner as `169f3eadf`, and the corrected
disk-accounting runner as `5a739076e`, all before the corresponding final
executions.  Final runs were pinned to CPU 4 on an Intel Xeon Platinum 8573C
and were uncontended.

| final run | exit | wall | user | system | peak RSS | other CPU |
|:--|--:|--:|--:|--:|--:|--:|
| canonical-v2 | 0 | 14.763046 s | 14.097107 s | 0.656189 s | 335,468 KiB | 0.03 s |
| independent-v2 | 0 | 14.518494 s | 14.082220 s | 0.430390 s | 335,452 KiB | 0.05 s |

The isolation ledger also preserves two successful pre-correction runs.  They
produced the same mathematical result but carried a misleading zero-valued
disk field; they are retained rather than erased.  The final independent JSON
is byte-identical to the canonical JSON.

| artifact | bytes | SHA-256 |
|:--|--:|:--|
| `PROTOCOL.md` | 8,544 | `d790b7f66965accb21138b758ae8ffd8937bade56cbbf8aea6d42d6be0342b4e` |
| `p256_dickson_walsh_transport.rs` | 50,309 | `21d18b894cd4ccf0c65710ec4ab554a2127dfc5d699221fc9d6d3739f339a803` |
| `walsh-transport-result.json` | 6,215 | `add31b33b52cc74103ba01d149fc4b8be85cb98090388daf7a593dd904aadf76` |
| `isolation.jsonl` | 8,709 | `cc0b50bd9f157469af3dbcd1611925bddb230ccffd4f9f372447d180c54e55d5` |

## Validation history

- `cargo test --bin p256_dickson_walsh_transport`: 5 passed.
- `cargo test --lib indexed_view_preserves_the_frozen_dump`: 1 passed,
  3,570 filtered.
- strict target clippy passed with allowances only for the repository's
  pre-existing `nonminimal_bool` diagnostics and stale lint name.
- `git diff --check`: passed.
- The first independent-reference test correctly failed with 10 apparent
  false negatives: the direct reference counted pairs of identity outputs
  again as collisions, whereas identity rows are a separate exhaustive class.
  After making the classes disjoint, optimized and independent sets agree.
- An accidentally broad initial test command attempted to link the repository
  bin fleet and exhausted the debug-build filesystem.  Only generated
  `target/debug` artifacts were removed; the scoped tests above then passed.

## Reproduction

```bash
cargo test --bin p256_dickson_walsh_transport
cargo test --lib indexed_view_preserves_the_frozen_dump
cargo clippy --bin p256_dickson_walsh_transport -- \
  -D warnings -A clippy::nonminimal-bool -A unknown-lints
cargo build --release --bin p256_dickson_walsh_transport --bin isolated_bench
target/release/isolated_bench run --wait --cpus 4 \
  --label round313-canonical-v2 \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_walsh_transport_round313_20261014/isolation.jsonl -- \
  target/release/p256_dickson_walsh_transport \
  --round2 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_dickson_terminal_degree_round2_20261004/factor-base-result.json \
  --round19 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_s17_multilevel_selector_round19_20261005/selector-result.json \
  --round312 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_dickson_additive_rank_round312_20261013/dickson-additive-rank-result.json \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_walsh_transport_round313_20261014/walsh-transport-result.json
```

## Next obligation

Complete linear characters on the Dickson fibre now join local additive atoms
as a certified zero-rank direction.  More selector tuning cannot lower the
394.425280-times-rho free-oracle floor.  The next credible experiment must
preregister a genuinely different complete family--for example a nonlinear,
target-coupled character convolution with independently replayable equations--
and must first demonstrate nonzero primitive log rank.  If it cannot, it is
not eligible for `ecbench` and cannot be credited toward parity.
