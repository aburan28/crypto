# P-256 half-turn involution atoms, round 311: result

## Verdict

The Round-309 half-turn involution admits an exact target factorization, but
its two atom channels are an invertible change of variables rather than a
smaller factor base.  The complete actual-P-256 census constructs 32,862
plus atoms and 32,862 minus atoms.  All 65,724 are distinct up to negation,
the channels do not overlap each other or the parent base, and none of the
complete scalar screens `{2,3,4,8,16}` finds even a partial match.  The
independent logarithm dimension therefore falls by exactly zero.

For each involution cycle, the identities

```text
A = P + P_pi,              B = P - P_pi,
A + B = 2P,                A - B = 2P_pi
```

and both cross-multiplied coordinate identities replay exactly.  The target
also splits exactly as

```text
sum a_i A_i = (Q+G)/2,     sum b_i B_i = (Q-G)/2.
```

This correction matters for accounting, but not enough for rho parity.  A
favourable target-only arity-20-plus-20 diagnostic is `3.649761` times rho
and requires `2^136.603866` materialized bytes.  It is not S17 and assumes
the atom logarithms are already known.  Collecting 138,031 usable atom rows
has a deliberately favourable charged lower projection of `2^137.537319`
group-operation equivalents, or `571.578319` times rho before target
decomposition; the optimistic cold total is `575.228080` times rho.  The
projected cost per usable row is `2^120.462686`, and structured residual
degree remains unset.

No promotion gate passes beyond correctness.  No actual unplanted full-depth
P-256 relation was attempted, as preregistered.  Round 311 is a reproducible
negative result: exact target separability does not compress the unknown-log
dimension, and relation collection is the dominant obstruction.

## Requirement-to-evidence status

| requested direction / gate | status | Round-311 evidence or gap |
|:--|:--:|:--|
| actual ICV1 P-256 factor-base construction | complete | `FB1h9361efda2bb1`, 65,724 columns in repository wide format |
| all parent columns remain variable | complete | all 32,862 two-cycles produce both A and B atoms; aggregate transform rank 65,724 |
| exact factorization and group replay | verified | zero invalid points, rebuild failures, group-identity failures, coordinate failures, or target-split failures |
| actual dimension reduction | failed | 65,724 unique atom columns from 65,724 parents; reduction exactly zero |
| collision or known-log transport | failed | zero internal duplicates, cross-channel overlaps, parent overlaps, and scalar matches for 2, 3, 4, 8, 16 |
| zero false positives / negatives | verified | complete F7 census, 960,400 states, zero of each and zero factorization/rank discrepancies |
| exact replay of every reported relation | vacuous but honest | zero actual relations reported; the labelled target is explicitly a control |
| structured degree `<=5` | failed / unset | no reduced residual system exists from which to establish a degree |
| complete cost below `2^120` | failed | favourable cold lower projection is `2^137.546502` operations |
| cost per usable relation below `2^103` | failed | favourable charged projection is `2^120.462686` per row |
| peak materialized storage below `2^50` | failed | target lists alone are `2^136.603866` bytes |
| complete collection, sparse algebra, and recovery | failed | collection is projected only; omitted costs are listed in the JSON |
| no discarded branch counted exhaustive | verified | toy census is complete; P-256 collection projections are labelled non-exhaustive lower bounds |
| rho parity | failed | optimistic target-only stage `3.649761` times rho; optimistic complete total `575.228080` times rho |
| unplanted P-256 relation | correctly not attempted | promotion gates fail |

## Actual atom census

The runner parses and hash-checks the complete Round-309 wide dump, rebuilds
the `x -> -x` pairing independently, and constructs both atom points on the
actual P-256 group.

| quantity | value |
|:--|--:|
| parent columns / two-cycles | 65,724 / 32,862 |
| plus atoms / unique up to negation | 32,862 / 32,862 |
| minus atoms / unique up to negation | 32,862 / 32,862 |
| union unique up to negation | 65,724 |
| internal duplicates, either channel | 0 |
| cross-channel overlaps | 0 |
| atom-to-parent overlaps | 0 |
| infinity or invalid points | 0 |
| cycle, group, or coordinate replay failures | 0 |
| aggregate local-transform rank | 65,724 |
| independent-log dimension reduction | 0 |

The scalar screens perform a complete set comparison up to negation after
multiplication by each registered small scalar.  Every one of the 20 directed
channel counts is zero; no complete or partial transport was found.

The actual-curve work ledger records 65,724 atom additions, 131,448
identity-replay additions, and 328,620 small-scalar additions.  The separate
deterministic target control uses two additions and three scalar
multiplications and has zero split-replay failures.  It receives no relation
credit.

The emitted base is:

| field | value |
|:--|:--|
| factor base | `FB1h9361efda2bb1` |
| columns / signed points | 65,724 / 131,448 |
| point-stream SHA-256 | `9ce858c07f411f1f39f85c6d720a742176f20a730b95c8ebb8e6671e77067e6d` |
| wide-dump SHA-256 | `57de26ddb4722ad5962b67c05d3074bcb14a52bd86b12583ce0fae30ee037511` |

## Complete F7 control

The toy reference exhausts 400 nonzero projective log profiles and all
`7^4=2,401` coefficient tuples for each profile, for 960,400 states.  Direct
and factorized outputs agree for every state.  The inverse atom transform,
forward/reverse row rank, and augmented rank also agree everywhere.

| operation | count |
|:--|--:|
| direct dot terms | 7,683,200 |
| atom dot terms | 3,841,600 |
| rank matrices | 1,920,800 |
| factorization discrepancies | 0 |
| inverse-transform discrepancies | 0 |
| rank discrepancies | 0 |
| false positives / false negatives | 0 / 0 |

## One-unit boundary

`S` means group-addition-equivalent operations divided by `sqrt(n)`.  Pollard
rho remains fixed at `S=1.3`.  Only correctness census work was executed;
every cost below other than rho is a projection or imported boundary.

| variant | S | ratio to rho | evidence class |
|:--|--:|--:|:--|
| Pollard rho | 1.3 | 1.000000 | fixed repository reference |
| Round-294 registered-base free collision oracle | 512.752864 | 394.425280 | imported optimistic lower projection |
| two separated atom targets only | 4.744689 | 3.649761 | target-only projection, atom logs assumed |
| collection with free paired-output samples | 371.525908 | 285.789160 | favourable non-exhaustive lower projection |
| collection with two charged group sums | 743.051815 | 571.578319 | favourable non-exhaustive lower projection |
| optimistic cold total | 747.796505 | 575.228080 | lower projection; no recovered DLP |

The collection formula is the exact imported K-th collision mean
`T_K=sqrt(2n)*Gamma(K+1/2)/Gamma(K)`, using `K=69,016` for both channels.
That deliberately omits the additional stopping time needed to fill both
quotas, coefficient generation, distinct-support rejection, verification,
duplicate/rank loss, solver setup, and recovery.  It is favourable to the
candidate and is not an exhaustive search.

Sparse algebra contributes a projected 367,197,621,936 operations
(`2^38.417766`), far below relation collection and unable to change the
decision.  The target lists contain `2^130.246314` entries across both
targets, cost `2^136.603866` bytes at 82 bytes per entry, and incur
`2^137.603866` bytes of write-plus-read traffic if materialized.

## Correction to the Round-309 boundary

Round 309's fixed S17 product mean `2^-271.270334` remains exact.  The present
factorization does not turn those S17 rows into cheap relations.  It defines
a different direct-product family: two arity-20 atom equations that can
expand to as many as 80 original point terms.  Its two targets separate, so a
generic product-group target charge is not the right target-stage accounting
for that different family.  The replacement target-stage bound is still
above rho and, crucially, assumes a full set of 65,724 atom logs.

Thus the useful correction is scoped: separability lowers one target-only
stage, while the invertible atom relabelling preserves the log dimension and
leaves complete collection hundreds of times above rho.

## Isolation and artifact integrity

The protocol was committed as `77737a4c2` and the runner as `6c2d05901`
before execution.  All runs were pinned to CPU 4 of an Intel Xeon Platinum
8573C host and were uncontended.

| run | exit | wall | user | peak RSS | other CPU |
|:--|--:|--:|--:|--:|--:|
| canonical-v1 | 0 | 257.442546 s | 256.269330 s | 238,188 KiB | 0.80 s |
| independent-v1 | 0 | 306.210842 s | 304.840250 s | 238,188 KiB | 0.52 s |
| independent-fixedpath-v1 | 0 | 263.126881 s | 262.165408 s | 238,176 KiB | 0.33 s |

The first independent result records a different atom-dump path, so its JSON
path field, semantic hash, and byte count differ; its 33 MB atom dump is
byte-identical.  The fixed-path independent pass makes the stronger check:
both the complete result JSON and atom dump are byte-identical to the
canonical artifacts.

| artifact | bytes | SHA-256 |
|:--|--:|:--|
| `involution-atom-result.json` | 13,604 | `9033d28dca7c42a0464ae11e5eb06052c0501cc79533b96a13dd148e7dfd0953` |
| `eigenspace-atom-factor-base.json` | 33,021,588 | `57de26ddb4722ad5962b67c05d3074bcb14a52bd86b12583ce0fae30ee037511` |
| `isolation.jsonl` | 7,533 | `b4c2ce571f8ca7c1e539f57f9296663e732a0a728b064fa2380955d65bbc4d67` |
| `PROTOCOL.md` | 5,710 | `4680790b08ebde4b3be60ff5c0a49f819db78649ebd66ae3a4e5660a790d42bd` |

Result semantic evidence SHA-256:
`dac123ec4cedd4b5cb0615760c666ea7348b064d1e84ddf95d511d094da8ec8b`.

## Reproduction

```bash
cargo test --bin p256_involution_atom
cargo clippy --bin p256_involution_atom -- -D warnings
cargo build --release --bin p256_involution_atom --bin isolated_bench
target/release/isolated_bench run --wait --cpus 4 \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_involution_atom_round311_20261013/isolation.jsonl \
  --label p256-involution-atom-round311-canonical-v1 -- \
  target/release/p256_involution_atom \
  --round294 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_variable_arity_round294_20261007/variable-arity-result.json \
  --round309 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_dickson_permutation_round309_20261008/dickson-permutation-result.json \
  --factor-base research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_dickson_permutation_round309_20261008/half-turn-factor-base.json \
  --round310 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_permutation_dichotomy_round310_20261008/permutation-dichotomy-result.json \
  --atom-base-out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_involution_atom_round311_20261013/eigenspace-atom-factor-base.json \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_involution_atom_round311_20261013/involution-atom-result.json
```

## Next obligation

Further tuning of this atomization is not a route to parity: its local
transform is full-rank and the complete actual census found no collision or
small-scalar transport.  A credible continuation must introduce and verify a
genuinely non-invertible global log transport or a non-generic collector that
produces many independent rows per charged event.  It must price its retained
domain, recover the target, and establish structured residual degree at most
five; local target separability alone is insufficient.
