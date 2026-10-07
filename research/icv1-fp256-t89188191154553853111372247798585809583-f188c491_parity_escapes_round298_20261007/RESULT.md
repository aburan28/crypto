# P-256 multi-row and transport parity escapes, round 298: result

## Verdict

Neither registered escape reaches rho parity.

For the 164-class Dickson union, one collision event is the only event count
whose free-oracle constant is below rho: its ratio is 0.964088, while two
events already cost 1.446132 times rho.  A one-event attack must therefore
obtain all 164 independent rows from one target bucket, requiring at least 165
distinct normalized decompositions when the rows arise as collision
differences.  Before quotienting the global-negation duplicates, that is at
least 330 raw folded preimages.

The exact 109-of-164 signed domain has only

```text
D/n = 1.0051098993154817
```

decompositions per P-256 target on average.  Under an explicitly labelled
random-map null, the union bound that any canonical target bucket reaches 165
distinct normalized rows is below `2^-721.185`; the corresponding Poisson
expected-bucket diagnostic is `2^-727.636`.  These are model diagnostics, not a measurement or proof that
the elliptic subset-sum map is random.

Complete toy enumeration confirms the rank mechanism.  At mean occupancies
0.067 and 1.501 distinct rows per canonical target, maximum homogeneous rank
is only one and four.  Full column rank first appears at mean distinct
occupancy 35.377, and the 14-column system saturates only in the extremely
dense regime with mean distinct occupancy 858.479.  All 1,069,136 raw rows
replay exactly with zero failures and collapse to 534,568 distinct normalized
rows.  Raw relation multiplicity does not supply evidence for the required
sparse P-256 rank-164 event.

The typed transport audit also closes the measured families.  P-256
base-field endomorphisms act on `E(F_p)` only as scalars; the exact one-class
scalar orbit is generic rho with a known anchor and costs 1.446505 times rho
with an unknown anchor.  The certified degree-11 isogeny walk preserves the
order-`n` subgroup but does not reduce the factor-base quotient and worsens its
measured boundary.  Covers, extension-field lifts and Jacobian
correspondences remain explicit unknown obligations because no executable map
was supplied; they receive no improvement credit and are not declared
nonexistent.

No P-256 relation was attempted and promotion is rejected.

## Exact parity requirement

The exact collision law gives:

| independent collision events | free-oracle `S/rho` | parity |
|--:|--:|:--:|
| 1 | **0.964088** | yes, before every omitted cost |
| 2 | 1.446132 | no |
| 164 | 13.920747 | no |

Consequently an event rank of 82, which would halve the 164-class requirement,
still leaves two events and already fails parity.  The necessary target is one
event of rank 164.  Since the first preimage establishes a bucket and each
additional preimage can supply at most one independent collision difference,
the bucket occupancy must be at least 165.

That occupancy is after canonicalizing the target modulo negation and
deduplicating coefficient rows.  The raw signed enumerator produces each such
row twice, so the equivalent raw folded occupancy floor is 330.

The fixed P-256 domain and order are

```text
D = 116383775127750444427048016195109203089116454475801022361283226351238392053760
D/2 = 58191887563875222213524008097554601544558227237900511180641613175619196026880
n = 115792089210356248762697446949407573529996955224135760342422259061068512044369
```

The conditional mean occupancy of a nonempty bucket under the same random-map
null is only 1.585358.  The null-model tail is useful as a falsification target:
a real structured screen must exhibit and replay exceptional rank, not simply
assert that the Dickson support is nonrandom.

## Complete multi-row rank experiment

The experiment uses Round 296's deterministic support on
`E/F_1151: y^2=x^3-3x+241`.  Targets are folded modulo negation; coefficient
rows are negated whenever the target representative is negated.  Augmented
rank includes the target coefficient.  Homogeneous rank uses exact row
differences, or direct zero-sum rows for the identity bucket.  The toy group
has 1,192 points and 597 canonical targets under negation.  Raw rows are
replayed, then global-negation duplicates are removed before occupancy and
rank are computed; means include empty canonical targets.

| `m,B` | raw domain | distinct rows | nonempty targets | distinct mean / 597 | max distinct occupancy | max augmented rank | max homogeneous rank | full-column-rank buckets | replay failures |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| 3,5 | 80 | 40 | 37 | 0.067 | 2 | 2 | 1 | 0 | 0 |
| 5,8 | 1,792 | 896 | 487 | 1.501 | 5 | 5 | 4 | 0 | 0 |
| 7,11 | 42,240 | 21,120 | 597 | 35.377 | 64 | 12 | 11 | 560 | 0 |
| 9,14 | 1,025,024 | 512,512 | 597 | 858.479 | 975 | 15 | 14 | 597 | 0 |

Every occupancy and rank histogram is preserved in the result JSON.  The
enumerator charges 9,530,096 logical candidate group additions, 90,592
precomputed transition additions, and 1,069,136 exact raw-row replays.  The
raw and distinct occupancy histograms are both retained, including zero-width
buckets; ranks use only the distinct normalized rows.

Ranks are computed modulo `2^61-1`.  This is an exact rank certificate for
both the rational and P-256 coefficient fields, not a probabilistic modular
rank check.  Entries are in `{-2,-1,0,1,2}`, dimensions are at most 15, and
the largest Hadamard determinant bound is `2^40.652`, below both the rank
modulus and the P-256 subgroup order.  Thus nonzero minors cannot disappear
or appear through modular wraparound.

## Typed transport assessment

The transfer skill's main instructions were available, but its referenced
methodology and assessment template could not be retrieved.  The runner emits
the required obligation statuses manually and preserves that limitation.

| path | typed map | subgroup certificate | quotient result | advantage | exploration boundary |
|:--|:--|:--|:--|:--|:--|
| base-field endomorphism | `a+b*pi` on P-256 | `pi(P)=P`; action is scalar `[a+b]` | exact scalar orbit has `K=1` | refuted: generic claw/rho or 1.446505 times rho | complete for `End(E)` acting on `E(F_p)` |
| degree-11 isogeny walk | separable isogeny between order-`n` subgroups | `gcd(n,11)=1`, so kernel intersection is trivial | no quotient reduction in screened bases | refuted in measured family: 494.612 times rho | 4,096 edges certified; full-depth inventory staged, not exhaustive |
| extension/cover/Jacobian | no supplied executable correspondence | unknown | unknown | no credit | open obligation, not a nonexistence result |

An isogeny can transport a DLP while merely changing its coordinates.  It is
not a relation among independent destination factor-base columns.  Round 28's
best model adds 73 columns, retains the universal degree-four S3 support and
does not reduce `K`; pulling membership back through the rational degree-11
map adds map and denominator constraints.

For a future cover or Jacobian claim, the minimum missing evidence is explicit
source and destination models, fields and dimensions, map formulas, induced
subgroup homomorphism, kernel and exceptional-point audit, recovery map, exact
factor-base action, and paired construction/conversion/verification/recovery
costs against rho.

## Gates

| gate | status | evidence |
|:--|:--:|:--|
| four toy domains complete | pass | exact `2^m*C(B,m)` counts and complete histograms |
| exact replay | pass | 1,069,136 / 1,069,136 rows; zero failures |
| exact rank transfer | pass | largest determinant bound `2^40.652 < 2^61-1 < n` |
| P-256 rank-164 event measured | **fail / unset** | no P-256 event; null tail is only a model diagnostic |
| certified non-generic transport | **fail** | scalar action or coordinate-only isogeny in measured families |
| parity at or below rho | **fail** | neither escape supplies the required event or quotient |
| structured residual degree at most 5 | **fail** | Round 296 observes productive degree 6 |
| zero FP/FN on complete checked instances | pass | exhaustive toy buckets and replays agree |
| usable relation below `2^103` | **fail / unset** | no promoted P-256 relation oracle |
| complete collection below `2^120` | **fail / unset** | no rank-164 event or new quotient |
| projected storage below `2^50` | **fail / unset** | no promoted implementation |
| promoted | **no** | parity-critical gates fail |

## Resources and reproduction

The corrected canonical run and independent replay were isolated, uncontended
and pinned to CPU 4 of the AMD EPYC 9V74 host.  The corrected canonical run
used 1.343144 s wall, 1.286448 s user, 0.056041 s system and 135,908 KiB peak
RSS.  The independent run used 1.526782 s wall and emitted a byte-identical
result.  The six-row isolation ledger preserves three earlier successful runs,
one failed correction run, and the two final runs.  The failed run stopped on
the exact `2*distinct=raw` invariant: identity rows had been oriented, but rows
mapping to the toy curve's nonidentity 2-torsion point had not.  Extending the
same orientation rule to every self-negative target resolved the invariant.
The three older runs predate distinct-row accounting; the first additionally
used `bits(n)` rather than `log2(n)` in the displayed random-map tail.  None of
those four superseded outcomes is used for the corrected claims.

- canonical result: 50,776 bytes, SHA-256
  `8718bbe26751ff0191164b0665618226a97556c64c6e2feeef72468a38206965`;
- semantic evidence SHA-256:
  `f5c634bff7357c4e422811b2f2b22ec5bdbc228e6ab0a7acf7d4bd8d37b7a49e`;
- six-run isolation receipt: SHA-256
  `6a4d95377dd61182cabb37b844fa42357a486fea277a050c1a805a2c527c1f88`;
- independent replay: byte-identical, SHA-256
  `8718bbe26751ff0191164b0665618226a97556c64c6e2feeef72468a38206965`.

```bash
cargo test --release --bin p256_parity_escape_screen
cargo clippy --release --bin p256_parity_escape_screen -- -D warnings
cargo build --release --bin p256_parity_escape_screen --bin isolated_bench
target/release/isolated_bench run --wait --cpus 4 \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_parity_escapes_round298_20261007/isolation.jsonl \
  --label p256-parity-escapes-round298-canonical-v3-normalized -- \
  target/release/p256_parity_escape_screen \
  --round297 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_scalar_stabilizer_round297_20261007/scalar-stabilizer-result.json \
  --round296 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_rr_selector_round296_20261007/rr-selector-result.json \
  --round25 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_scalar_orbit_round25_20261006/scalar-orbit-result.json \
  --round28 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_isogenous_dickson_round28_20261006/isogenous-dickson-result.json \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_parity_escapes_round298_20261007/parity-escape-result.json
```
