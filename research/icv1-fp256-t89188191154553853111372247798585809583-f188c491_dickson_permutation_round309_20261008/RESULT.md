# P-256 Dickson half-turn permutation lift, round 309: result

## Verdict

The Dickson half-turn produces a real same-variable two-row construction, but
its product target is far too sparse to help P-256 index calculus.

The complete registered-base rebuild finds an exact closed intersection under
`tau(x)=-x mod p`:

```text
derived factor base   FB1h6255ce9746fe
columns               65,724
two-cycles            32,862
excluded parent rows  65,734
mapping failures      0
```

The full 33,021,374-byte factor base is stored in the repository's
`ecbench.factor_base_dump/v1-wide` format.  Its canonical FB1 SHA-256 is
`6255ce9746fe1869c4ea08c17f12a4a5c11fc283bc3f6f4ae9acd4794339ed5f`.

Unlike Round 308's independent-coordinate lift, the auxiliary row is a
permutation of the **same** unknown factor-base logs.  A planted product hit
on `(Q,G)` therefore gives useful augmented rank two with no new log block.
All P-256 equations replay exactly and the rank is two over both `2^61-1` and
the P-256 subgroup order in both row orders.

The domain cost blocks promotion.  The exact signed distinct-column S17
domain has only `2^240.729666` rows.  Its mean number of hits is
`2^-15.270334` on one ordinary P-256 target but only **`2^-271.270334`** on
the product target `(Q,G)`.  The first arity whose domain reaches `n^2` is
40; its exact balanced 20+20 enumeration materializes
**`2^279.001098` entries**, or **`2^285.358650` bytes** with the stated
82-byte record.  That particular exhaustive baseline is
`2^151.122587` times rho.  A memoryless generic product walk can reduce
storage and the list overhead, but retains an order-`n` rather than
order-`sqrt(n)` time scale.

The derived base cannot inherit the parent's local degree-four measurement,
so structured degree remains unset.  No unplanted full-depth P-256 relation
was attempted.

## Requirement-to-evidence status

| requested direction / gate | status | Round-309 evidence or gap |
|:--|:--:|:--|
| different factor base with global symmetry | constructed | `FB1h6255ce9746fe`, 65,724 columns, 32,862 exact half-turn cycles |
| repository ICV1 and factor-base formats | verified | fixed ICV1 slug; complete `ecbench.factor_base_dump/v1-wide` artifact |
| two rows on one original log vector | verified as planted control | augmented rank two; no auxiliary log block; exact group replay |
| non-labelled actual product-target relation | failed | zero actual relations reported; S17 expected product hits `2^-271.270334` |
| zero false positives / negatives | verified | complete toy and actual-base controls report zero |
| structured degree `<=5` | not established | parent degree four is not transferred to the intersection base |
| complete cost below `2^120` and rho | failed | exact 40-term balanced baseline has list exponent 279.001098; generic memoryless product search remains order `n` |
| cost per usable relation below `2^103` | failed | neither exact balanced nor generic product search approaches the gate |
| peak materialized storage below `2^50` | failed for exact baseline | `2^285.358650` bytes; memoryless storage alone can pass but time cannot |
| complete sparse linear algebra and recovery | unset | no qualifying relation oracle or degree certificate |
| no discarded branch counted exhaustive | verified | actual-base and toy censuses are complete; P-256 image costs are labelled projections |
| unplanted P-256 relation | correctly not attempted | yield, degree, time, storage, and recovery gates fail |

## Exact factor-base construction

The native builder reconstructs all 131,458 columns of
`FB1h2f8621cda105`, checks its FB1 digest and point-set digest, and then
rebuilds it independently.  On the depth-18 torus fibre, the exponent
half-turn multiplies the torus point by `-1`, hence sends trace `x` to `-x`.
The derived base keeps exactly the columns for which both traces lift to
P-256 points.

| quantity | exact value |
|:--|--:|
| parent columns | 131,458 |
| derived columns / signed points | 65,724 / 131,448 |
| two-cycles / fixed points | 32,862 / 0 |
| excluded parent columns | 65,734 |
| closure / bijection / involution failures | **0 / 0 / 0** |
| point and mapped-coordinate failures | **0** |
| derived point-set SHA-256 | `b622786d559db0a2e4e73c21468b77214466c793811f69bf61a6cf5ae0687e5a` |
| permutation SHA-256 | `19d125650bf56bd3c688c961f904ad02277aad6ccf1dec6c37a9d36efa7d425c` |

This is a new factor base, not a relabeling of all parent columns.  Its degree
must be measured independently.

## Complete toy census

The complete `F_7`, width-four control uses `pi=(0 1)(2 3)` and enumerates all
400 projective log profiles and all 2,401 coefficient rows per profile.

| profile class | profiles | product image | bucket occupancy | `(d=2,1)` hits |
|:--|--:|--:|--:|--:|
| `+1` eigenspace | 8 | 7 | 343 | 0 |
| `-1` eigenspace | 8 | 7 | 343 | 0 |
| rank two | 384 | 49 | 49 | 18,816 |

All 960,400 states, 7,683,200 dot-product terms, and 18,816 target-hit rank
checks completed with zero bucket, replay, rank, false-positive, or
false-negative discrepancy.  Restricting to an eigenspace makes the product
image cheap only by collapsing its useful coordinate rank to one.

## P-256 information certificate

The 17-column deterministic control solves two coefficients so that

```text
sum_i c_i P_i       = Q = [d]G
sum_i c_i P_pi(i)   = G.
```

The resulting augmented rows `(c,-1,0)` and `(pi(c),0,-1)` have rank two over
both coefficient fields and in both row orders.  They use one common factor-
base log vector and introduce no auxiliary unknowns.  The two exact group
replays used 18 constructed points, 52 scalar multiplications, 13,282 scalar
bits, and 34 point additions, with zero failures.  This is a scalar-labelled
information control, not a constructed relation on `FB1h6255ce9746fe`.

## Domain and cost projection

| quantity | exact or projected value | interpretation |
|:--|--:|:--|
| signed S17 domain | `2^240.729666` | exact `2^17*C(65724,17)` |
| S17 ordinary-target mean | `2^-15.270334` | about one hit per 39,521 targets |
| S17 product-target mean | `2^-271.270334` | effectively no `(Q,G)` hit |
| first ordinary-cover arity | 19 | exact domain threshold `>=n` |
| 19-term balanced-list width | `2^148.249278` | exact 9+10 baseline |
| first product-cover arity | 40 | exact domain threshold `>=n^2` |
| 40-term balanced-list width | `2^279.001098` | exact 20+20 baseline |
| 82-byte list storage | `2^285.358650` bytes | fails materialization gate |
| one write plus one read | `2^286.358650` bytes | disk-traffic projection |
| balanced-baseline ratio to rho | `2^151.122587` | projection, not measurement |

The `2^279.001098` figure prices the exact balanced distinct-column baseline;
it is not asserted as a universal lower bound.  Generic memoryless collision
methods can trade the list away, but a uniform product target still has a
256-bit rather than 128-bit search exponent.  Either accounting fails the
promotion threshold by a wide margin.

## Decision and next obligation

Class: **negative with an information-structure advance**.  The half-turn
permutation is the first screened product lift here that gives two useful rows
without adding unknown logs.  It does not improve the executable boundary
because its joint target space and relation domain are fatal.

The next viable target is now precise:

> Find a same-base non-scalar correspondence whose two useful rows lie on a
> searchable image of size about `n`, not `n^2`, without collapsing to one
> eigendirection, and retain structured degree at most five.

Any such claim must explain why the smaller joint image is not an affine-log
transport or a hidden direct-DLP witness.

## Isolation and artifact integrity

The protocol was committed as `7dabcdddc`, the implementation as `425e9d9dc`,
and a post-run lint-only divisibility expression was corrected in `e8e797b42`.
The v1 receipts are preserved; canonical-v2 and independent-v2 correspond to
the final committed source.  All four runs were uncontended, and every run
emitted byte-identical factor-base, result, and assessment files.

| run | exit | wall | user | system | peak RSS | contention |
|:--|--:|--:|--:|--:|--:|:--|
| canonical-v1 | 0 | 112.473719 s | 112.351549 s | 0.108168 s | 303,780 KiB | none |
| independent-v1 | 0 | 112.545868 s | 112.407656 s | 0.134469 s | 303,760 KiB | none |
| canonical-v2 | 0 | 115.407440 s | 115.275865 s | 0.125780 s | 303,908 KiB | none |
| independent-v2 | 0 | 111.193284 s | 111.035712 s | 0.153617 s | 303,784 KiB | none |

| artifact | bytes | SHA-256 |
|:--|--:|:--|
| `half-turn-factor-base.json` | 33,021,374 | `677f5ddbe0f2ea1e0c2bd29272d76461acc0d0197d11843d3ebf77d91eaf15ac` |
| `dickson-permutation-result.json` | 8,204 | `1fa17da4d11345f2dce55d6704355d193e081eae6912d3b45a9605e2efa32bc1` |
| `transfer-assessment.json` | 3,465 | `f756692af527fa73d6875a8131b71c163717685d7ff423db0506aed3344ec3d8` |
| `isolation.jsonl` | 11,664 | `4a11e0bee6776087a42c97f04fd1a95eaf9907ae1a5a5f330d7b09bc0c5ad9b0` |
| `PROTOCOL.md` | 7,819 | `a425c826de3b519fb88d405d8dd6ab052d8a8ce867a494890cb0d19dd28ef795` |

Result semantic evidence:
`e53af26ce2bef69d7ab98b4e30e7f90966e3efa9cb678f0e3a880367f3431159`.
Transfer-assessment semantic evidence:
`3a25fab3bfd548dee4af0a16fdcb0f95747b17a1ac8bdf197d2267f3461591b9`.

The transfer workflow's methodology and assessment-template companion
resources were unavailable, and that limitation is recorded in the assessment.

## Reproduction

```bash
cargo test --bin p256_dickson_permutation
cargo clippy --bin p256_dickson_permutation -- -D warnings
cargo build --release --bin p256_dickson_permutation --bin isolated_bench
target/release/isolated_bench run --wait --cpus 4 \
  --label p256-dickson-permutation-round309-canonical-v2 \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_dickson_permutation_round309_20261008/isolation.jsonl -- \
  target/release/p256_dickson_permutation \
  --round2 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_dickson_terminal_degree_round2_20261004/factor-base-result.json \
  --round22 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_factor_base_symmetry_round22_20261006/symmetry-result.json \
  --round39 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_affine_log_orbit_round39_20261006/affine-log-orbit-result.json \
  --round308 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_product_lift_round308_20261008/product-lift-result.json \
  --factor-base-out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_dickson_permutation_round309_20261008/half-turn-factor-base.json \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_dickson_permutation_round309_20261008/dickson-permutation-result.json \
  --assessment research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_dickson_permutation_round309_20261008/transfer-assessment.json
```
