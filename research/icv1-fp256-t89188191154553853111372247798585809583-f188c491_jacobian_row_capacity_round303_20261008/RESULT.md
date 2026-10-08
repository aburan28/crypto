# P-256 Jacobian row capacity, round 303: result

## Verdict

A homomorphic same-target cover/Jacobian route cannot meet the one-event
rank-164 parity gate below **genus 165**.

For a separable cover `pi:C -> E`, differences between points in one target
fibre map into `Prym(pi)`, whose dimension is `g-1`.  If that Prym contains
`m` copies of the P-256 elliptic isogeny factor, then `m<=g-1`.  Round 25
certifies that rational endomorphisms of each P-256 factor act as known scalars
on `E(F_p)[n]`, so one factor supplies at most one independent order-`n` log
coordinate after scalar transports are quotiented.  Conditionally on the rows
factoring through these rational homomorphisms,

```text
same-target homogeneous row rank <= m <= dim Prym(pi) = g-1.
```

The fixed factor base has 164 independent columns.  Round 298 shows that one
collision event is `0.9640878` times rho before omitted costs, but two events
already cost `1.4461317` times rho.  Therefore the first homomorphic capacity
that can even pass the event-count boundary requires all of:

```text
genus(C)                         >= 165
dim Prym(pi)                     >= 164
P-256 multiplicity in Prym(pi)  >= 164
distinct normalized fibre points >= 165
degree(pi)                       >= 165
```

At the minimum genus, Riemann--Hurwitz requires total ramification degree 328.
An unramified high-degree elliptic isogeny stays genus one and has zero Prym
dimension, so degree alone is not an escape.

The 256-sheet abstract character control reaches rank 164 over both
`2^61-1` and the P-256 subgroup order.  This validates the capacity accounting;
it does **not** construct a genus-165 curve, a rational fibre, a relation
system, or an attack.  Every promotion gate depending on such a construction
therefore fails, structured degree remains unset, and no unplanted P-256
relation was attempted.

## Requirement-to-evidence status

| requested direction / gate | status | Round-303 evidence or gap |
|:--|:--:|:--|
| pursue many independent rows per event | bounded for homomorphic cover rows | same-target rank is at most the P-256 multiplicity in the Prym, hence at most `g-1` |
| exact global transport accounting | verified | Round-25 scalar action is hash pinned; one P-256 factor contributes at most one quotient coordinate |
| exact character/rank controls | verified | 11 complete matrices, two moduli, forward and reversed sheet order, zero discrepancies |
| one-event rank 164 | capacity only at genus 165 | 256-sheet abstract matrix has rank 164; no curve realizes it |
| explicit qualifying P-256 cover | not established | genus `>=165`, degree `>=165`, Prym multiplicity `>=164` remain construction obligations |
| zero false positives / negatives on checked controls | verified | all 44 rank computations equal their declared character count |
| structured residual degree `<=5` | not established | no relation system exists for the capacity-only model |
| complete collection below `2^120` | failed / unset | no cover, rational fibre yield, solver, or recovery pipeline |
| usable relation below `2^103` | failed / unset | no P-256 relation oracle |
| storage below `2^50` | failed / unset | the 1.28 MiB control is not the projected attack storage |
| complete cost at or below rho | failed | the only below-rho row is an unconstructed capacity threshold |
| actual unplanted P-256 relation | correctly not attempted | construction, degree, cost, and storage gates fail |

## Boundary table

One unit is the optimistic collision-event ratio to Pollard rho.  Every cover
row is a capacity lower boundary with a free oracle and omits construction,
fibre search, algebra, recovery, and linear algebra.

| variant | maximum homogeneous rank / event | events | optimistic ratio to rho | status |
|:--|--:|--:|--:|:--|
| Pollard rho | n/a | n/a | **1.000000** | complete reference |
| genus 2 homomorphic cover | 1 | 164 | **13.920747** | a single Prym direction does not reduce the event count |
| genus 83 capacity | 82 | 2 | **1.446132** | still above rho |
| genus 164 capacity | 163 | 2 | **1.446132** | one row short; still above rho |
| genus 165 capacity | 164 | 1 | **0.964088** | capacity only; no cover or attack |
| registered 17-term `FB1h2f8621cda105` | 1 | 164 | **394.425280** | unchanged free-oracle comparison |

Class: **negative**.  No executable boundary moved.

## Exact capacity profiles

Each row completely enumerates the smallest power-of-two sheet set supporting
the declared nontrivial Walsh characters.  The sheet-zero row is subtracted,
and the homogeneous matrix is ranked independently over both coefficient
fields and again in reversed sheet order.

| genus | Prym/rank cap | events | ratio to rho | sheets | rank over both fields | P-256 matrix bytes |
|--:|--:|--:|--:|--:|--:|--:|
| 2 | 1 | 164 | 13.920747 | 2 | 1 | 32 |
| 3 | 2 | 82 | 9.835955 | 4 | 2 | 192 |
| 5 | 4 | 41 | 6.944477 | 8 | 4 | 896 |
| 9 | 8 | 21 | 4.955602 | 16 | 8 | 3,840 |
| 17 | 16 | 11 | 3.567258 | 32 | 16 | 15,872 |
| 33 | 32 | 6 | 2.609816 | 64 | 32 | 64,512 |
| 65 | 64 | 3 | 1.807665 | 128 | 64 | 260,096 |
| 83 | 82 | 2 | 1.446132 | 128 | 82 | 333,248 |
| 129 | 128 | 2 | 1.446132 | 256 | 128 | 1,044,480 |
| 164 | 163 | 2 | 1.446132 | 256 | 163 | 1,330,080 |
| 165 | 164 | 1 | 0.964088 | 256 | 164 | 1,338,240 |

All four ranks in every row—two fields times two sheet orders—equal the
declared cap.  The runner materialized 44 matrices, 4,556 rows, and 548,936
entries, then charged 2,656 pivot scans, 2,656 inversions, 6,134,196 modular
multiplications, and 5,968,800 modular subtractions.  No pivot swaps were
needed under either deterministic ordering.

## Registered split-cover control

For the Round-302 genus-two cover,

```text
J(H) ~ E x E'
```

and the Prym is the complementary elliptic curve `E'`.  The exact `[n]R != O`
certificate proves that `E'(F_p)` has no rational P-256 order-`n` factor.
Thus its P-256 Prym multiplicity and same-target sheet-rank gain are both zero,
matching Round 300's quotient-rank gain of zero.

A hypothetical genus-two self-gluing with Prym isogenous to P-256 would raise
the ceiling only to one homogeneous row per fibre event.  It would still need
164 events and retain the `13.920747`-times-rho free-oracle boundary.  This is
why another degree-two cover is not the target; parity requires 164 Prym
directions simultaneously.

## Scope and open obligation

The result applies to rows obtained from same-target fibre differences through
rational homomorphisms from a smooth cover Jacobian to P-256 factors.  It is
not a universal lower bound on encodings or relation mechanisms that do not
factor through such homomorphisms.

The weakest surviving obligation is consequently precise: construct and
replay a genus-at-least-165, degree-at-least-165 P-256 cover whose Prym contains
at least 164 rational P-256 order-`n` factors, or supply a genuinely
non-homomorphic mechanism outside this dimension bound.  Either route must
still meet degree-at-most-five, full relation/recovery cost, storage, and exact
replay gates.

## Isolation and artifact integrity

The canonical run and final independent replay were isolated and uncontended.
The first independent replay was correct and byte-identical but marked
contended because another process consumed 0.07 CPU seconds; its receipt is
preserved and not used as the final isolation claim.

| run | exit | wall | user | system | peak RSS | contention |
|:--|--:|--:|--:|--:|--:|:--|
| canonical-v1 | 0 | 0.379063 s | 0.378956 s | 0 s | 4,084 KiB | none |
| independent-v1 | 0 | 0.449038 s | 0.448922 s | 0 s | 4,072 KiB | yes; preserved |
| independent-v2 | 0 | 0.478966 s | 0.476257 s | 0.002604 s | 4,072 KiB | none |

The canonical and both independent invocations emitted byte-identical result
and assessment JSON.

| artifact | bytes | SHA-256 |
|:--|--:|:--|
| `jacobian-row-capacity-result.json` | 18,217 | `6b023f47fbb2c57eb79f870f1c25be48e0b4230d974c3823caf76784a095caa0` |
| `transfer-assessment.json` | 4,303 | `17c13bd7309ae8399c1ccd83b83887d68ce95d15f9e1b9fd1ced1bd466f7c62f` |
| `isolation.jsonl` | 7,409 | `6b5d3e05c7c74e4191f41882aafe33fdd18a9541aead71f3719adb7afc94bc50` |
| `PROTOCOL.md` | 7,461 | `19ba69e7a7df73e4618c3aceefc4b4f8126595ad8e7bf40f8b41c02527dd6e5f` |

The result semantic-evidence SHA-256 is
`f3dfdf11427c3c0b4c7a9ddd913e8da2f28a4264b5659ade3768c2f5de995b6b`.
The transfer assessment semantic hash is
`a11f0c269c2c98c5ea6a617b3f760cf5df939098627079f3e236dc7b48437e46`.

The transfer workflow's referenced `references/methodology.md` and
`assets/assessment-template.json` resources remained unavailable.  The typed
graph, obligations, controls, scope, and weakest open obligation are emitted
directly and the limitation is recorded in the assessment JSON.

## Reproduction

```bash
cargo test --bin p256_jacobian_row_capacity
cargo clippy --bin p256_jacobian_row_capacity -- -D warnings
cargo build --release --bin p256_jacobian_row_capacity --bin isolated_bench
target/release/isolated_bench run --wait --cpus 4 \
  --label p256-jacobian-row-capacity-round303-canonical-v1 \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_jacobian_row_capacity_round303_20261008/isolation.jsonl -- \
  target/release/p256_jacobian_row_capacity \
  --round25 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_scalar_orbit_round25_20261006/scalar-orbit-result.json \
  --round298 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_parity_escapes_round298_20261007/parity-escape-result.json \
  --round300 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_cover_fiber_round300_20261007/cover-fiber-result.json \
  --round302 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_complementary_quotient_round302_20261008/complementary-quotient-result.json \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_jacobian_row_capacity_round303_20261008/jacobian-row-capacity-result.json \
  --assessment research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_jacobian_row_capacity_round303_20261008/transfer-assessment.json
```
