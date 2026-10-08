# P-256 product-lift information and collision cost, round 308: result

## Verdict

A direct two-coordinate product lift does **not** turn one collision into two
useful equations for the original P-256 recovery system.

For

```text
Phi(c) = (sum_i c_i P_i, sum_i c_i R_i) in G x G,
```

a full collision does replay two coordinate equalities.  The information
accounting has three distinct outcomes:

1. if `R_i=[a]P_i`, the second equality is a known scalar multiple and the
   useful rank is one;
2. if the auxiliary logs are independent and unknown, the expanded block
   rank is two but a new log block has been introduced, and projection onto
   the original unknowns has rank one; and
3. if the auxiliary labels are known, their equality is a coefficient filter,
   not a new equation in the unknown P-256 logs.

The complete order-five control checked 7,447,750 coefficient pairs and
279,000 full product collisions with zero false positives, false negatives,
or replay failures.  The deterministic 17-column P-256 control replayed all
four linked/independent coordinate equations.  Ranks agree over `2^61-1` and
the P-256 subgroup order in both row orders: linked rank one, independent
expanded rank two, and original-system projection rank one.

For comparison, a uniform balanced full collision in `G x G`, even under the
optimistic simultaneous-negation quotient, costs about **2^255.825748**
samples.  That is **2^127.947236 times Pollard rho**.  Crediting the event with
its one useful original-system row gives an optimistic lower bound of
**2^272.900381** operations for 138,031 rows, before duplicate loss, actual
relation yield, sparse linear algebra, or recovery.  A memoryless walk can
avoid the materialized-list storage, but it cannot remove the order-`n` time
exponent.

The direct product route therefore fails the time, useful-rank, degree, and
complete-attack gates.  The executable lower boundary remains **13.920747
times rho** for the minimum-width control and **394.425280 times rho** for
`FB1h2f8621cda105`.  No unplanted full-depth P-256 relation was attempted.

## Requirement-to-evidence status

| requested direction / gate | status | Round-308 evidence or gap |
|:--|:--:|:--|
| pursue many independent rows per event | refuted for direct product lifts | linked rank is one; independent block rank two projects to original-system rank one |
| exact group replay | verified | four P-256 coordinate equations and 558,000 toy collision equations replay; zero failures |
| zero false positives / negatives | verified | complete 7,447,750-pair toy census reports zero of each |
| two useful non-labelled original-system rows | failed | no screened case supplies a second row on the original log variables |
| structured residual degree `<=5` | not established | no promoted coupled relation system exists; degree remains unset |
| complete collection below `2^120` | failed | optimistic 138,031-row product-event lower bound is `2^272.900381` before omitted work |
| usable relation below `2^103` | failed | optimistic product event is `2^255.825748` samples for one useful original row |
| peak materialized storage below `2^50` | split / complete gate failed | 82-byte birthday list is `2^262.183300` bytes; 512-byte memoryless state passes storage alone but fails time |
| parity with rho | failed | generic balanced product event is `2^127.947236` times rho |
| 138,031 rows and sparse linear algebra | lower bound only | 138,031 includes the registered 5% rank allowance; duplicates, yield, and sparse LA are omitted, so complete cost is higher |
| no discarded branch counted exhaustive | verified | the complete toy census is exhaustive; the P-256 cost is explicitly a model, not measurement |
| actual unplanted P-256 relation | correctly not attempted | useful-rank, degree, time, and complete-attack gates fail |

## Exact complete toy census

The preregistered cell is the cyclic group of order five with three columns.
All 31 normalized nonzero projective profiles were paired in all 961 ordered
ways.  Every one of the 125 coefficient states and every unordered state pair
was checked.

| quantity | exact value |
|:--|--:|
| projective profile pairs | 961 |
| proportional / independent pairs | 31 / 930 |
| coefficient states evaluated | 120,125 |
| unordered coefficient-pair tests | 7,447,750 |
| first-coordinate collisions | 1,441,500 |
| full product collisions | 279,000 |
| proportional / independent full collisions | 46,500 / 232,500 |
| independent first-only collisions | 1,162,500 |
| replayed collision equations | 558,000 |
| coordinate / bucket / replay discrepancies | **0 / 0 / 0** |
| false positives / false negatives | **0 / 0** |

The linked representative has rank one in both row orders.  The independent
representative has rank two in its six-column block system but rank one after
projection to the original three columns.  The runner performed 720,750
dot-product terms and 1,441,500 explicit cyclic additions; logical peak
materialization was 625 bytes.

## Deterministic P-256 certificate

The runner hash-selected 17 nonzero primary labels and 17 nonzero auxiliary
labels modulo the exact P-256 subgroup order, then materialized all 51 points:
17 primary, 17 linked auxiliary, and 17 independent auxiliary points.  All
are valid public points.

| system | rows | columns | rank over both fields and orders | useful original-system rank | status |
|:--|--:|--:|--:|--:|:--|
| linked `R_i=[a]P_i` | 2 | 17 | **1** | **1** | second replay is a known scalar multiple |
| independent expanded blocks | 2 | 34 | **2** | **1** | rank gain comes with a new unknown log block |
| independent, original projection | 2 | 17 | **1** | **1** | auxiliary row projects to zero |
| known-label selector | n/a | 17 | n/a | **0 new rows** | condition is computable from coefficients |

The four P-256 group replays used 61 scalar multiplications, 15,552 scalar
bits, 7,890 nonzero scalar bits, and 10 point additions.  Twelve rank matrices
used 544 materialized entries, 16 inversions, 408 modular multiplications, and
68 modular subtractions.  The witnesses are scalar-labelled controls and
receive no construction credit on either actual factor base.

## Product-collision projection

The exact simultaneous-negation state counts are

```text
original: (n+1)/2
  57896044605178124381348723474703786764998477612067880171211129530534256022185

product:  (n^2+1)/2
  6703903961849550000561278353995505841769572391168607857098878760322432797127854223500633464330763947064984370529374005771557885306222791227620809912304081
```

The optimistic birthday estimate is `sqrt(pi)/2*n`, or `2^255.825748`
samples.  Relative to `1.3*sqrt(n/2)` rho it is
`3.2806207776272584e38 = 2^127.947236` times rho.

| projection | log2 cost / size | gate status |
|:--|--:|:--|
| one product event, credited one useful original row | 255.825748 operations | fails per-row `<103` |
| 138,031 useful rows with 5% rank allowance | 272.900381 operations | fails complete `<120`; optimistic lower bound only |
| materialized list, 82-byte records | 262.183300 bytes | fails storage `<50` |
| one write plus one read of that list | 263.183300 bytes traffic | fails |
| memoryless product walk | 512-byte state | storage passes alone; time unchanged |

Conditioning directly on a known auxiliary linear constraint may avoid
materializing the full product list, but it still supplies zero new equations
in the unknown primary logs.  Conversely, retaining an independent unknown
auxiliary coordinate supplies a second block equation only by adding the
corresponding unknown block.  Neither interpretation earns two-row credit.

## Boundary table

| variant | useful original rows / event | ratio to rho | status |
|:--|--:|--:|:--|
| Pollard rho | n/a | **1.000000** | complete reference |
| free two-primitive oracle | 2 | **1.446505** | already above rho; unconstructed |
| minimum-width `FB1hc72514a2a8d3` | 1 | **13.920747** | executable free-oracle boundary unchanged |
| registered 17-term `FB1h2f8621cda105` | 1 | **394.425280** | unchanged free-perfect-oracle comparison |
| generic balanced product lift | 1 | **3.280621e38** | order-`n` model; fails |

Class: **negative**.  No executable boundary or structured degree improved.

## Scope and next obligation

This result is exact for direct product lifts whose auxiliary coordinate is a
known scalar transport, an independent unknown cyclic-group coordinate, or a
known-label selector.  It is not a universal theorem about nonlinear or
non-homomorphic correspondences.

The weakest remaining obligation is now narrower:

> Construct a typed non-homomorphic event whose two exactly replayed equations
> both survive projection to the original P-256 recovery variables, introduce
> no unresolved auxiliary log block, and cost no more than original-group rho.

Any candidate must then provide its actual relation system, structured degree,
yield, 138,031-row collection cost, sparse linear algebra, recovery, and exact
group replay.  Local degree alone is not promotion evidence.

## Isolation and artifact integrity

The preregistration was committed as `991d0fe86850f44ebee34c995cfa2b98f139f066`
before the runner was committed as
`d5b5065bbc91575be75600415f377313dc39a732`.  Canonical and independent
isolated runs were uncontended and emitted byte-identical result and assessment
JSON.

| run | exit | wall | user | system | peak RSS | contention |
|:--|--:|--:|--:|--:|--:|:--|
| canonical-v1 | 0 | 0.333924 s | 0.329731 s | 0.004053 s | 2,256 KiB | none |
| independent-v1 | 0 | 0.341273 s | 0.341162 s | 0 s | 2,256 KiB | none |

Artifact SHA-256 values:

| artifact | bytes | SHA-256 |
|:--|--:|:--|
| `product-lift-result.json` | 12,245 | `21e8684f780835a822208061a962c843b097aa810b889bc2a09e2dcb85b88ed4` |
| `transfer-assessment.json` | 4,620 | `4a11fd61cd3c816e85eb9d27ed609e5d87fcc92810a364d10cfdba8ca4a62750` |
| `isolation.jsonl` | 5,162 | `dfc90bbd52a740c26ae0ce814f4e15b052a681c701b4ccd824fd3147f62efcc2` |
| `PROTOCOL.md` | 9,315 | `309ddbc53d8ad212adad296bc8baecc848e9a90edf28623ae693cbaf32ef9b38` |

The result semantic-evidence SHA-256 is
`6320f53c8719793b983cfc21f0e123e7c100f14ff43a2711f3f32edbef4417cd`.
The transfer assessment semantic hash is
`09b396da4801ef808aad3837b0f024e2a3ad52eae8f484f13a64a25789370686`.

The transfer workflow was used for the typed-map and obligation audit.  Its
referenced `references/methodology.md` and
`assets/assessment-template.json` resources were unavailable in this
environment, so the limitation is recorded in the assessment JSON.

## Verification and reproduction

The scoped test binary completed three tests with zero failures.  The normal
repository commands are:

```bash
cargo test --bin p256_product_lift
cargo clippy --bin p256_product_lift -- -D warnings
cargo build --release --bin p256_product_lift --bin isolated_bench
target/release/isolated_bench run --wait --cpus 4 \
  --label p256-product-lift-round308-canonical-v1 \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_product_lift_round308_20261008/isolation.jsonl -- \
  target/release/p256_product_lift \
  --round39 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_affine_log_orbit_round39_20261006/affine-log-orbit-result.json \
  --round302 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_complementary_quotient_round302_20261008/complementary-quotient-result.json \
  --round303 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_jacobian_row_capacity_round303_20261008/jacobian-row-capacity-result.json \
  --round307 research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_event_provenance_round307_20261008/event-provenance-result.json \
  --out research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_product_lift_round308_20261008/product-lift-result.json \
  --assessment research/icv1-fp256-t89188191154553853111372247798585809583-f188c491_product_lift_round308_20261008/transfer-assessment.json
```

The local workspace's shared Cargo target was read-only and the remaining
quota could not hold another full dependency rebuild.  The recorded runs used
Rust 1.98.1 to compile the committed runner directly against the repository's
existing release dependency artifacts; this execution-environment detail does
not change the source, inputs, or deterministic outputs.
