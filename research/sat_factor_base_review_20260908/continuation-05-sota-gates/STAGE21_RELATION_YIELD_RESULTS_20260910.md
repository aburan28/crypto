# Stage 21: production relation-yield result

The frozen production run completed at source commit
`cf1a564757f9d39e21d8db8a211fc77b027c3c37`. It measured the exact
two-summand relation yield of one target-independent algebraic factor base on
public synthetic Koblitz targets. The terminal run and its outer process receipt
are sealed, and the control-plane verifier replayed the derivable payload
relationships and hashes.

The result is still a **candidate pending independent payload replay**. The
verifier did not independently rematerialize the factor base, replay all
9,165,621 curve additions, or re-add the retained point witnesses. It therefore
sets `scientific_measurement_admitted: false` and
`measurement_admission_status: pending_independent_payload_replay`.

## Measured yield and controls

| Arm | Targets | Exact pair-table hits | Misses | Interpretation |
|---|---:|---:|---:|---|
| Natural hash-to-subgroup | 256 | **163** | 93 | Primary yield estimate |
| Planted SAT | 64 | **64** | 0 | Witness-control pass |
| Exact-table proven UNSAT | 64 | **0** | 64 | Exact-miss-control pass |

The natural two-summand hit rate is

\[
\hat p = \frac{163}{256}=0.63671875.
\]

The two-sided 95% Wilson score interval is
`[0.576185047549112, 0.6932099917094386]`. The observed trials per
relation are `256/163 = 1.5705521472392638`.

Natural targets come from the frozen domain-separated BLAKE3 stream, exact
uniform-affine decoding, and public cofactor projection. No target scalar or
factor-base discrete-log label is constructed. Accepted targets are distinct,
so the sample is technically without replacement. The 256 targets are about
`0.000122146` of the 2,095,852-point nonidentity subgroup population; no
finite-population correction is applied. The Wilson interval treats the frozen
BLAKE3 stream as a pseudorandom Bernoulli sample from this stated distribution.

The planted arm is selected by priority over canonical factor-base pairs.
Targets with several decompositions receive several priority opportunities, so
this arm is multiplicity-weighted. Its 64/64 result validates construction,
lookup, subgroup membership, and witness handling; it is not a natural
prevalence estimate. The proven-UNSAT arm establishes absence only for two
summands over this exact materialized factor base, with repeated points allowed.
It is not higher-summand or SAT-solver UNSAT evidence.

## Public factor-base decision

The production run repeated public discovery on both curve variants before the
yield process. Selection used no target, relation yield, solver timing, subgroup
enumeration, or discrete-log labels.

| Curve | Selected divisor | Rational points | Projected points | Signed orbits before projection | Projected signed columns | Selected for yield |
|---|---|---:|---:|---:|---:|:---:|
| `K_0 / F_(2^23)` | `[0,2]` | 4,281 | 4,279 | 95 | 93 | yes |
| `K_1 / F_(2^23)` | `[0,1]` | 3,957 | 3,957 | 87 | 86 | no |

The selected `K_0` base has dimension 12, 4,096 abscissae, divisor
polynomial bitmask `5279`, group order 8,383,412, prime subgroup order
2,095,853, and cofactor 4. Its factor-base hash is
`9b5bb635f1c505cb871a9e8723d2b79109688606cc42323ab92e05b8beaf4b58`;
the predicate hash is
`563f00a4ac14ee2714fd0b27b64d2fb2beab8d3d9e7edbd5d13384516b2f38c0`.

## Exact reference oracle

The producer enumerated every canonical pair `i <= j` exactly once:

| Quantity | Value |
|---|---:|
| Factor-base points | 4,281 |
| Required canonical pairs `F(F+1)/2` | 9,165,621 |
| Enumerated pairs | 9,165,621 |
| Unique point sums | 5,575,848 |
| Duplicate pair sums | 3,589,773 |
| Hash-table capacity | 14,680,064 |

The canonical pair transcript hash is
`c7a8ab74fc4ed9ea3832244f99b92e87b75ea268ffbcc89fb72fc15cd069ccae`.
The terminal result binding is
`a06e2e096049b43d9faeb9f42df54fe56d3c346c63032b558db34dca85c80ab0`.

## Internal stage timing

These are monotonic in-process wall timers. Nested fields are not additive: in
particular, planted selection occurs inside the canonical pair-table scan.

| Stage | Seconds |
|---|---:|
| Curve and subgroup construction | 0.001199291 |
| Factor-base predicate and materialization | 0.067709708 |
| Projected-column census | 0.057080833 |
| Natural target generation | 0.010826209 |
| Cofactor-class construction | 1.664421041 |
| Canonical pair table plus planted-priority selection | **68.818071959** |
| Planted construction-witness and subgroup verification | 0.019228000 |
| Proven-UNSAT target screening | 0.010400083 |
| Covariate extraction | 0.004093917 |
| Natural exact queries / witness verification | 0.000046417 / 0.001247458 |
| Planted exact queries / witness verification | 0.000006083 / 0.000455708 |
| Proven-UNSAT exact queries / witness verification | 0.000006958 / 0.000001458 |
| **Yield producer end to end** | **70.661434458** |

The public-discovery internal end-to-end timers were 0.502843041 seconds for
`a=0` and 0.277381875 seconds for `a=1`.

## Fresh-process resource accounting

| Process | Core-s | Process wall-s | Peak RSS bytes |
|---|---:|---:|---:|
| Fresh locked release build | 107.348907 | 118.334764625 | 1,054,294,016 |
| Public discovery, `a=0` | 0.432730 | 0.931791917 | 6,242,304 |
| Public discovery, `a=1` | 0.259528 | 0.334594459 | 6,012,928 |
| Relation-yield producer | 61.992819 | 70.899883458 | 293,339,136 |
| **Four charged children** | **170.033984** | **190.501034459 summed** | **1,054,294,016 maximum** |

The inclusive whole-driver receipt records 163.834223 user-seconds,
7.273478 system-seconds, **171.107701 core-seconds**, **193.640836375 seconds
elapsed wall**, and **1,054,294,016 bytes peak process RSS**. It encloses the
child processes and must not be added to the child sum. The outer-minus-child
difference is 1.073717 core-seconds and 3.139801916 wall-seconds.

`single_core_seconds` in the process receipts is a legacy alias of total user
plus system CPU. Measured single-core elapsed time is unavailable. Peak RSS is
the maximum individual process high-water mark, not a sum. Scientific producers
were sequential and requested one worker, but the Cargo parent and compiler
child may overlap; simultaneous aggregate build-tree RSS is unavailable.
Conflicts are inapplicable because the exact pair table is not a
conflict-driven SAT solver.

## Operation counts and retained-size bounds

The producer reports these high-level operations:

- 9,165,621 pair-table group additions;
- 4,281 factor-base projection scalar-multiplication calls;
- 4,281 cofactor-class scalar-multiplication calls and 4,281 class negations;
- 2,296,723 planted pair-priority hashes;
- 64 final planted subgroup checks and 64 construction-witness re-additions;
- 8,832 Frobenius applications for target covariates;
- natural generation: 522 hash candidates, 266 affine-decode rejections,
  256 cofactor projections, zero oracle lookups during selection, 256 final
  lookups, and 163 witness re-additions;
- proven-UNSAT generation: 462 hash candidates, 241 affine-decode rejections,
  221 cofactor projections and selection lookups, 157 rejected pair-table hits,
  64 retained exact misses, and 64 final lookups.

Scalar-multiplication internals and factor-base-constructor internals are
covered by process time but are not relabeled as counted group additions.

| Retained payload lower bound | Bytes |
|---|---:|
| Factor-base payload | 42,291 |
| Canonical factor-point clone | 29,967 |
| Live pair-table key/witness payload | 89,213,568 |
| Pair-table capacity key/witness payload | 234,881,024 |
| Selected target coordinate payload | 2,688 |
| **Sum using table-capacity payload** | **234,955,970** |

These are encoded or flat key/value payload lower bounds. They exclude Rust
object headers, allocator metadata, hash-table control bytes, and unrelated
capacity. They are not heap measurements and do not replace the external RSS
receipt.

## Custody and replay

- Frozen protocol SHA-256:
  `6ab87d1c2b98e36858b506fbb62c43be920b0e4a43150c860024d806e379ca25`
- Committed production tree SHA-256:
  `3e6a5f41e51f2d7dbee2ad030b9d019a3b6c05974db267f4221e2165b03a263d`
- Source binding SHA-256:
  `9ed7f379515117928cabbc085af225b98fcedb56b8aec2a33b354c856eb8c253`
- Host binding SHA-256:
  `87030247ea5a1481c58b3fb1785015fa0ecb3ca7f0707fbe66b68def5119711d`
- Yield binary SHA-256:
  `d65509158c9f80becdc355b1e286a4760c08072cf2c49b278b207763a0e022dc`
- Discovery binary SHA-256:
  `189368f194209283e81e6a7c35c6305322c76dfd68f39540782c7a8370f1b150`
- Yield result SHA-256:
  `d69165c4cf6207f03611a08d62ee66dd5c600ae6926ff9d6e21178ba50487f42`
- Run summary SHA-256:
  `83aca6d143ee90e977582044511236032bc2259e8afdc5608cd561487ac1adce`
- Raw run-seal SHA-256:
  `5bc7d97265832cdb0e2ed6eb3e085fa1009a8d54e9c4816def4337988372f51b`
- Canonical run-inventory SHA-256:
  `66e56bbd4b7f35486e87b8af479e2dd320b84f361c3769303e3f71c576d7fb21`
- Outer metrics SHA-256:
  `bf7f582fe9c3a86da4717a65fe9f26bae080581a6d95adf2937c05b9b6260402`
- Verification SHA-256:
  `8b542621d6bb02bfeb7529c1ec6f490b0d4e9efbc9d245319a6ad976767abfb3`
- Raw verification-seal SHA-256:
  `968339f5e57f945ab0fd2bc700b57b5e5dfd6a80cc4539a8fe1b039a8211d342`

The verifier recomputed the policy, predicate, ordered-covariate-row, arm
target, arm witness, and terminal result-binding BLAKE3 commitments. It did not
replay:

1. factor-base point materialization and the ordered factor-base hash;
2. all canonical pair additions and the canonical-pair transcript hash;
3. exact curve re-addition of retained witness indices.

Those three items are the remaining independent payload replay required before
scientific admission. The verification directory has no separate process-meter
receipt, so verifier CPU, wall, and RSS are not included in the production
resource totals.

## Later CI-only changes

The production packet is bound to `cf1a5647`. Four later commits through
current head `dfa36f4d865301afec1c654de90b0e5a00640e2b` improve only CI failure
visibility and Rust toolchain-shim dispatch. Their diff changes the workflow,
custody runner, and runner tests. There is no diff after `cf1a5647` in
`examples/koblitz_relation_yield_bridge.rs` or `src/**`; the producer and its
cryptographic mathematics are unchanged.

## Admission and claim boundary

The following states remain authoritative:

- `scientific_measurement_admitted: false`;
- `measurement_admission_status: pending_independent_payload_replay`;
- `external_portable_verification_satisfied: false` because archived source and
  executable identities retain original absolute checkout paths and there is no
  archive-relative source-snapshot replay;
- `independent_external_reproduction_satisfied: false`;
- `full_cost_gate_passed: false`;
- Koblitz index-calculus SOTA: false.

Full-cost blockers remain: no measured single-core elapsed time; no simultaneous
aggregate Cargo build-tree RSS; prior dependency acquisition in the Cargo cache
is outside the receipt; no relation matrix or modular linear algebra; no
scalar-hidden end-to-end IC run; no matched automorphism-optimized Pollard-rho
arm; and no independent external reproduction or novelty review.

The supported result is a completed finite public-synthetic n=23 exact
two-summand relation-yield candidate with successful planted and exact-miss
controls. It is not a SAT-solver speed result, relation-matrix or linear-algebra
result, scalar-hidden DLP recovery, rho comparison, scaling law, external
reproduction, novelty finding, deployed-key result, or Koblitz index-calculus
state-of-the-art claim.

The compact machine-readable companion is
`stage-21-relation-yield-result-summary-20260910.json`.
The self-hashed `stage-21-relation-yield-result-seal-20260910.json` binds this
report, that summary, the measured source commit, and the raw production run,
outer-receipt, and verification hashes. It does not substitute for the pending
independent factor-base and canonical-pair replay.
