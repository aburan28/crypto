# Cold compact-orbit source pipeline: n37 admission result

**Decision: the preregistered full-rank gate was not met.** All eight unseen
public Q were recovered in both measured rounds by both IC arms, and all 48
measured answers (including strong rho) passed the runner's separate checks and
exact replay. The framework stops when the target scalar is uniquely pinned;
its relation matrix had rank **19–38 of 43 unknowns** (42 factor-base logs
plus the target scalar). This is a valid one-target ECDLP solve, but it is
not a full-rank factor-base log solve. The protocol therefore does not admit
a complete IC/rho route verdict or an ECC2K-130 transfer claim. The eager
table's 66,822 group additions alone still give a separate, unconditional
single-target setup floor for this exact n37 support.

The frozen [protocol](PROTOCOL.md), [spec](SPEC.json) and exact
[`ecbench plan` output](PLAN.txt) were pushed to
[PR #1283](https://github.com/aburan28/crypto/pull/1283) before this session.
The plan fixed 72 executions: eight registered ICV1 n37 public Q, three
arms, one warmup and two measured rounds. No target, support, trial cap,
unit or arm was changed after looking at results. The initial parity
failure and source-order repair are retained in [PARITY_FAILURE.md](PARITY_FAILURE.md).

## Reproducibility and correctness

| Item | Frozen record |
|---|---|
| Curve | `icv1-f2m37-tm534059-32aad96b`, `r=230603167`, signed-Frobenius order 74 |
| Source support | 42 columns, 3,108 signed points; every point, label, representative and BLAKE3 hash matched `SOURCE42.jsonl` in the release test |
| Source file SHA-256 | `0a32de24a5680ff46baf9543e8bf8e32447323491e2a2adedce2548c14a25f75` |
| Base BLAKE3 | `8423b135df3515b0e284126d9eeb95f71e3901ae12908b3cad8eb03d1f31d4bc` |
| `Cargo.lock` SHA-256 | `b28d3c2d81146a40d00df85f45c6460d9e2bc5c26875307a06cfd15efff99365` |
| Binary SHA-256 / commit | `66de6c2ffd5f36131d5d117932b0c9ca826a399f6c5229083498760f7facfb1a` / `e6a67bac2b6aa7e199a705bbd3f22794a0de764e` (clean) |
| Host | macOS aarch64, 14 logical CPUs, rustc 1.93.1, environment class `ECBENV2hd82681268e96`; all runs L0 |
| Session | [`ECBS1hc08f177f1f05`](sessions/n37_cold_v1/session.json), 72/72 verified including 24 warmups |
| Audit | [receipt](AUDIT.json), SHA-256 `7538c5af482b9528d67ed4a45cd038fa60eb90d6c3f12edc9db453480c69c479`; `ok=true`, 48/48 measured replays identical, zero problems |
| Post-lint replay | [second receipt](AUDIT_POST_LINT.json), SHA-256 `4755f70d7f532cd3081d35235e0c1f09ba99b4162c7f18fd387729bb79043cad`, auditor binary SHA-256 `f41c4c2db7498d4df5319bc3b62aad6176b04502256300f95c1f82cc850cb368`; 48/48 identical after the semantics-preserving Rust 1.98 `is_multiple_of` lint fix |

The two IC arms recovered the same scalar on all 16 measured (Q, round)
pairs. Each pair had the same factor-base identity, relation and matrix
counters, and group additions, doublings and scalar multiplications in every
phase. The only cost difference was the counted oracle's native setup:
66,780 canonicalisations, 111,888 Frobenius maps and 111,888 lookups,
yielding **8,210.630988 GAE** per cold table. Both table builds performed
66,822 group additions, stored 64,467 entries and used 42 representatives.
Every measured IC row had a verified scalar, 19–38 independent relation
rows, zero dependent rows, zero folded-table mismatches and no exhaustion.
Across the counted arm's 16 rows, all 491 relation trials yielded a relation.
The stop at `matrix.pinned(d_col)` in `collect_and_solve_with` explains the
short rank; no full-rank continuation mode was used.

## Charged operations on the exact same Q

`S = GAE / sqrt(r)`. These are **lower bounds**, because the legacy
half-trace lift, inversion, key scan, field checks, hash/index operations
and storage have no pinned common price; strong rho likewise leaves its
canonicalisations and table operations unpriced. L0 wall times are
descriptive only. The table and saved comparisons are in
[TABLE.txt](TABLE.txt) and the
[session comparisons](sessions/n37_cold_v1/comparisons/).

The complete IC unpriced set recorded on every run is
`inversions_uncharged`, `abscissae_scanned_uncharged`,
`curve_equation_checks_uncharged`, `factor_base_vec_bytes_floor_uncharged`,
`factor_point_inserts_uncharged`, `hash_probes_uncharged`,
`legacy_as_solves_uncharged`, `lift_field_multiplies_uncharged`,
`lift_field_squares_uncharged`, `min_key_field_squares_uncharged`,
`point_lifts_uncharged` and `point_negations_uncharged`.
Strong rho separately records `canonicalisations_uncharged`,
`partition_hashes_uncharged`, `table_inserts_uncharged` and
`table_queries_uncharged`. These names and their exact counts remain in
each sealed [record](sessions/n37_cold_v1/records.jsonl); the vector-byte
count is a memory floor, not an operation to add to GAE.
The folded table's normal-basis constructor and hash-table allocation do
not have dedicated operation counters in the current oracle, so the
recorded setup cost is also a lower bound beyond the listed native units.
`ecbench` reconstructs the deterministic base once more to report its
identity after the measured pipeline; that metadata work is excluded from
algorithmic GAE but appears in process wall time. No wall-time result is
claimed from this L0 session.

| Arm | Measured answers | Mean S [95% interval] | S / generic floor | S / same-Q rho | Class |
|---|---:|---:|---:|---:|---|
| Strong signed-Frobenius rho | 16/16 | 0.285 [0.240, 0.338] | 1.957 | 1.000 | reference, bounded |
| IC, legacy table accounting | 16/16 | 4.892 [4.874, 4.911] | 33.580 | 17.155 | accounting control, bounded |
| IC, counted table accounting | 16/16 | 5.433 [5.415, 5.451] | 37.291 | **19.050** [16.097, 22.663] | accounting, full-rank gate failed, bounded |

The counted/legacy IC ratio is 1.110515 [1.110107, 1.110934], with
identical algorithmic traces. The counted IC/rho paired ratio ranges
**10.796–32.149** across all 16 measurements. Exact per-Q GAE, phase
costs, scalars, ratios and ranks are in [PAIRED.json](PAIRED.json):

| Public workload | Counted IC / rho, rounds 1 / 2 | IC rank, rounds 1 / 2 |
|---|---:|---:|
| `W5b6a0d3da5cf` | 16.695 / 14.102 | 31 / 27 |
| `W5c6fc71fcae9` | 18.831 / 21.519 | 23 / 27 |
| `W5d12b17dadbb` | 15.817 / 24.463 | 34 / 38 |
| `W651d55346dc8` | 18.965 / 17.461 | 35 / 37 |
| `W790c41bebda8` | 10.796 / 21.909 | 29 / 32 |
| `W8a727ef60282` | 32.149 / 21.978 | 29 / 34 |
| `W928ce8e9f736` | 17.030 / 27.313 | 31 / 32 |
| `Wb555dea2e333` | 18.674 / 31.332 | 19 / 33 |

The source constructor scanned 105 raw abscissae, returned 89 lifted
points, passed 87 cofactor-projected candidates through the orbit key,
and selected 42 distinct orbits. Its counted group work plus pinned
Frobenius maps cost 4,302.554148 GAE. The counted oracle setup cost
75,032.630988 GAE, and mean online relation, linear-algebra and recovery
work cost 3,169.529466 GAE per Q. The recorded factor-base vector-storage
floor is 149,184 bytes; process peak RSS was 9,088–9,552 KiB for the
counted IC arm and 5,376–5,488 KiB for strong rho. Both figures are
measurements at this n37 size, not a memory scaling law.

The **66,822 setup additions alone** exceed the maximum recorded strong-rho
GAE on these 16 pairs (7,640.493947) by 8.746×, before native setup,
base construction, queries or rank. The minimum counted IC lower bound
(81,530.754545 GAE) exceeds that rho maximum by 10.671×. Those are
structural and sample-specific charged-operation inequalities for the
eager, one-target table; they do not repair the failed full-rank admission
condition, price missing native work, or establish wall-time speed.

## Next decision

Keep early target pinning as a separate, verified single-target variant;
the framework's current stopping rule is useful and should not be silently
called a full-rank solve. Preregister a full-rank continuation that keeps
collecting independent relations after the target is pinned, with a new
disjoint public-Q set and an explicit 43-column rank target. In parallel,
preregister a true shared-table batch comparison against matched
automorphism-aware **batched** rho. The arithmetic from these single-Q
means would place an *independent-rho* setup-amortisation crossing near
69 Q (`79,335.185136 / (4,330.868947 - 3,169.529466)`), but that is
only a hypothesis generator: the batch method, rank, memory, unpriced
native work and batched-rho denominator have not been measured. A new
protocol must fix the batch target count before outcomes. n41 and n53
remain counting-capacity follow-ups; no n131 inference is admitted without
the n83 confidence gate and explicit transfer assumptions.
