# Source audit of the N83 a=1/K600 online lead

Requested evaluation: the retained a=1, 600-column online IC lead over the 52.928-bit subgroup, not the a=0 81-bit primary sweep. Audit performed against aburan28/crypto commit a8194d616ed8c04489ec46177d369f1b1670360a. GitHub main was independently fetched and still resolved to that commit (2026-10-12T00:35:23Z). This audit does not run a benchmark or promote a performance claim.

## Immediate experimental implications

1. The historical runner calls koblitz_rho_fixture 83 1 signed_frobenius 1 packed. Current reference rules explicitly reject packed as the primary vs_rho baseline: it stores every step, scans the signed-Frobenius orbit in polynomial coordinates and performs one inversion per step. The strong reference uses distinguished points, normal-basis canonicalization and batched inversion across 32 single-thread lockstep lanes. The historical 8.181 median does not measure the current strong control.

2. The existing fixture CLI cannot run its strong arm at n=83: main lines 1748–1774 assert backend=="packed" for n>63. The library supports wide-field strong rho. A new exact public-point adapter is necessary rather than changing only the CLI word packed to strong.

3. Two n=83 K1 constructors are different models. KoblitzCurve::new(1,83) constructs the legacy sparse-search modulus [0,2,4,7] and generator (6614244651398184860292015,8868260410029249391539529). known_n83_k1() hard-codes modulus [0,1,2,45] and another generator. Both have the same subgroup order. Never compare their raw polynomial coordinates without an explicit field/point isomorphism.

4. Current IC source includes batched wide inversions, optional sparse canonicalization and optional root Bloom lookup added after the October 7 archived runs. A new observation is a new source-bound run. It must not be reported as merely retiming the archived binary.

5. The wide IC producer exits zero even when targets_failed>0 (line 2393). The harness must reject missing relations, group_verified!=true, targets_solved!=1, non-full rank and unequal recovered scalars. Exit status alone is not a success gate.

6. Set KIC_TARGET_THREADS=1 for a defensible exclusive-wall phase split. Parallel target scanning sums relation-check nanoseconds from concurrent successful workers and subtracts them from wall, then clamps PDP at zero. Those are not exclusive wall phases.

## Exact workload binding

Curve: y²+xy=x³+x²+1 over GF(2^83), irreducible polynomial z^83+z^7+z^4+z²+1.
Subgroup order: 8569786107849059; cofactor: 1128547018.
Generator: (6614244651398184860292015,8868260410029249391539529).
Archived Q: (2655898345537166180532960,3881526307991673265292615).
Retained base: experiments/koblitz-single-target-n83-20261006/base_n83_K600.jsonl, 600 signed orbit columns, 99,600 points.
Declared base set digest: d4a4da26d6e7cb09367de922d9806631fe49989f70d378cae76f347f8d76a9f7.
The header carries no source commit/source SHA-256/producer revision. Freeze the actual file bytes and independently replay the point labels; do not treat its declared digest alone as byte validation.

The source code uses a u128 field word but u64 subgroup scalar. That is appropriate for this 53-bit subgroup, not the 81-bit a=0 subgroup.

## IC preparation and public target access

examples/koblitz_orbit_dlp_fast.rs, wide::run, lines 1931–2393:

- Reads the base header and target file before setup. Public target lines are JSON arrays; any other line is parsed as a legacy known-answer scalar. The new harness must require exactly one array of canonical decimal strings and forbid legacy scalar input.
- Constructs the curve and checks header subgroup order and exact polynomial low terms. Reconstructs point/label arrays. For retained input it trusts the declared base_hash rather than recomputing it.
- Builds normal basis and S3 solver, executes 2047 root/Frobenius self-check probes, then builds the compact root index.
- Acquires factor logs using guided relations [a]G-R_j independent of the supplied target. Parallel rank uses scheduling-independent per-column scalars. It waits for all workers, performs echelon insertion and solves the full-rank system.
- Preparation accesses target array length for capacities but not target coordinate values. Its actual curve/base/index/log calculations are target independent.
- No known-answer validation sidecar is opened in the native producer. The historical Python launcher opens fixture.json for curve metadata, but passes only BASE, TARGETS, rank seed and output path to IC; public-point rho input is an environment array without a scalar. Inherited environment remains a concern: the launcher does not clear KIC_TARGET_THREADS, KIC_ROOT_BLOOM, KIC_SPARSE_CANONICAL or dump knobs. Use a clean declared environment.
- The native wide producer does not export factor logs or rank rows. KIC_DUMP_RANK is implemented only in the narrow path. Current final-scalar correctness is strong, but independently replaying every rank row/log needs an added diagnostic export or a separate calibrated replay contract.
- Parallel rank_failures remains zero even when a per-column loop has failed extraction attempts. rank_attempts includes those attempts; rank_probes accumulates successful relation probes only because failed extraction returns None without its probe count. Do not label these counters complete operation work. Wall includes the retries.
- The default parallel rank path can retry indefinitely. Use an outer process wall and hard memory guard that retains UNKNOWN, failure and partial output.

The archived full reusable preparation is recorded as 305045.905083 ms; README's approximately 4.8-minute value describes rank acquisition only. New preparation cost must be measured on the current executable and actual core envelope.

## IC five-phase timer audit

lines 2252–2322:

1. query_started begins after Q is selected/constructed. query_hash_started encloses a point-derived scan-start hash. This is T_target_query.
2. target_pdp_started encloses extract128_parallel.
3. lift128 lines 1693–1720 accumulates relation-check time ONLY for its final addition/equality against Q. The first two additions, point-choice lookup and loops remain under PDP. The phase is therefore a narrow final-lift check, not the entire group-witness validation routine.
4. Single-thread PDP = extraction wall minus this accumulated nested check duration; that produces a nonoverlapping partition with small timer gaps. Parallel checks are concurrent sums, so subtracting them from wall is invalid as an exclusive-wall decomposition.
5. target_descent times summing orbit coefficients times factor logs modulo r.
6. target_recovery_check times a fresh fast [d]G==Q verification.
7. target_ms spans the full online interval. The phase sum differs from target_ms by under 0.001 ms on each historical run. The original replay permits max(0.02ms,online*1e-8) tolerance. Preserve the residual explicitly; do not claim exact equality.
8. JSON serialization and target-row output follow target_ms. They are excluded from online; targets_total/process_total include output. State the same output boundary for both arms.
9. Input curve/subgroup validation is assumed from the public fixture; the IC target stage does not independently check arbitrary Q membership before its timer. The new strong worker includes target membership/subgroup validation in its online query phase; it is explicitly reported and can be separately compared if the matched protocol declares prevalidated fixtures.

Root-index probes count returned partner candidates reaching canonicalization/lookup, not all S3 quadratic calls or group additions. Record those counters in their own unit.

## Rho setup boundary resolved for new source-bound runs

examples/koblitz_rho_fixture.rs lines 1444–1464 and 1535–1547:

The wide packed jump table explicitly calls raw_random_state128, which computes aG+bQ for every jump. Its timed setup is therefore target dependent. Under this source it belongs inside packed online. The wide backend/curve construction precedes that setup timer and is target independent.

Historical derived reports say setup dependence was unresolved because receipts did not bind the source. This source audit resolves the current implementation, not the archived binary retroactively.

Strong rho differs: jumps(seed) are known multiples of G only (koblitz_strong_rho.rs lines 1090–1104) and can be prepared before Q. solve() performs Q-dependent walk starts (cG+Q), all collision work and a candidate [d]G==Q test before return (lines 1109–1219). Exclude its curve/normal-basis/G-only jump preparation; include solve and a separately timed generic scalar replay.

The default strong policy is lanes=32, dp_bits=4, step_cap_factor=2000, escape_fruitless=false. These are vectorized single-thread lanes, not 32 operating-system workers. For the same online single-core envelope, strong and IC target_threads=1 use one CPU. IC's 12-worker reusable rank preparation remains separately disclosed and cannot be treated as a matched single-core cold cost.

## Correctness receipt limitations

Historical verify_n83_runs.py:
- independently implements GF(2^83) arithmetic;
- checks G/Q curve membership, order annihilation, [d]G==Q, equal scalar and equal public Q in both arms;
- checks status, scalar-input flags and timer sum tolerance;
- requires the historical target's exact 8,845,441 probe count;
- does NOT independently replay the four factor-base-point relation, every factor log or the 600 rank equations;
- rho_phase_sum_matches compares a sum to the same sum expression; it is an arithmetic tautology rather than validation of timer boundaries.

Do not run this legacy verifier unchanged on a fresh target or current batched source: its hard-coded probe count is a historical receipt check, not a general correctness gate. New verifier should read the exact supplied Q and replay group membership, rG=O, rQ=O, [d]G=Q, the actual four-point relation and phase sums, with explicit source/binary hashes.

The launcher treats ru_maxrss as bytes because it was written for macOS. On Linux ru_maxrss is KiB. Scale correctly or use an audited native wait4 runner. Producer peak_rss_bytes() executes ps and reads current RSS, not lifetime peak, despite the key name; prefer the outer process peak receipt.

## Minimal feasible evaluation

First freeze one exact public point and source/binary/input hashes, validate its field/subgroup offline, define an equal online CPU and common hard memory envelope, and retain preparation as a separate quantity.

Run a current-source IC arm on the exact retained K600 object with KIC_TARGET_THREADS=1 and explicitly fixed rank workers/optional knobs. Run the new strong public-point arm on the identical Q. A packed arm is optional historical context only. One verified pair can directly falsify the specific claim that current IC is faster than the current strong control on that Q; it cannot establish a target-distribution speedup.

For a runtime claim, prospectively freeze at least five paired rounds with interleaved/randomized order, A/A noise controls and qualifying host receipts. Use independent new targets to evaluate extraction-tail variation after the initial one-target question. Retain capped outcomes and successful/failed counts. Do not infer a fresh-target CI from re-running the archived Q.

The independent exact-point strong worker prepared here is n83_strong_point_eval.rs. It uses KoblitzCurve::new(1,83), verifies the legacy modulus/r/h/G, instantiates StrongRhoG<WideGf2,u64>, refuses scalar or multi-target inputs, reports preparation separately, includes Q validation/solve/generic replay online and exits nonzero for invalid/capped/unverified outcomes. Its small n=7/n=17 tests compare the narrow and wide strong trajectories and charges exactly. No heavy compile or benchmark was performed by this audit agent.

Native isolated_bench can reserve/monitor an arbitrary external worker command. It does not turn that worker into an ecbench sealed run or certify the current L2 online claim/physical-host tier. Full ecbench integration can remain a separate work item; the bounded run is still valid exploratory evidence with honest receipts.

## Source changes since archived run

GitHub path history retrieved:
- 0cfba41136ff3b8a7ff235cddd04f5b5a1da249f, 2026-10-11: batched wide IC inversions and Gf2_128 fold reduction.
- 4b6e09b78f010d72de866b7eed2236782c47f656, 2026-10-11: optional S3 root Bloom lookup.
- c74bb3684e7939131c5c132573f1324f178d593c, 2026-10-11: optional sparse canonicalization.
- e478b952a8908ab9d0b4d74ff916c87a41b695ec, 2026-10-11: native n53 compact-orbit comparison/replay.
All fresh records must cite their actual compiled source, not copy the archived October 7 binary hash.

## Primary source pins

Snapshot commit: a8194d616ed8c04489ec46177d369f1b1670360a.
Git blobs:
- examples/koblitz_orbit_dlp_fast.rs: 0cf5be61913c5a925f637c87dc82670cb020ba3b.
- examples/koblitz_rho_fixture.rs: 38e9223ff89d7f967d4bf9f0ddd290a66f25e5a7.
- experiments/koblitz-single-target-n83-20261006/run_n83_pairs.py: 546163f88b15e694b4fdd4ee256ef293c8818e37.
- experiments/koblitz-single-target-n83-20261006/verify_n83_runs.py: 0d34426ecb4171d89b6cd1dc49d44bb7eeabf8ee.
- src/cryptanalysis/koblitz_fast_arith.rs: 3ec1da4e6c6b101ebde0b8ae81862895fdcd5566.
- src/cryptanalysis/koblitz_index_calculus.rs: 82f63978252bd9f50c1e48d7f3dcccdf61e03a13.
- src/cryptanalysis/koblitz_strong_rho.rs: 88917d4dda83ef1bcc87c0b4d83cb3305ed6a6d0.
- src/cryptanalysis/semaev_decomp.rs: d6cae7366bfca318698e26d1df953749a1069324.
- src/cryptanalysis/wide_gf2m.rs: 649789b6c8aa06a6d7cb0de33b52d4cddf9fa50b.

Local exact source copies and source_manifest.json identify SHA-256 hashes. No remote repository was mutated and no paid compute was requested.

