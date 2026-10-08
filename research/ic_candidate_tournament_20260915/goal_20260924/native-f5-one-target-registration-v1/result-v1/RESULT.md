# Audited F5 end-to-end control on one disclosed n17 point

The sole preregistered F5 target execution **recovered and independently
verified** the discrete logarithm `14668` for the supplied public point
`[61889,74818]` on `EC1N17Ce1hdfbf24105ef5`. The original frozen checker
returned `PASS_NATIVE_SOURCE_BOUND_TARGET_EXECUTION_AUDIT`, with
`source_bound_execution_admitted=true` and `verified_recovery=true`. This is a
complete F5 solve on this disclosed educational point. It is **not** a fresh
paired tournament result or a CPU speedup claim.

| Evidence | Original outcome |
| --- | ---: |
| Registration seal | `fcc320ea595781d94ab9cfc9461759c8f639fc58505489c1c8e4715e28be7d2e` |
| Original target worker calls | 1, consumed before launch; no retry |
| Nested original F5 preparation auditor calls | 1, before online timing |
| Target attempts | 2: one proved negative, then one verified witness |
| Recovered scalar | 14668, independently replayed against the target point |
| Source-attested online interval | 6,501.300958 ms |
| Outer worker wall | 33,006.260708 ms, including preparation recheck and setup |
| Memory peak / operation counts | unknown / unknown |
| Same-point rho / online speedup | not measured / unknown |

The five exclusive target-dependent phases close exactly to the source-attested
online interval; the unsuccessful first PDP attempt is charged to this target.

| Online phase | Time (ms) |
| --- | ---: |
| Target query | 40.232959 |
| Target PDP, including the proved negative | 6,436.693790 |
| Target relation check | 24.355208 |
| Target descent | 0.009250 |
| Target recovery check | 0.009751 |
| **Total** | **6,501.300958** |

The controller's exact publication preflight, frozen source/build identity,
one-use consumed claim, worker/preparation-auditor PID ledgers, stdout/stderr,
timeouts and process-group drains all passed the original audit. The checker
reconstructed both seeded `aG+bQ` queries, their PDP verdicts, the factor-base
logs, scalar recovery and phase closure. It executed no native child itself.
The original execution and audit are copied byte for byte here, including the
first negative attempt; `execution-original/terminal.json` binds their file
inventory. The immutable source/build and complete archived capsule were
published in `../publication-v1/` **before** the execute call. The original
consumed claim is `consumed-original.json`.

A separate post-run correctness control used the checked repository Sage
launcher and a different elliptic-curve implementation to verify
`14668·[43693,23339] = [61889,74818]` and both subgroup-order checks. It
returned `PASS_INDEPENDENT_SAGE_SCALAR_REPLAY`; its script, output and checked
Sage runtime receipt are retained here. This replay is outside online timing
and does not change the original frozen audit or its performance class.

The original preparation was the separately sealed F5 512-query natural
panel: 129 verified witnesses, 127 accepted rows, rank 29/29 and all 29 logs
independently replayed. That preparation cost is target-independent and is not
hidden inside the 6.5013-second online interval. The worker's outer wall also
includes its preparation audit and setup; it is not an online comparison
denominator. The 53.606685-second independent target audit ran after the worker
and outside online timing.

The host was an unisolated Apple M4 Pro, so the wall observation is
**exploratory**. The target was deterministically generated without a scalar,
but the protocol does not certify it as previously unseen across all archives;
the audit keeps `fresh_paired_qualification=false`. It also keeps
`candidate_id`, `workload_id`, `run_id`, operation counts, memory peak,
`online_wall_ns` as a promoted headline, and `online_speedup` null. A frozen
same-point strong rho arm, a completed SAT pipeline, canonical identities and
host-level isolation/noise evidence are still required for comparative claims.
No consumed preparation or target capsule may be rerun to fill those gaps.
