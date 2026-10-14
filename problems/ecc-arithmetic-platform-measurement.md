# ECC arithmetic cross-platform performance is unqualified

- **Status:** `open`
- **Kind:** empirical bottleneck (Mac contention and unnamed, unisolated Modal CPU); implementation gap (remaining prime, point, Ed25519, Metal, CUDA and HIP paths); access limit (new X25519 Modal source bundle)
- **Affected work:** [ECC arithmetic plan](../research/ecc_arithmetic_20261013/PLAN.md), [architecture map](../research/ecc_arithmetic_20261013/ARCHITECTURE.md), [Apple M4 Pro run report](../research/ecc_arithmetic_20261013/reports/run-20261013T1837Z/report.md), [Modal x86-64 report](../research/ecc_arithmetic_20261013/reports/modal-20261013/REPORT.md), [X25519 validation](../research/ecc_arithmetic_20261013/x25519_validation/README.md), [P-521 validation](../research/ecc_arithmetic_20261013/p521_validation/README.md)
- **Next action owner:** ECC arithmetic implementation and measurement lane
- **Last updated:** 2026-10-13

## Problem

We cannot establish a fastest-in-class or controlled wall-time advantage for
the requested ECC arithmetic portfolio. The Mac host was heavily contended;
the completed Modal x86-64 panels lack a CPU model, isolation and A/A control,
and their absolute times shifted between workers. Matched implementations
and comparators are absent for several requested curve and GPU paths.

The requested scope includes NIST P-192/224/256/384/521, ECC2K-130 over
GF(2^131), Curve25519 and Ed25519, general binary and prime fields, and
Koblitz Frobenius scalar arithmetic on ARM64, Intel, Metal and a major GPU
backend. The current experiment covers only generic binary-field multiply
for four field degrees on Apple M4 Pro and two Modal x86-64 workers. The exact ECC2K curve, subgroup,
point operation and walk boundary were not timed.

## Evidence and uncertainty

- **Observed:** The [final paired CSV](../research/ecc_arithmetic_20261013/paired_mac_final.csv) has 24 AB/BA pairs for each of four fields. GF(2^131) baseline times span 27,826,500–932,235,083 ns per 500,000 products; candidate times span 16,801,958–694,301,375 ns. The unchanged four-limb control also fluctuates. The host reported load averages above 100 during timing and above 300 later. All candidate [focused tests](../research/ecc_arithmetic_20261013/candidate_tests_final.stdout) passed (17/17).
- **Observed Modal:** After explicit user authorization, the first [Modal receipt](../research/ecc_arithmetic_20261013/modal_result.json) failed at an incorrect remote source path. A corrected launcher preserved the same six Rust inputs and completed [two x86-64 workers](../research/ecc_arithmetic_20261013/reports/modal-20261013/REPORT.md): 34 tests passed on each, 384 arm rows completed, and all paired checksums matched. The workers reported PCLMULQDQ but `cpu_model: unknown`, no reserved-core receipt, and substantially different absolute times between runs, including the unchanged four-limb control. The [raw CSV](../research/ecc_arithmetic_20261013/modal_paired.csv) and [graph](../research/ecc_arithmetic_20261013/modal_paired.svg) retain every round.
- **Observed X25519 implementation:** A fixed radix field and ladder passed [1,000 local field rows and 100 full-operation rows](../research/ecc_arithmetic_20261013/x25519_validation/README.md) against independent integer checkers, plus two RFC 7748 vectors in the standalone driver. A smaller [source-bound Cargo harness](../research/ecc_arithmetic_20261013/source_tests/README.md) passed 13 focused tests across the 25519 and P-521 modules. This has no timing result or machine-code side-channel certificate. The large-library release build [stopped under host load](../research/ecc_arithmetic_20261013/x25519_validation/CARGO_RELEASE_ATTEMPT.md) before emitting a test result.
- **Observed P-521 implementation:** A nine-word Mersenne field passed [1,000 local rows](../research/ecc_arithmetic_20261013/p521_validation/README.md) against independent Python integer checks of canonicalization, addition, subtraction, multiplication and inversion. The source is exposed as a module but not integrated into P-521 point or scalar arithmetic. No timing or compiled-code side-channel claim follows.
- **Observed access limit:** [Automatic approval review](../research/ecc_arithmetic_20261013/x25519_validation/MODAL_REVIEW_DENIAL.md) rejected the new eight-file X25519 source export before a worker started; the earlier seven-file binary-field authorization did not extend to this bundle. Conductor's local conflict/status endpoints initially timed out, then accepted the exact P-521 source and module scopes. Neither infrastructure event is mathematical evidence.
- **Observed infrastructure:** The first full-library release build hit [system-volume ENOSPC](../research/ecc_arithmetic_20261013/full_lib_tests.stderr). The second was stopped at load above 300 on 14 cores; its [log](../research/ecc_arithmetic_20261013/full_lib_tests_attempt2.stderr) is retained. A sparse checkout attempt lacked compile-time fixtures. The [sandboxed full suite](../research/ecc_arithmetic_20261013/full_lib_tests_attempt4.stdout) compiled and ran: 4,044 passed, nine local-socket tests failed with permission denied, and 114 were ignored. The [collaboration group](../research/ecc_arithmetic_20261013/pollard_collab_tests.stdout), including all nine, passed 37/37 with loopback access; a full rerun was [stopped](../research/ecc_arithmetic_20261013/full_lib_tests_attempt5_interruption.txt) under load above 150. An [x86-64 cross-check](../research/ecc_arithmetic_20261013/x86_cross_check_attempt2.stderr) compiled. A sandboxed Modal call [could not connect](../research/ecc_arithmetic_20261013/modal_cli_attempt1.stderr); [automatic approval review](../research/ecc_arithmetic_20261013/MODAL_REVIEW_DENIAL.md) then rejected source export until the user explicitly authorized it. That access condition has been cleared.
- **Derived or proved:** The new two-limb hardware product uses three rather than four carryless products; the three-limb version uses six rather than nine. This is a source-level operation count, not a runtime or whole-curve gain.
- **Suspected explanation:** Scheduling contention explains the large timing spread, but no profiler or isolation receipt separates it from thermal, frequency or cache effects.
- **Missing capability:** An isolated ARM host, an identified and isolated x86-64 host, specific approval for the new X25519 source export, matched Metal/CUDA/HIP kernels, fixed-limb prime-field implementations for P-192/224/384, P-521 point integration, Ed25519 fixed-field integration, constant-time validation for secret-scalar paths, and same-host production comparators. The existing portable binary-field fallback has operand-dependent control flow and table access.
- **What this evidence does not establish:** A qualified wall-time speedup, faster scalar multiplication, a faster ECC2K walk, or a worldwide record. Missing measurements are not regressions or negative mathematics.

## What would unblock it

- [ ] Run the frozen binary-field A/A and A/B panels on a quiet, pinned ARM64 core and an x86-64 host with recorded instruction features; retain all rounds, checksums, build hashes and a measurement-quality receipt meeting the repository's controlled comparison gate.
- [ ] Finish the X25519 release library and protocol tests; verify the compiled secret-scalar path and run isolated same-boundary A/A and A/B panels.
- [ ] Implement and differentially validate the listed prime, Ed25519, Koblitz, Metal and CUDA/HIP paths, then run same-boundary latency and batch panels against named production comparators on each relevant host.
- [x] Open the P-521 file scope after Conductor accepts a conflict check and lease; exact source and module claims were accepted as T-2518 and T-2520. Point integration and validation remain open above.
- [ ] For any secret-scalar implementation, remove or isolate operand-dependent control flow and memory access, and validate the complete scalar path against its constant-time contract.
- [ ] For any worldwide claim, obtain a matched comparison on the same hardware and operation boundary, with uncertainty and independent replay; pass the applicable machine-readable claim and review gates.

The first condition would establish only a CPU binary-field result. The
remaining conditions are necessary to address the full request. None alone
establishes a global winner.

## Next discriminating step

- **Action or proof obligation:** Rerun the same source-bound binary-field panel on isolated, identified ARM64 and x86-64 hosts with A/A controls; begin the fixed-limb Curve25519 and NIST-prime correctness path independently.
- **Inputs and controls:** Keep the 256 input pairs, seed, field polynomials, 500,000 products and 24 AB/BA rounds; use the unchanged four-limb path as an additional control, and check all outputs against the bitwise reference.
- **Execution limits and stop condition:** Use the fixed panel counts and the existing machine headroom and measurement-quality rules. Stop timing if the host is contended or correctness fails; preserve the failed receipt. This is not a campaign spending cap.
- **Outcomes and next actions:** A qualified gain proceeds to same-curve point and scalar tests. A qualified regression sends the implementation back for instruction/profile analysis. High A/A noise or infrastructure failure repeats only after isolation or access changes. A field-level result never promotes a walk or worldwide claim.
- **Recheck trigger:** Dedicated identified measurement hosts become available, the X25519 source bundle is explicitly approved for remote export, or one of the missing backend implementations is validated.

## Resolution and history

2026-10-13: Opened with the Mac raw rounds, correctness tests, build failures,
and Modal denial. Unresolved. Completing this record or merging its PR does
not meet the original performance target.

2026-10-13: User authorized the frozen GitHub and Modal exports. The original
Modal launcher's path error was preserved; two corrected x86-64 field panels
completed. The access condition cleared, but the controlled-comparison and
broader implementation exit criteria remain open.

2026-10-13: Added a locally validated X25519 fixed-radix field and ladder.
The new source bundle's Modal export was rejected before launch, and repeated
Conductor timeouts prevented a P-521 source claim. Those distinct access
conditions and the unverified tests remain open.

2026-10-13: Conductor subsequently accepted exact P-521 source and module
claims. A local fixed-width P-521 field then passed 1,000 independent integer
checks including inversion. The P-521 point/scalar and measurement gates remain
open. The X25519 source-export condition remains open.
