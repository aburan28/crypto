# ECC arithmetic cross-platform performance is unqualified

- **Status:** `open`
- **Kind:** empirical bottleneck (host contention); implementation gap (unbuilt prime, Metal, CUDA and HIP paths); access limit (source export to Modal rejected by automatic review)
- **Affected work:** [ECC arithmetic plan](../research/ecc_arithmetic_20261013/PLAN.md), [architecture map](../research/ecc_arithmetic_20261013/ARCHITECTURE.md), [Apple M4 Pro run report](../research/ecc_arithmetic_20261013/reports/run-20261013T1837Z/report.md)
- **Next action owner:** ECC arithmetic implementation and measurement lane
- **Last updated:** 2026-10-13

## Problem

We cannot establish a fastest-in-class or even a controlled wall-time
advantage for the requested ECC arithmetic portfolio because the measured
Mac host was heavily contended, the Intel/cloud source-upload run was
rejected by automatic approval review, and matched implementations and
comparators are absent for several requested curve and GPU paths.

The requested scope includes NIST P-192/224/256/384/521, ECC2K-130 over
GF(2^131), Curve25519 and Ed25519, general binary and prime fields, and
Koblitz Frobenius scalar arithmetic on ARM64, Intel, Metal and a major GPU
backend. The current experiment covers only generic binary-field multiply
for four field degrees on Apple M4 Pro. The exact ECC2K curve, subgroup,
point operation and walk boundary were not timed.

## Evidence and uncertainty

- **Observed:** The [final paired CSV](../research/ecc_arithmetic_20261013/paired_mac_final.csv) has 24 AB/BA pairs for each of four fields. GF(2^131) baseline times span 27,826,500–932,235,083 ns per 500,000 products; candidate times span 16,801,958–694,301,375 ns. The unchanged four-limb control also fluctuates. The host reported load averages above 100 during timing and above 300 later. All candidate [focused tests](../research/ecc_arithmetic_20261013/candidate_tests_final.stdout) passed (17/17).
- **Observed infrastructure:** The first full-library release build hit [system-volume ENOSPC](../research/ecc_arithmetic_20261013/full_lib_tests.stderr). The second was stopped at load above 300 on 14 cores; its [log](../research/ecc_arithmetic_20261013/full_lib_tests_attempt2.stderr) is retained. A sparse checkout attempt lacked compile-time fixtures. The [sandboxed full suite](../research/ecc_arithmetic_20261013/full_lib_tests_attempt4.stdout) compiled and ran: 4,044 passed, nine local-socket tests failed with permission denied, and 114 were ignored. The [collaboration group](../research/ecc_arithmetic_20261013/pollard_collab_tests.stdout), including all nine, passed 37/37 with loopback access; a full rerun was [stopped](../research/ecc_arithmetic_20261013/full_lib_tests_attempt5_interruption.txt) under load above 150. An [x86-64 cross-check](../research/ecc_arithmetic_20261013/x86_cross_check_attempt2.stderr) compiled but did not run on Intel. A sandboxed Modal call [could not connect](../research/ecc_arithmetic_20261013/modal_cli_attempt1.stderr). [Automatic approval review](../research/ecc_arithmetic_20261013/MODAL_REVIEW_DENIAL.md) rejected a network-enabled retry because it would send the study source and `f2m.rs` to Modal without specific source-export authorization. No remote benchmark completed.
- **Derived or proved:** The new two-limb hardware product uses three rather than four carryless products; the three-limb version uses six rather than nine. This is a source-level operation count, not a runtime or whole-curve gain.
- **Suspected explanation:** Scheduling contention explains the large timing spread, but no profiler or isolation receipt separates it from thermal, frequency or cache effects.
- **Missing capability:** An isolated ARM host, an authorized and isolated x86-64 worker, matched Metal/CUDA/HIP kernels, fixed-limb prime-field implementations for the remaining curves, constant-time validation for secret-scalar paths, and same-host production comparators. The existing portable binary-field fallback has operand-dependent control flow and table access.
- **What this evidence does not establish:** A wall-time speedup, an Intel or GPU result, faster scalar multiplication, a faster ECC2K walk, or a worldwide record. Missing measurements are not regressions or negative mathematics.

## What would unblock it

- [ ] Run the frozen binary-field A/A and A/B panels on a quiet, pinned ARM64 core and an x86-64 host with recorded instruction features; retain all rounds, checksums, build hashes and a measurement-quality receipt meeting the repository's controlled comparison gate.
- [ ] Resolve specific authorization for exporting the minimal source bundle to Modal, or use a preauthorized worker with equivalent source and receipt custody. A rejected export must not be routed indirectly.
- [ ] Implement and differentially validate the listed prime, X25519/Ed25519, Koblitz, Metal and CUDA/HIP paths, then run same-boundary latency and batch panels against named production comparators on each relevant host.
- [ ] For any secret-scalar implementation, remove or isolate operand-dependent control flow and memory access, and validate the complete scalar path against its constant-time contract.
- [ ] For any worldwide claim, obtain a matched comparison on the same hardware and operation boundary, with uncertainty and independent replay; pass the applicable machine-readable claim and review gates.

The first condition would establish only a CPU binary-field result. The
remaining conditions are necessary to address the full request. None alone
establishes a global winner.

## Next discriminating step

- **Action or proof obligation:** Rerun the same source-bound binary-field panel on an isolated ARM64 host with an A/A control, then on an authorized x86-64 host; begin the fixed-limb Curve25519 and NIST-prime correctness path independently.
- **Inputs and controls:** Keep the 256 input pairs, seed, field polynomials, 500,000 products and 24 AB/BA rounds; use the unchanged four-limb path as an additional control, and check all outputs against the bitwise reference.
- **Execution limits and stop condition:** Use the fixed panel counts and the existing machine headroom and measurement-quality rules. Stop timing if the host is contended or correctness fails; preserve the failed receipt. This is not a campaign spending cap.
- **Outcomes and next actions:** A qualified gain proceeds to same-curve point and scalar tests. A qualified regression sends the implementation back for instruction/profile analysis. High A/A noise or infrastructure failure repeats only after isolation or access changes. A field-level result never promotes a walk or worldwide claim.
- **Recheck trigger:** Dedicated measurement host becomes available, source-export approval changes, or one of the missing backend implementations is validated.

## Resolution and history

2026-10-13: Opened with the Mac raw rounds, correctness tests, build failures,
and Modal denial. Unresolved. Completing this record or merging its PR does
not meet the original performance target.
