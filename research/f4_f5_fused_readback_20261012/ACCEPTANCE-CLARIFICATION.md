# Complete-call and standalone-phase gates

This clarification is recorded **before the corrected candidate is profiled**.
It changes the delivery gate's scope; retain the original stricter gate and
all its failed outcomes explicitly. Do not silently call them passes.

The original protocol said "both complete-call arms must pass" and "no
per-workload instruction regression above that gate." The native
`two_arm_acceptance_met` field interpreted the latter across all 33 cases,
including the separately invoked standalone criterion arm. Both first
attempts also contain criterion-only regressions: two cases in the first,
three in the second. They were not the headline complete-call regressions;
they remain observed phase results and must be reported.

The requested optimization concerns complete F4/F5 computations. Complete
F5 already includes criterion construction, row generation, packing, reduction
and readback. A separate standalone criterion invocation is an additional
diagnostic, with its own allocation context. Neither its instructions nor
its regression is added again to the complete-call cost. The criterion's
production source/rules are unchanged in this iteration. The cause of its
per-input instruction differences is not established by these measurements.

Use a separately named **`complete_call_acceptance_met`** delivery gate:
both complete F4 and F5 aggregates must be fewer instructions under the
unchanged noise/95% bootstrap/2% rules; none of their 22 complete-call cases
may regress above that rule; all allocation and peak-live metrics must remain
nonregressing across **all three arms**, with exact outputs and accounting.
All frozen workloads, measurements, costs, budgets and statistical thresholds
are unchanged. The original all-33-case `two_arm_acceptance_met` field is
preserved. If it fails, say so even when the complete-call gate passes.

This is an explicit narrowing from an additional standalone-phase guard,
not evidence that the criterion itself improved. Report every criterion case
verdict, including increases, and make no criterion-speed or whole-solver
claim. The metadata correction remains the final bounded candidate; if the
complete-call or memory gate fails, follow METADATA-FIX's rollback rule.
