# ECC2K cycle escape correctness repair

Baseline: e07fd850596eb326a16542e034ec95832cc99e88.

This is a correctness change, not a performance claim. The baseline already
contains rule v2 (pairwise exclusions and both four-term tau relations).
The old September 21 review predates those changes. Remaining concerns are
history-dependent exits, unbounded standalone defaults, and overbroad claims.

Hypothesis: a history hint may safely avoid work only if it cannot choose an
exit. Validate the hinted raw cycle (at most eight steps), select a single
orbit-invariant anchor among vertices whose cyclic histories trigger the hint,
and advance the branch only at that anchor. Other hints leave the raw step
unchanged. This should make arrivals at different phases of the same detected
cycle leave through the same edge, allowing a bounded delay to coalescence.
It is not a promise of immediate equality for states with different histories.

Frozen validation: the F131 two-cycle containing [1184]P, all cyclic entry
phases and empty/full histories in a finite synthetic model, both tau four-term
families, and the existing cycle-rule/covariance and restart tests. Test unrelated
false hints explicitly. Host replay and packed selection must agree. Reject old
table checkpoints/corpora by versioning.

Success: no divergent exits in the deterministic regressions, host/packed
agreement, finite standalone guard, and clear residual-cycle limitations.
Stop: an inconsistent exit, failed coefficient replay, or device mismatch.
GPU throughput and full-DLP speedup remain unmeasured; historical table-walk
performance does not certify the new rule. No fleet rollout is part of this task.

Implementation refinement: the raw-cycle probe must decline to redirect any
cycle containing a distinguished point, otherwise it could skip a report.
This is covered at every vertex of the synthetic cycles. The same threshold
is supplied to reference replay and packed selection.
