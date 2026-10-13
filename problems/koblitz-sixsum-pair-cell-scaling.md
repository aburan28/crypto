# Koblitz six-sum pair-cell scaling

- **Status:** `open`
- **Kind:** empirical bottleneck; implementation gap
- **Affected work:** [higher-arity/rho study](../research/koblitz_higher_arity_rho_20261013/README.md), ECC2K-130 relation extraction
- **Next action owner:** higher-arity extraction and ecbench wide-pipeline maintainers; unassigned
- **Last updated:** 2026-10-13

## Problem

The measured six-sum extraction path cannot approach matched strong-rho cost
because it retains the same pair-cell enumeration while its two-S3 replacement
charges three group additions per cell, and the native IC framework cannot yet
run the required n=83 confidence gate.

The bounded discriminator is Koblitz `a=0`, n=41,
`icv1-f2m41-tm2308219-7f48b14a`, subgroup order 549756390943, a
`compact-orbit-scan` base with 20 columns, one public target, one warmup and
three measured rounds. The affected acceptance criterion is competitive full
cold `S = GAE/sqrt(r)` against `rho.signed_frobenius_strong`, followed by a
registered n=83 gate before any ECC2K-130 transfer discussion.

## Evidence and uncertainty

- **Observed:** [session ECBS1h4f09e1612522](../research/koblitz_higher_arity_rho_20261013/sessions/n41-local-l0/session.json) verified all 20 executions and [replayed](../research/koblitz_higher_arity_rho_20261013/sessions/n41-local-l0/audit-replay-all.json) all 15 measured runs. B3 and B2 found the same 58 relations and visited the same 5,109,566 pair cells. B3 charged 15,328,698 additional point additions and measured `S=7.652`, `55.447x` rho's `S=0.138`; S3/lookups remain unpriced. An isolated same-cell probe did make B3 extraction 6.44x faster than B2.
- **Derived or proved:** on this implementation and workload, replacing S4 by the two-S3 chain changes constant-factor field work but not pair-cell count; charging three additions per visited cell explains the counted regression. No asymptotic lower bound is proved.
- **Suspected explanation:** useful progress requires sharing the constructed `R +/- s1` points across cells/targets or changing the collision/relation representation so materially fewer cells are visited. This is a design hypothesis, not an established construction.
- **Missing capability:** `ecbench fb` reports `no factor-base plug-in runs on a wide (n > 62) curve yet`; `ic.pipeline` has no wide dispatcher, so B3 cannot run at n=83.
- **What this evidence does not establish:** no scaling exponent, n=83 result, ECC2K-130 transfer, or mathematical obstruction; n=41 wall results are diagnostic and the full IC costs are lower bounds.

## What would unblock it

- [ ] Implement a verified extraction/relation representation that reduces or amortizes the charged additions per accepted relation, then beat the current B3 cold operation count on the same frozen n=41 public workload without adding unpriced work.
- [ ] Add full-width factor-base, oracle, coefficient and linear-algebra support for registered n=83 Koblitz workloads, with exact scalar replay and no truncation.
- [ ] Complete the preregistered n=83 comparison against same-target strong rho and show competitive full cold `S` with every material native unit priced.

The first condition is a local discriminator. The latter two are necessary for
the requested confidence gate; even all three would leave n=83-to-n=131
transfer and independent replication obligations.

## Next discriminating step

- **Action or proof obligation:** prototype a row/batch traversal that reuses each materialized `R +/- s1` across multiple `s2` cells, and measure charged additions per verified relation against B3.
- **Inputs and controls:** the frozen n=41 curve, public target, K=20 base, diagonal cell set, B2/B3 relation verifier, and strong rho reference from the linked protocol.
- **Execution limits and stop condition:** follow the existing 900-second per-run protocol; stop on correctness disagreement, resource refusal, completed audited session, or the protocol timeout.
- **Outcomes and next actions:** a lower full cold count advances to the wide implementation gate; only lower wall with unchanged/higher counts remains an engineering improvement; disagreement is a correctness failure; timeout is bounded inconclusive evidence.
- **Recheck trigger:** a reusable-point/batched-cell implementation or full-width IC pipeline lands.

## Resolution and history

2026-10-13: opened after the two-S3-chain implementation and audited matched-rho
session. The earlier broad blocker is narrowed to unchanged pair-cell work plus
the missing wide IC path. Unresolved.

