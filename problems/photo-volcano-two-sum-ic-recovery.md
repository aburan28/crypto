# Strict two-sum relation yield does not yet support random-corpus IC recovery

- **Status:** open
- **Kind:** empirical bottleneck and implementation gap; no proved impossibility
- **Affected work:** [frozen IC phase](../research/photo_volcano_restart_20261017/IC_PROTOCOL.md), [audited diagnostics](../research/photo_volcano_restart_20261017/reports/quantitative_1452/README.md), and full three-method completion of the photographed-volcano experiment
- **Next action owner:** experiment maintainer
- **Last updated:** 2026-10-17

## Problem

The implemented strict signed two-sum factor base has not produced the independent random relations needed to recover any of the 100 random target scalars on any of the 46 vertices.

The source has p=72057980207847001, a=18456521102804593, b=19958176276690439 and prime subgroup order r=18014494963342201. The workload is the frozen ten generators × ten paired cryptographically sampled scalars transported over the complete selected rational degree-2 component (46 vertices, 90 directed maps, five levels). This is not a census of all prime-degree isogenies. The acceptance criterion is genuine relation collection, factor-base logarithms, target descent and independently checked recovered ECDLP answers, with all work retained and measured. Generic rho/BSGS answers and planted controls do not satisfy it.

## Evidence and uncertainty

- **Observed:** 14,720 actual Gröbner systems were run: 9,200 random setup probes, 4,600 target probes and 920 planted positive controls. Random setup and target relation yields are zero, setup rank is zero, and full random-corpus IC recoveries are zero. All 920 positive controls pass. Singular independently checks equality of every returned ideal. See the committed polynomial replay archive and [phase report](../research/photo_volcano_restart_20261017/reports/quantitative_1452/README.md).
- **Observed:** a separate deliberately planted full-rank-12 pipeline passes polynomial solving → relation reconstruction → linear algebra → public point equation. It is a correctness control and contributes no main-corpus recovery.
- **Derived:** the measured signed two-sum support has 289 group elements including the identity, hence 288 nonidentity elements. For a target uniform over the nonidentity prime subgroup, support probability is 288/(r−1), approximately 1.5987×10^-14. Independent uniform probes under that law have expected waiting time approximately 6.255×10^13 trials for one relation. This is a support calculation, not an executed search or an impossibility proof.
- **Derived:** transporting the same factor base and target through an isogeny injective on this subgroup preserves every additive decomposition. These 46 copies cannot be counted as independent relation-yield experiments.
- **Missing capability:** an executable, independently verified relation collection and target-descent configuration that yields sufficient random relations on these fixed full-subgroup instances. The present native polynomial solver implements the small strict two-sum configuration; no larger viable configuration is established here.
- **What this evidence does not establish:** no lower bound against other IC constructions, no mathematical impossibility, no hardness separation between isogenous subgroups, and no claim that a planted control is a random solve. This gap does not stop the running generic solvers or their field-counter replay.

## What would unblock it

- [ ] On the exact frozen random workload, collect and independently check enough unknown-log relations for full required matrix rank; retain candidate generation, failures, arithmetic, elapsed/CPU time and memory costs.
- [ ] Complete genuine target descent and recover each requested scalar without feeding rho/BSGS answers or private sampling coefficients to the IC algorithm. Independently verify every public point equation and preserve every attempt.

These conditions would establish actual IC recovery on this workload. Qualified repeated isolated timing, canonical ecbench admission and strongest-reference comparison would remain separate obligations.

## Next discriminating step

- **Action:** specify a new additive decomposition/factor-base configuration and its exact enumerable support or justified relation law before expanding collection. Retain the existing phase unchanged; any executable extension gets a separately frozen protocol and controls.
- **Inputs and controls:** the same public 4,600 destination P/Q instances; the current support enumeration and positive/negative polynomial controls; no scalar witness input to the solver.
- **Execution limits and stop condition:** use the existing host/resource and archival contracts. Finish each declared finite diagnostic batch and retain its results; do not silently invent a spending cap or discard failed candidates. No larger relation batch is claimed as already executed.
- **Outcomes:** useful rank and genuine descent permit a full frozen IC run; zero or inadequate yield remains a bounded negative measurement; infrastructure failure is recorded separately. None substitutes generic ECDLP solves for IC.
- **Recheck trigger:** a concrete new decomposition/collector or independently demonstrated random-relation configuration for this exact subgroup.

## Resolution and history

2026-10-17: record the completed polynomial/ideal audits and exact small-support explanation. Full rho is complete; baseline BSGS and native arithmetic replays continue independently. Unresolved. Writing or merging this record does not complete the requested experiment.
