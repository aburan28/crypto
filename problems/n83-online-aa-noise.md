# N83 cloud online A/A noise blocks finalist selection

- **Status:** `investigating`
- **Kind:** Empirical measurement noise; its cause is unknown.
- **Affected work:** [N83 range screens](../research/n83_online_range_20261016/RESULT.md), [mirrored rescreen](../research/n83_online_rescreen_20261016/PROTOCOL.md), and the five-point Tier 1 selection gate.
- **Next action owner:** N83 IC experiment campaign, Conductor T-3274.
- **Last updated:** 2026-10-16.

## Problem

The exact-root N83 `a=1` one-point online finalist cannot be selected from
the first two cloud screens because their predeclared K2543/root22 A/B
intervals differ by 24.70% and 25.93%, above the 10% gate. Both screens
used the same exposed public point, producer source SHA-256
`3ee6bc11399b382e72789b152043e5a67abc5b783c4e998ad8c534d60baebc03`,
fresh processes and exact-base replay. The selection status remains
`BOUNDED_UNKNOWN`; the five sealed Tier 1 targets stay unmeasured until a
replacement screening gate passes.

## Evidence and uncertainty

- **Observed:** [Run 38895743960](../research/n83_online_range_20261016/analysis_38895743960.json) measured 18.570 and 14.487 ms for K2543/root22 A/B (24.70%). [Run 38896107047](../research/n83_online_range_20261016/analysis_38896107047.json) measured 18.327 and 14.120 ms (25.93%). Both analyses passed independent replay of completed base, relation, scalar and rho rows. K3072/root11 measured 7.384 and 7.282 ms, while K2543/root11 measured 13.896 and 13.882 ms, respectively. K3584/root11 fixed reached the 900-second cap in both. Raw ZIP hashes, analysis hashes and byte-identical local replay logs are in the [result record](../research/n83_online_range_20261016/RESULT.md).
- **Derived or proved:** No cause of the A/B difference is established. The two observed differences exceed the frozen 10% gate arithmetically.
- **Suspected explanation:** Host scheduling, frequency or order effects could contribute; the current records do not separate these mechanisms.
- **Missing capability:** A repeat with adjacent controls, mirrored order and the same exact source, point and resource receipts.
- **Inference limit:** The observed 7 ms K3072 interval is a measured screen row; the two noisy control pairs do not admit it as the selected fresh-point pipeline.

## What would unblock it

- [ ] The [mirrored rescreen](../research/n83_online_rescreen_20261016/PROTOCOL.md) independently replays all eight arms and same-point rho, verifies equal base and exact-index identities within each K, and records no resource-gate violation.
- [ ] Both adjacent A/B and C/D online pairs for each K differ by at most 10% of their pair means; one candidate's slowest of four intervals is faster than the other's fastest by more than the largest adjacent A/A absolute difference.
- [ ] The winning base dump and exact IC1 manifest replay and hash before `selection.json` is versioned, and the five sealed targets remain untouched until that record is fixed.

These conditions admit a finalist for the separate Tier 1 panel. They do
not themselves measure that panel's five same-point IC/rho pairs.

## Next discriminating step

- **Action:** [Run 38905299528](https://github.com/aburan28/crypto/actions/runs/38905299528) executes the frozen mirrored K3072/root11 versus K2543/root11 screen.
- **Inputs and controls:** Same exposed point, source and rank seed; four fresh processes per candidate in adjacent A/A pairs and mirrored order; same-point rho; independent base, relation and scalar replay.
- **Limits:** 900 seconds and 15.25 GiB process RSS per IC arm, 120 seconds for rho, one Ubuntu 24.04 cloud allocation; retain each timeout and resource stop.
- **Outcomes:** A passing strict separation permits exact-manifest selection. A noisy, incomplete or overlapping result keeps selection unknown and calls for a Tier 1 instance screen with the same gates before any sealed-point run. An infrastructure failure remains an infrastructure status.
- **Recheck trigger:** The mirrored run's archived raw artifact and independent local replay.

## Resolution and history

2026-10-16: Opened after two independently replayed screens failed the same
10% A/A gate. The mirrored rescreen was frozen and launched before its
result was known. Unresolved.
