#set page(paper: "a4", margin: 1.6cm)
#set text(size: 9.5pt)
#set table(stroke: 0.4pt)

= lit2026 W5 run report: run-20261011T2225Z

*Run.* `scripts/run_arms.sh` (arms) and `scripts/diag.sh` (diagnostics), started 2026-10-11 22:25 UTC; `costmodel` re-recorded 2026-10-12 after a JSON-format fix in the probe's output (same algorithm). Candidate commit `d22fb3b52` (branch `lit2026/w5-lane-batched-f4`); original arm: `f6_ic_probe` built at `4b6e09b78` (shared baseline worktree). Host: shared Apple M4 Pro (14 cores), 1-minute load 82.5–388.1 over the 100 arm samples (median 147.1). Isolation *L0*, unpinned. Raw records: `raw/`; Rust-derived summary: `summary.json` (`f4_lane_probe summarize`); graphs drawn from them by `../../plot.py`.

*Workloads.* The four `f6_ic_probe` cells of W4 and the eager-rewrite study, all Koblitz curves named by ICV1 slug, prime-order subgroup chosen by `KoblitzCurve::new`, standard-subspace factor base, $m = 3$, `SolverEngine::InheritedF4 { max_degree: 3 }`, targets $[k]G$:
K1/2#super[17] dim 6 (`icv1-f2m17-tm101-00378d4e`, 8 targets); K1/2#super[23] dim 7 (`icv1-f2m23-tm5197-1f85e9e1`, 4); K0/2#super[23] dim 8 (`icv1-f2m23-t5197-69e76b73`, 4); K0/2#super[31] dim 9 (`icv1-f2m31-tm90707-c95f16f5`, 1 target, budget 2000). Engines `f4` (inherited F4) and `f6p` (F6-IC with the pair-sum index), `--reps 3`, 5 rounds, arm order rotated each round.

*Arms.* `orig` = baseline binary; `prev` = candidate binary, switches off; `gate` = `KIC_F6_GATE_FIRST=1`; `lanes` = `KIC_F4_LANES=k` (k = the cell's target count: 8, 4, 4, 1) with `KIC_F4_INHERIT_ROOT=rref`; `rref` = `KIC_F4_INHERIT_ROOT=rref` alone (the sequential control the lane roots reproduce row for row). All 100 runs exited 0.

== Stage unit (word_xors, deterministic: identical in all 5 rounds of every arm)

#table(columns: 7, align: (left, left, right, right, right, right, right),
  [cell], [engine], [orig = prev], [gate], [lanes], [rref], [orig / gate],
  [K1/2^17], [f4], [15 015 618], [15 015 618], [48 972 167], [14 098 582], [1.0000],
  [K1/2^17], [f6p], [12 635 914], [8 061 140], [46 474 171], [11 600 586], [1.5675],
  [K1/2^23], [f4], [53 273 191], [53 273 191], [160 671 638], [44 709 930], [1.0000],
  [K1/2^23], [f6p], [51 918 433], [30 697 846], [159 273 715], [43 312 007], [1.6913],
  [K0/2^23], [f4], [158 512 775], [158 512 775], [268 819 024], [123 273 212], [1.0000],
  [K0/2^23], [f6p], [124 949 801], [71 049 603], [235 173 139], [89 627 327], [1.7586],
  [K0/2^31], [f4], [281 583 026], [281 583 026], [190 512 003], [190 512 003], [1.0000],
  [K0/2^31], [f6p], [278 237 634], [156 423 046], [187 048 754], [187 048 754], [1.7788],
)

The lanes arm's stage unit includes each lane's share of the bit-sliced root elimination (64-bit word equivalents). Its tree part (stage unit less root eliminations) equals the `rref` arm's in every cell and engine (`lanes_tree_equals_rref_tree: true`). K0/2#super[31] has one target, so its lanes arm builds no batch and equals `rref`.

== Times (L0)

Process CPU (user + system) per round, previous / arm, median [min–max] over 5 rounds: gate f6p 2.22 [2.11–2.28] (K1/2^17), 2.49 [2.35–2.51] (K1/2^23), 2.24 [2.16–2.29] (K0/2^23), 2.13 [2.04–2.52] (K0/2^31); gate f4 0.98–1.04 (medians). rref / lanes CPU medians 0.83–0.97 (lanes slower) where a batch was built. The original binary does not report CPU time. Wall per round, original / arm, is in `arms.svg`; the original/previous null control ranges 0.24–4.60 per round, so wall admits no claim at this load.

#figure(image("arms.svg", width: 100%), caption: [Arms, run-20261011T2225Z: stage-unit ratio original / arm (top), wall original / arm (middle, log scale, median and min–max of 5 rounds), process CPU previous / arm (bottom). L0, load 82–388.])

== Identity

Per target and arm in every round: hit, exhausted, reductions, splits, propagations, infeasible, geometric refutations and witnesses, partial refutations, group additions, pair-index builds and lookups — identical across all 5 arms and 5 rounds in all 8 cell/engine pairs (`verdicts_and_counters_identical: true`). The one relation (K1/2^17, target 3, base indices [21, 48, 57]) is identical across the candidate binary's four arms (`relations_identical_across_candidate_arms: true`); the original binary does not print indices.

== Diagnostics

#figure(image("lanes.svg", width: 100%), caption: [Root phase only (`diag-lanes-*.jsonl`): k targets' degree-3 roots, sequential against lane-batched, 5 interleaved repetitions. Left: lane word-equivalents / sequential root stage units; middle: wall lanes / sequential (median); right: share of pivot columns pivoted by one row in every lane, and pivot passes per column per lane.])

#figure(image("lockstep.svg", width: 85%), caption: [Splitting-tree nodes (assignments) reached by two or more targets, and the share of the stage unit spent at them (`diag-lockstep*.jsonl`); K0/2^31 with 4 targets.])

#figure(image("costmodel.svg", width: 100%), caption: [Degree-3 root Macaulay rows and columns, measured / generic Boolean count (`diag-costmodel.jsonl`).])

Profile (`diag-profile.jsonl`, one pass, wall timers inside the solver): basis specialisation takes 82–97 % of the solver's build time and root builds 1.0–16.8 % in the 8 cell/engine pairs.

== Paths

`raw/r<round>-<cell>-<arm>.jsonl` (first line: `uptime`, UTC), `raw/*.status`, `raw/diag-{costmodel,lockstep,lockstep-k0-31-4targets,profile,lanes-echelon,lanes-rref}.jsonl`, `summary.json`, `arms.svg`, `lanes.svg`, `lockstep.svg`, `costmodel.svg`.
