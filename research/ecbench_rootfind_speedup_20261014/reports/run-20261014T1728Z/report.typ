#set page(paper: "a4", margin: 1.8cm)
#set text(size: 9.5pt)
= Run report: sessions C and D (stride-renew targets), 2026-10-14 17:28Z

*Study:* `research/ecbench_rootfind_speedup_20261014` (Conductor T-1). *Protocol:* `AMENDMENT_1.md`, committed before these runs. *Binary:* `ecbench` SHA-256 `9a45707c6a59a426…`. *Host:* 4 vCPU Intel Xeon 2.10 GHz, shared cloud container; `--cpus none --allow-busy`, isolation *L0* on every record. *Runner:* `run_amendment1.sh`, which first re-audited sessions A and B with this binary (ok, 768 and 144 records).

#table(columns: 5, [session], [spec], [curve], [targets × rounds], [records],
  [C `ECBS1h97bbbf34584f`], [`claim-koblitz-rootfind-stride-renew-n13.json`], [`icv1-f2m13-t181-515ee569`, r = 2003], [32 × 3 (+1 warm-up)], [896 / 896 verified],
  [D `ECBS1hb8d9ed954ea3`], [`claim-koblitz-rootfind-stride-renew-n19.json`], [`icv1-f2m19-t797-b6cf2467`, r = 130873], [16 × 3 (+1 warm-up)], [256 / 256 verified])

Audits (`verify --replay-all`): C ok, 896 records; D ok, 256 records. No run failed or timed out.

== Measured S (GAE / √r), mean over measured runs, 95 % interval as `ecbench table` recorded

#table(columns: 4, [arm], [C (n = 13)], [D (n = 19)], [priced],
  [rho-strong], [14.289 [14.148, 14.437]], [3.131 [3.038, 3.224]], [lower bound],
  [rho-frob], [1.927 [1.810, 2.122]], [—], [lower bound],
  [d5-walk / d6-walk], [23.508 [22.545, 24.504]], [17.552 [16.129, 18.902]], [C full; D lower bound],
  [d5-stride / d6-stride], [13.048 [12.525, 13.594]], [23.426 [20.989, 25.808]], [C full; D lower bound],
  [d5-renew / d6-renew], [12.483 [12.141, 12.823]], [14.657 [13.217, 16.134]], [C full; D lower bound],
  [d4-renew], [11.005 [10.297, 11.762]], [—], [full],
  [rho-strong-aa], [14.289 [14.148, 14.437]], [—], [lower bound])

== Paired comparisons (`ecbench compare`, candidate ÷ baseline, cluster interval)

- C: d5-renew ÷ d5-walk = 0.5310 [0.506, 0.557], 96 pairs on 32 workloads.
- C: d4-renew ÷ rho-strong = 0.7702 [0.720, 0.824] (diagnostic only: rho-strong has unpriced work).
- C: d5-renew ÷ d5-stride = 0.9567 [0.910, 1.004].
- C: d5-renew ÷ rho-strong = 0.8736 [0.847, 0.900] (diagnostic only).
- C: rho-strong-aa ÷ rho-strong = 1.0000 [1.000, 1.000].
- D: d6-renew ÷ d6-walk = 0.8351 [0.741, 0.943], 48 pairs on 16 workloads (diagnostic only: both arms leave lookups and row operations unpriced on this uncalibrated curve).
- D: d6-renew ÷ d6-stride = 0.6257 [0.545, 0.713] (diagnostic only).
- D: d6-renew ÷ rho-strong = 4.6813 [4.191, 5.202]; d6-walk ÷ rho-strong = 5.6058 [5.122, 6.082].

Wall medians are L0 and descriptive only.

#figure(image("s-by-arm.svg", width: 100%), caption: [Mean S with the recorded 95 % interval per arm, sessions C (`ECBS1h97bbbf34584f`) and D (`ECBS1hb8d9ed954ea3`); circles fully priced, triangles lower bounds.])

Raw records: `sessions/claim-renew-n13/records.jsonl`, `sessions/claim-renew-n19/records.jsonl`; tables and comparisons under `tables/`.
