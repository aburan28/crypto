#set page(paper: "a4", margin: 1.8cm)
#set text(size: 9.5pt)
= Run report: sessions A and B (stride targets), 2026-10-14 17:00Z

*Study:* `research/ecbench_rootfind_speedup_20261014` (Conductor T-1). *Protocol:* `PROTOCOL.md`, committed in `07dd2d044` before these runs. *Binary:* `ecbench` SHA-256 `62bb5363f4dfca85…` (built at `07dd2d044`). *Host:* 4 vCPU Intel Xeon 2.10 GHz, shared cloud container; `--cpus none --allow-busy`, isolation *L0* on every record. *Runner:* `run.sh`.

#table(columns: 5, [session], [spec], [curve], [targets × rounds], [records],
  [A `ECBS1h66ed5f6ccf76`], [`claim-koblitz-rootfind-stride-n13.json`], [`icv1-f2m13-t181-515ee569`, r = 2003], [32 × 3 (+1 warm-up)], [768 / 768 verified],
  [B `ECBS1ha187c9428699`], [`check-koblitz-rootfind-stride-n19.json`], [`icv1-f2m19-t797-b6cf2467`, r = 130873], [16 × 2 (+1 warm-up)], [144 / 144 verified])

Audits (`verify --replay-all`): A ok, 768 records; B ok, 144 records (`audits/audit-*.json`); both re-audited ok by the Amendment 1 binary (`audits/reaudit-*.json`). No run failed or timed out.

== Measured S (GAE / √r), mean over measured runs, 95 % interval as `ecbench table` recorded

#table(columns: 4, [arm], [A (n = 13)], [B (n = 19)], [priced],
  [rho-strong], [15.260 [14.156, 18.247]], [3.112 [3.023, 3.205]], [lower bound],
  [rho-frob], [1.870 [1.774, 2.060]], [—], [lower bound],
  [d5-walk / d6-walk], [22.996 [22.186, 23.900]], [18.353 [16.517, 20.286]], [A full; B lower bound],
  [d5-stride / d6-stride], [13.436 [12.814, 14.326]], [21.036 [18.201, 23.838]], [A full; B lower bound],
  [d4-stride], [11.375 [10.522, 12.298]], [—], [full],
  [rho-strong-aa], [15.260 [14.156, 18.247]], [—], [lower bound])

== Paired comparisons (`ecbench compare`, candidate ÷ baseline, cluster interval)

- A: d5-stride ÷ d5-walk = 0.5843 [0.548, 0.629], 96 pairs on 32 workloads.
- A: d4-stride ÷ rho-strong = 0.7454 [0.618, 0.845] (diagnostic only: rho-strong has unpriced work).
- A: d5-stride ÷ rho-strong = 0.8805 [0.722, 0.991] (diagnostic only).
- A: rho-strong-aa ÷ rho-strong = 1.0000 [1.000, 1.000].
- B: d6-stride ÷ d6-walk = 1.1462 [0.960, 1.363], 32 pairs on 16 workloads (diagnostic only).
- B: d6-stride ÷ rho-strong = 6.7601 [5.853, 7.648]; d6-walk ÷ rho-strong = 5.8979 [5.242, 6.568].

Recorded counters (means over measured runs): A d5-walk 118.9 trials, rank 9.17, rows 9.18; A d5-stride 200.6 trials, rank 12.14, rows 14.94; B d6-walk 3,998.7 trials, rows = rank 16.28; B d6-stride 6,347.5 trials, rank 23.41, rows 29.56. Wall medians are L0 and descriptive only.

#figure(image("s-by-arm.svg", width: 100%), caption: [Mean S with the recorded 95 % interval per arm, sessions A (`ECBS1h66ed5f6ccf76`) and B (`ECBS1ha187c9428699`); circles fully priced, triangles lower bounds.])

Raw records: `sessions/claim-n13/records.jsonl`, `sessions/check-n19/records.jsonl`; tables and comparisons under `tables/`.
