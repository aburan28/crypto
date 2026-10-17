# Volcano ECDLP census: Phase A status (setup, census, calibration)

Plan: [PLAN.md](PLAN.md). Frozen choices: [FROZEN.json](FROZEN.json). Phase B inputs: `specs/` (79 specs) and `runs/manifest.json`.

| item | where |
|:--|:--|
| calibration (b = 18..26, 9 arms, 8 targets, 3 rounds, all audits OK) | `calibration/` (`summary.md`, `sessions/`) |
| search plan, seeds, search log (every rejected seed) | `census/search_plan.json`, `census/seeds.json`, `census/search_log.jsonl` |
| census per volcano: `vertices.jsonl`, `edges.jsonl`, `transport/<src>.jsonl`, `walks/`, `diagram.{json,svg,png}` | `census/V1..V4/` |
| walk cross-check against `volcano-rho --curve` | `census/crosscheck.txt` |
| transported-instance identity (ecbench target vs ca-carried [k_i]G_src) | `checks/identity.jsonl` |
| tool patches | `volcano-rho-min-n.patch` (scratch worktree of the #2030 branch); `ca` changes on cryptanalysis branch `claude/volcano-census-export-20261017` (`--census-out`, `--edges-out`, `--ks`, test `tests/test_volcano_census.sh`) |

Not built: the all-l class census (`volcano-cli`, #2001): the binary is not in this tree or the #2030 branch. Conductor is not installed on this machine.

V4 (twist of V1) is censused but not measurable: its order is 2^6 5^2 17 1583, so the largest prime subgroup has r = 1583. No V4 specs were written; see the report to the coordinator.

V3 (j = 0 crater) has 16 vertices (1,1,2,4,8): the three 2-isogenies out of j = 0 land on one j. `ca` cannot step onto j = 0, so family U on V3 reaches 15 of 16 vertices (the crater is missing).
