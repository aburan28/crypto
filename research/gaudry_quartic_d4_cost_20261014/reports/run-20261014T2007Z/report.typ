#set page(paper: "a4", margin: 1.7cm)
#set text(size: 9.3pt)
#set table(stroke: 0.4pt)
#show table: set text(size: 8.3pt)
= Run report: ecbench session `ECBS1h4687f1d786c2` (k4-p29-d4cost), 2026-10-14T20:07Z

Spec `specs/SPEC-k4-p29-d4cost.json` (id `ECS1h9deb3f4c9890`). Command `ecbench run --cpus auto --wait --allow-busy`; binary sha256 `67b8cc95…3e22`, commit `2b49f38d1`. Host: Xeon \@ 2.10 GHz, 4 vCPU, one reserved core. The ladder ran beside the session on the other cores.

Workloads: `fp4_torsion(29, 1, 16)` (ICV1 id in `plan.json`; not registered), $r = 44263$, targets `Wc87e324798ff` and `W1016932647b2` (seed 202610134029). One round, alternate order. 8 executions requested, 8 completed, 8 verified.

#table(columns: 6, align: (left, left, right, right, right, left),
  [arm], [method], [target 0 $S$], [target 1 $S$], [wall (s)], [level],
  [rho-neg], [rho.negation], [1.065], [1.412], [0.01], [L2],
  [ic-d4-committed], [`gaudry-d4:macaulay=14,charge_unsolved=1`], [657,144,598], [722,837,732], [7,321 / 7,928], [L1],
  [ic-d4-chains], [`gaudry-d4:macaulay=26,weighted=1,` \ `chains=1,sparsest=1,replay=1,` \ `charge_unsolved=1`], [1,653,062], [1,779,983], [221 / 235], [L1],
  [rho-neg-aa], [rho.negation], [1.065], [1.412], [0.01], [L2],
)

Recorded counters for the IC arms on targets 0 / 1: 574 / 635 trials in both arms; 9 / 8 relations; unsolved systems 9 / 8 (committed) and 11 / 13 (new); shape replays 568 / 625 with 3 / 4 fallbacks; solver ops $1.341 dot 10^13$ / $1.475 dot 10^13$ (committed) and $3.373 dot 10^10$ / $3.632 dot 10^10$ (new) $FF_p$ multiplications, priced at 97 per addition.

`ecbench compare` (`comparisons/`): ic-d4-chains / ic-d4-committed = 0.0025 in ops (cluster interval over 2 targets); ic-d4-chains / rho-neg = $1.39 dot 10^6$ (1.26–1.55). No wall comparison was admitted: no pair reached L2. Audit: `audits/k4-p29-d4cost.audit.json` (`verify --replay-all`): OK, 8 replays reproduced, 0 problems.

#figure(image("s_by_arm.svg", width: 72%), caption: [$S$ per arm and target as recorded.])
