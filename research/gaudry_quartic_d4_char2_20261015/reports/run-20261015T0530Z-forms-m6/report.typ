#set page(paper: "a4", margin: 1.7cm)
#set text(size: 9.3pt)
#set table(stroke: 0.4pt)
#show table: set text(size: 8.3pt)
= Run report: ecbench session `ECBS1h03a1d62aba88` (char2-forms-m6), 2026-10-15T05:30Z

Spec `specs/SPEC-char2-forms-m6.json` (id `ECS1hdcae87bf3a40`). Command `ecbench run --cpus auto --wait --allow-busy`; binary sha256 `92610bdf…`, commit `666438375`, clean. Host: Xeon \@ 2.10 GHz, 4 vCPU, run CPU 3 reserved. Finished 08:00Z.

Workload: curve `icv1-f2m24-t1171-7fbf6d93` ($E(FF_((2^6)^4))$, $log_2 r = 20.19$, $h = 14$), targets `W1240e58de0c0` and `Wabf66d7f9859`. One round. 10 executions requested, 10 completed, 10 verified.

#table(columns: 5, align: (left, right, right, right, left),
  [arm], [target 0 $S$], [target 1 $S$], [median wall (s)], [level],
  [rho-neg], [1.017], [1.314], [0.008], [L1/L2],
  [ic-plain], [713,704], [722,233], [4,517], [L1],
  [ic-weighted], [38,542], [31,128], [489], [L1],
  [ic-frobenius], [4,501], [3,654], [19.7], [L1],
  [rho-neg-aa], [1.017], [1.314], [0.008], [L2],
)

Recorded: 5,397 / 4,353 trials and 14 / 12 relations in every IC arm. `ecbench compare` (cluster, 2 targets): frobenius/plain 0.0057 (0.0051–0.0063); frobenius/weighted 0.117 (0.1168–0.1174); weighted/plain 0.0485 (0.0431–0.0540); frobenius/rho-neg 3,499 (2,781–4,427). Audit: `audits/char2-forms-m6.audit.json` (`verify --replay-all`): OK, 10 replays reproduced, 0 problems. Records: `sessions/char2-forms-m6/`.

#figure(image("s_forms_m6.svg", width: 72%), caption: [$S$ per arm and target as recorded, session `ECBS1h03a1d62aba88`.])
