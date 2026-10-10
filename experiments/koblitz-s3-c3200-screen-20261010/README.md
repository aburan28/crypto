# Frozen N53 C1200/C3200 one-target paired screen

This screen pairs eight frozen N53 public points, T173–T180, between the
C1200 exact S3 reference and the C3200 adaptive quotient-root index. Each
arm starts with empty process caches, prepares its own complete base-log
system, then times one previously unseen target from first target-dependent
work through recovered scalar and native group check. The headline stage
comparison is paired `C1200_online_ms / C3200_online_ms`; root calls and
individual phases are diagnostics. Reusable setup, fixture construction,
input loading, and launch remain outside the online interval. The controlled
speedup field is `null` because this host lacks an exclusive CPU partition
receipt.

`SCREEN_FREEZE.json` binds both exact candidate manifests, source and binary
SHA-256s, target-point and workload hashes, relation seed, 14 preparation
workers, one target worker, 32 GiB RSS cap, 600-second per-run timeout, and
one repetition per arm and target. Arm order alternates AB/BA by target.
Every one-target solve has a fresh `run_one.py` helper process so
`RUSAGE_CHILDREN ru_maxrss` describes that solver rather than the maximum of
all previous solvers. `run_paired.py` writes each result and its compressed
raw report before starting the next run. It stops on the first failed gate
and preserves all completed and unstarted rows. Resume is available only for
a verified prefix of the frozen schedule.

The owner-set advancement threshold is an exploratory online geometric mean
of at least 2× across all eight verified pairs. The analyzer must also retain
raw failures, rank and relation costs, setup as a host-noise control, target
phase totals, and paired root-call ratios. A local ratio by itself is a
candidate-selection signal; an isolated host receipt is required to promote
it to a controlled CPU speedup result.

## Complete local result

All 16 one-target runs verified the fixture scalar and native group equality,
reached their frozen full preparation rank, recorded the five online phases,
and remained below 32 GiB peak RSS. The largest observed RSS was
30,930,141,184 bytes (28.81 GiB). `analysis.json` audits every receipt and
compressed report; `local/archive_receipt.json` binds all 67 raw files under
`local/paired_screen/` by SHA-256. Re-running the analyzer on that archived
directory reproduced `analysis.json` byte for byte.

| Point | C1200 online ms | C3200 online ms | C1200/C3200 online | C1200/C3200 S3 calls |
| --- | ---: | ---: | ---: | ---: |
| T173 | 0.266 | 0.533 | 0.499× | 14 / 14 |
| T174 | 0.540 | 0.498 | 1.084× | 692 / 80 |
| T175 | 0.490 | 0.607 | 0.806× | 708 / 142 |
| T176 | 0.254 | 0.328 | 0.774× | 42 / 42 |
| T177 | 0.379 | 0.909 | 0.417× | 422 / 242 |
| T178 | 0.707 | 0.700 | 1.010× | 816 / 108 |
| T179 | 1.285 | 3.638 | 0.353× | 2118 / 408 |
| T180 | 0.491 | 0.529 | 0.928× | 478 / 122 |

The target-resampled geometric mean online ratio is **0.682×**, bootstrap
95% interval **[0.514×, 0.884×]**, and the call ratio is **3.220×**,
interval **[1.836×, 5.522×]**. The online stage gate of 2× was not reached.
The unchanged C1200 reusable setup ranged from 0.58× to 1.13× its median;
C3200 ranged from 0.78× to 1.40×. Both fail the 10% setup noise
check, and the host has no exclusive partition receipt. Those wall ratios
are exploratory.

C3200 has 339,200 usable factor-base points and 3200 folded columns, versus
C2800's 296,800 points and 2800 columns. Its target S3-call count was
**identical to C2800 on all eight frozen points**. The extra 42,400 usable
points did not move the first verified relation in this panel. That exact
algorithmic result, together with the resource increase, directs the next
iteration to the target-check and large-index lookup kernels rather than
another small base expansion.
