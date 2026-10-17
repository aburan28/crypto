# N83 root11 fine-grid screen selects K3072 for the fresh-point panel

The independently replayed 128-point screen selected **K3072/root11** under
the frozen minimum-geometric-mean online rule. It solved all 128 public points
with **509,952 usable factor-base points**, **3,072 folded columns**, and a
**16.405 ms** geometric mean target interval. The exact candidate is
`IC1N83Ckb1a1fb509952PDP4rootRCguidedLAgaussTDdirectISO0hec7ae07d8b9a`.
K3072 was already a fixed arm, so the 24-point Tier 1 holdout resolves to
three distinct candidates: K2543, K2816, and K3072.

| K | Usable points | Folded columns | Geometric online mean, ms | Median online, ms | Complete shared-setup process, s | Geometric online ratio K/K3072 | Descriptive 95% interval |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 2752 | 456,832 | 2,752 | 20.796 | 27.618 | 529.493 | 1.268 | [0.931, 1.727] |
| 2816 | 467,456 | 2,816 | 17.690 | 20.418 | **367.939** | 1.078 | [0.813, 1.430] |
| 2880 | 478,080 | 2,880 | 20.326 | 22.193 | 383.417 | 1.239 | [0.901, 1.703] |
| 2944 | 488,704 | 2,944 | 26.956 | 29.794 | 386.406 | 1.643 | [1.240, 2.177] |
| **3072** | **509,952** | **3,072** | **16.405** | **23.821** | 491.832 | 1.000 | reference |

The ratio interval is the exponential of the paired mean log difference plus
or minus 1.979 sample standard errors across the 128 public points. These
intervals are descriptive after selecting K3072 on those same points. K3072
was faster than K2816 on 64 of 128 individual points; their paired geometric
ratio K2816/K3072 was 1.078. The holdout's exact paired selector and cloud
envelope determine the fresh-point result. The complete shared-setup process
times include base construction, root-index and log preparation, all 128
queries, and producer checks; the five processes ran sequentially on the
ordinary macOS ARM64 host with 48 GiB RAM recorded in the [host receipt](single_target/host.json).

The [full replay](results/batch_full_replay.json) reconstructed all five
factor bases and verified every four-point relation, fixture scalar, recovered
scalar, timing phase, and source/binary binding. The
[frozen selection](results/local_selection.json) and [per-arm records](results/)
retain the 640 target intervals and probes. The public batch target SHA-256 is
`72486b7d708023838db1104a4a5e75a04b915d6f6d4013057c934fb31000aa87`;
the answer-sidecar SHA-256 is
`3ceee12fb848cc261433f3cd9e7d297f961effede1974e8c95fa25d9bac68720`.
The executed IC source SHA-256 is
`8a2fe55dd90dc33892c197d27260182bf783fc8212fb485d8297daf2bbba8e12`.
Raw base dumps, target results, commands, environment receipts, and the
offline release binary remain under
`/Volumes/SSD990/crypto/artifacts/n83_root11_finegrid_20261017/`; the
full-replay receipt records individual raw-output digests.

The [24-point protocol](PROTOCOL.md) now runs the three exact candidates on
four distinct Tier 1 cloud allocations. It pairs each IC run with a
same-point rho solve, alternates execution order, records before/after
calibration, and applies exact paired sign-flip tests to the independently
replayed single-target online and complete cold intervals.
