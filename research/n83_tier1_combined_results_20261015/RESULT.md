# N83 sparse plus Bloom Tier 1 same-point result

The combined K3300/root11 sparse canonical selector and Bloom prefilter
completed five fresh, alternating, same-point IC/rho pairs on one GitHub
Actions cloud instance. All five independently replayed scalar answers and
four-point relations passed. The median paired rho/IC online ratio is
**35.103515×**; the five-target bootstrap 95% interval is
**[20.921769×, 1009.710039×]**. The wide interval reflects the distinct
target costs in this five-point panel; each ratio compares IC and rho on
the same public point.

| Fresh pair | Order | IC online (ms) | Rho online (ms) | Rho/IC | IC preparation (s) | Status |
| ---: | --- | ---: | ---: | ---: | ---: | --- |
| 1 | IC, rho | 41.635577 | 871.089943 | 20.921769× | 426.926 | `PASS_TIER1_PAIR` |
| 2 | rho, IC | 1.869869 | 1888.025501 | 1009.710039× | 414.430 | `PASS_TIER1_PAIR` |
| 3 | IC, rho | 47.471969 | 1666.432990 | 35.103515× | 417.643 | `PASS_TIER1_PAIR` |
| 4 | rho, IC | 28.522464 | 2518.838812 | 88.310702× | 413.798 | `PASS_TIER1_PAIR` |
| 5 | IC, rho | 108.052461 | 2597.706089 | 24.041156× | 417.829 | `PASS_TIER1_PAIR` |

The mean IC online interval is 45.510468 ms and the mean rho online
interval is 1908.418667 ms. Preparation is recorded per fresh IC process
and excluded from both online intervals. Each row has its own canonical
one-target workload and run ID in `workloads.json` and `measurements.json`.
The public points were frozen before execution, disjoint from earlier
Bloom-on and factor-base panels, and verified against the answer sidecar.

The exposed AMD EPYC 9V45 allocation had four affined CPUs, 16,373,452 KiB
visible memory, and a 14 GiB IC RSS guard. Peak observed IC RSS was
10.553 GiB. All five pair noise records report zero cgroup throttle and
OOM events and zero measured steal time. Calibration max/min was
1.120741, below the frozen 1.15 bound. This is Tier 1 cloud-instance
evidence; the host receipt preserves the exposed topology and cgroup data.

The candidate is
`IC1N83Ckb1a1fb547800PDP4rootRCguidedLAgaussTDdirectISO0h2aa3ee55ba07`:
547,800 distinct usable base points, 3,300 folded columns, 903,852,075
exact root-index keys, the sparse selector and three-hash Bloom prefilter.
The executed source is commit
`92e6f373312435bf89439b8a08bf49aa383989b2` with IC source SHA-256
`814afa03ae1c386e62e518be3c7a55fb108b09c41b05bef64ca1371787466c39`.

The raw producer is [Actions run 38683770144](https://github.com/aburan28/crypto/actions/runs/38683770144),
artifact `11950373884`, ZIP SHA-256
`6880939046d8980e188b0017d11ec1f2cbf25704b5757901de90e1d69d69473a`.
Independent analysis SHA-256 is
`b55cb95bfc9ec4379a9e7f80217d393f7d4abfc5a8720054bfd1ebc4cce4d466`.
The source-pinned CI workflow downloads these exact bytes, replays all
five pairs, and regenerates the one-target records on the current PR head.
