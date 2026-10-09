# Modal VM CPU telemetry result

The one frozen, source-free Modal VM probe produced a complete guest receipt and exposed process perf cycles/instructions, CPU affinity, topology, NUMA placement, steal ticks, and interrupt counts. Its predefined guest observability score was **5/7**. The host isolation gate remains false, so the controlled CPU speedup field is null.

| Item | Observation |
| --- | --- |
| Freeze commit | `0599bf540` |
| Sandbox | `sb-oaq34hMRw2UhcNoyMCJalI`; remote exit 0; later read-only `poll()` returned 0 |
| Requested VM | 16 physical CPU cores, 32 GiB memory, 120-second timeout, blocked network |
| Elapsed local receipt interval | 55.343 s; 3.002 s idle snapshot interval; 5.035 s busy snapshot interval |
| Guest CPU | AMD EPYC 9654 96-Core Processor; 32 visible and allowed logical CPUs; selected logical CPU 0 |
| Guest affinity/topology | Pin to CPU 0 held through busy interval; guest reports one thread sibling for that logical CPU |
| Guest NUMA | One visible node, CPUs 0–31; process memory policy `default`, 3,211 observed pages on node 0 after busy interval |
| Process perf | 18,291,454,531 cycles; 63,019,318,404 instructions; raw instructions/cycles = 3.4453; no perf open errors |
| Busy interval selected CPU | 0 reported steal ticks; 1,303 interrupts; 66 involuntary process context switches |
| Cgroup throttling | `cpu.stat` reported zero additional throttled microseconds; `cpu.max`, partition and effective exclusive cpuset were unavailable |
| PSI and frequency | System CPU/memory PSI and guest CPU frequency were unavailable |
| Result gate | Guest observability 5/7; host-wide exclusivity unverified; controlled speedup unknown |

`receipt.json` preserves all snapshots, including unavailable fields, and `analysis.json` contains the predetermined score and counter deltas. The freeze and remote source SHA-256 values in the receipt match `FREEZE.json`; remote exit and receipt integrity checks passed. The runner attempted Sandbox creation once and called termination in its `finally` block. No cryptanalysis source, curve fixtures, targets, candidate tables, or repository archive were transferred.

The raw process perf ratio is an availability diagnostic. This probe did not group the events or record their enabled/running time, so the number is not a calibrated workload IPC result. Guest CPU numbering and topology are observations inside the VM; they do not prove which physical host core, SMT sibling, or other tenant shares the hardware. The recorded interrupt and context-switch activity further prevents treating this one interval as an isolated timing baseline.

For the S3 panel, a Modal run can still serve as a functional or exploratory screen if its inputs are approved for transfer. A controlled ≥2× CPU online claim requires the repository's host-wide isolation receipt, paired target solves, and noise gates. A dedicated host with auditable exclusive cores remains the appropriate venue for that gate unless Modal supplies equivalent host-level evidence.
