# AutoLab configuration-search history

Eight development probes tested configuration changes against the original batch16 source/configuration. The frozen checker independently reverified 573 complete profile/native pairs and preserved three candidate timeouts. These development measurements do not decide promotion.

The enumeration candidate failed all three repetitions of one target; its aggregate cost and speedup remain unmeasured. Two probes also retained a broken-pipe error after their measurement result was saved. Those reporting failures do not erase the underlying completed measurements or the timeout evidence.

Primary unit: Valgrind 3.22 amd64 user-space instructions (Ir), from process startup through termination. All variants use the same four curve cells, eight public targets, fixed support and three repetitions. S = Ir / sqrt(subgroup order). The weak floor is K instructions for this full-rank K-column collector; it cannot establish an algorithmic advance.

| Experiment | Arm | S (Ir / sqrt(r)) | Cost / incumbent | Cost / rho | Cost / floor | Verified pairs | Class |
|---|---|---:|---:|---:|---:|---:|---|
| pair_table, batch 24 | incumbent | 568507 | 1 | 13.0047 | 2.04519e+07 | 24/24 | reference |
| pair_table, batch 24 | candidate | 615485 | 1.08263 | 14.0793 | 2.21419e+07 | 24/24 | engineering experiment |
| pair_table, batch 24 | rho | 43715.5 | 0.0768953 | 1 | unmeasured | 24/24 | reference |
| pair_table, batch 8 | incumbent | 568508 | 1 | 13.0047 | 2.04519e+07 | 24/24 | reference |
| pair_table, batch 8 | candidate | 529151 | 0.930772 | 12.1044 | 1.9036e+07 | 24/24 | engineering experiment |
| pair_table, batch 8 | rho | 43715.5 | 0.0768952 | 1 | unmeasured | 24/24 | reference |
| pair_table, batch 4 | incumbent | 568509 | 1 | 13.0047 | 2.04519e+07 | 24/24 | reference |
| pair_table, batch 4 | candidate | 525596 | 0.924516 | 12.0231 | 1.89081e+07 | 24/24 | engineering experiment |
| pair_table, batch 4 | rho | 43715.5 | 0.076895 | 1 | unmeasured | 24/24 | reference |
| pair_table, batch 2 | incumbent | 568511 | 1 | 13.0048 | 2.0452e+07 | 24/24 | reference |
| pair_table, batch 2 | candidate | 520445 | 0.915453 | 11.9053 | 1.87228e+07 | 24/24 | engineering experiment |
| pair_table, batch 2 | rho | 43715.4 | 0.0768946 | 1 | unmeasured | 24/24 | reference |
| pair_table, batch 1 | incumbent | 568508 | 1 | 13.0047 | 2.04519e+07 | 24/24 | reference |
| pair_table, batch 1 | candidate | 516610 | 0.908711 | 11.8175 | 1.85849e+07 | 24/24 | engineering experiment |
| pair_table, batch 1 | rho | 43715.7 | 0.0768954 | 1 | unmeasured | 24/24 | reference |
| pair_table, batch 1, window 8 | incumbent | 568511 | 1 | 13.0047 | 2.0452e+07 | 24/24 | reference |
| pair_table, batch 1, window 8 | candidate | 495388 | 0.871378 | 11.332 | 1.78214e+07 | 24/24 | engineering experiment |
| pair_table, batch 1, window 8 | rho | 43715.7 | 0.076895 | 1 | unmeasured | 24/24 | reference |
| pair_table, batch 1, window 4 | incumbent | 568508 | 1 | 13.0047 | 2.04519e+07 | 24/24 | reference |
| pair_table, batch 1, window 4 | candidate | 500771 | 0.880851 | 11.4552 | 1.80151e+07 | 24/24 | engineering experiment |
| pair_table, batch 1, window 4 | rho | 43715.6 | 0.0768954 | 1 | unmeasured | 24/24 | reference |
| enumerate, batch 1 | incumbent | 568508 | 1 | 13.0047 | 2.04519e+07 | 24/24 | reference |
| enumerate, batch 1 | candidate | unmeasured | unmeasured | unmeasured | unmeasured | 21/24 | engineering experiment |
| enumerate, batch 1 | rho | 43715.7 | 0.0768954 | 1 | unmeasured | 24/24 | reference |

[Frozen measurements and independent audit](measurements.json) contain source/binary checks, every receipt hash, full candidate configurations and retained errors. The audit used 21.76 host CPU-seconds, separately from benchmark container usage.

[Final campaign and current status](../README.md).
