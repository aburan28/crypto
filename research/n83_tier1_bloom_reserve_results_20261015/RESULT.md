# Five fresh N83 Bloom pairs pass the Tier 1 online speedup gate

The K3300/root11 Bloom candidate solved all five newly frozen public N83
subgroup points with independently replayed four-point relations and scalars.
Each same-point rho arm also independently recovered its scalar. All five
alternating pairs passed the Tier 1 correctness and noise gate on one GitHub
Actions Ubuntu 24.04 allocation. The median paired **rho online / IC online**
ratio was **104.23×**; the frozen five-pair target-bootstrap 95% interval was
**[53.12×, 220.74×]**. Every individual paired ratio exceeded one.

| Pair | Order | IC online | Rho online | Rho / IC online | Status |
| ---: | --- | ---: | ---: | ---: | --- |
| 1 | IC, rho | 21.790818 ms | 2271.284768 ms | 104.23× | PASS_TIER1_PAIR |
| 2 | Rho, IC | 13.190975 ms | 2911.733512 ms | 220.74× | PASS_TIER1_PAIR |
| 3 | IC, rho | 13.263169 ms | 948.840954 ms | 71.54× | PASS_TIER1_PAIR |
| 4 | Rho, IC | 6.981690 ms | 1539.165781 ms | 220.46× | PASS_TIER1_PAIR |
| 5 | IC, rho | 31.639253 ms | 1680.775037 ms | 53.12× | PASS_TIER1_PAIR |

The mean online times were 17.373181 ms for IC and 1870.360010 ms for rho.
These are one-target online intervals measured after each IC process prepared
its base, exact root index, Bloom bitset, and factor-log system. The IC
preparation costs, recorded separately for the five cold processes, were
680.676, 636.842, 637.017, 632.928, and 635.197 seconds. The largest
observed IC RSS was 11,329,925,120 bytes, within the 14 GiB guard.

The candidate is
`IC1N83Ckb1a1fb547800PDP4rootRCguidedLAgaussTDdirectISO0hc3b6db4ff8c0`:
547,800 actual usable subgroup points, 3,300 folded columns, base hash
`a4ba33098d2776520b98a326585a2301f8ccd5bc005df7c1a70fbe8496c61d8b`,
and 903,852,075 exact root keys. The host exposed an AMD EPYC 7763, four
selected execution CPUs with affinity set successfully, `cpu.max=max`, and
`memory.max=max`. The five pairs recorded zero cgroup throttling events and
zero steal ticks on those CPUs. Before/after calibration's global
maximum/minimum ratio was 1.0387, below the frozen 1.15 bound. These are
Tier 1 cloud-instance observations under the recorded allocation and noise
envelope.

The producer was GitHub Actions run `38678729913`, artifact `11949772579`.
Its raw ZIP SHA-256 is
`82be7c3f51faf55fc5d7d35ac63052b485f981267721e70650d46fae7ad9c097`;
the independent analysis SHA-256 is
`80850576f40572195c169b6f94126b190194c5e2aac10ad1f88d6110c28ef740`.
The executed source commit is `b70a2611fb8bc130c708ac91f26981592b2c1f78`,
and the IC source SHA-256 is
`814afa03ae1c386e62e518be3c7a55fb108b09c41b05bef64ca1371787466c39`.
The five-point public fixture SHA-256 is
`ef616c5254a9d883e962024b9f4bd7f794cfe5576c73a93b099ac5be6fe45668`.
The raw ZIP is durably archived at
`/Volumes/SSD990-2/n83-archive-20261015/raw/n83-tier1-bloom-38678729913.zip`.
The five canonical one-target workload and measurement rows beside this
report retain the exact online phases, preparation, resource and replay
receipts for each target.
