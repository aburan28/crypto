# Balanced process order leaves K3300 and K3350 unresolved

The four-process comparison completed all 16,000 target answers and independently
replayed their scalars and four-point relations. The per-target contrast,
mean K3350 minus mean K3300, was **−0.429354 ms**, with a 95% contiguous
50-target block bootstrap interval of **[−2.057542, +1.227720] ms**. The
frozen selection rule records `BOUNDED_UNKNOWN` for the K3300/K3350 ranking.

| Fresh process order | K | Mean online ms per point | Setup s | Full process s | Peak RSS GiB | Replay |
| :--- | ---: | ---: | ---: | ---: | ---: | :--- |
| First | 3300 | 41.831879 | 577.367 | 803.092 | 10.053 | `PASS_INDEPENDENT_REPLAY` |
| Second | 3350 | 39.772037 | 552.940 | 768.030 | 10.198 | `PASS_INDEPENDENT_REPLAY` |
| Third | 3350 | 39.638486 | 551.181 | 766.047 | 10.199 | `PASS_INDEPENDENT_REPLAY` |
| Fourth | 3300 | 38.437351 | 531.613 | 740.534 | 10.050 | `PASS_INDEPENDENT_REPLAY` |

K3300's first-minus-second process contrast was **+3.394528 ms**, with a
95% interval of **[+3.280912, +3.507754] ms**. K3350's corresponding
contrast was **+0.133551 ms**, interval **[+0.069623, +0.198504] ms**.
The K3300 process shift is larger than the balanced mean difference. Both
bases remained within the 14 GiB RSS guard and recorded zero cgroup CPU
throttling and OOM kills. The cloud allocation exposed four AMD EPYC 9V74
CPUs. Online target intervals exclude the target-independent setup shown
separately above.

The common workload contains 4,000 fresh public N83 subgroup points, each
reused in all four arms. K3300 has 547,800 distinct usable base points and
3,300 folded columns; K3350 has 556,100 and 3,350. Source commit
`d5e6fc717bac432efba63564ba12081d464eebb6` built executable SHA-256
`318d2310ea4e42d21606159dbaf2a8afed38de8d4e2c7948fc3d969c612b7f71`.
The common workload ID is `2fc6da18140d`. Its K3300 and K3350 candidate
IDs are
`IC1N83Ckb1a1fb547800PDP4rootRCguidedLAgaussTDdirectISO0h8135b807be2f`
and
`IC1N83Ckb1a1fb556100PDP4rootRCguidedLAgaussTDdirectISO0h08423225ff50`.
[Actions run 38713270151](https://github.com/aburan28/crypto/actions/runs/38713270151)
passed source preflight, the four bounded processes, independent replay,
and raw upload. Artifact `11964534821` has ZIP SHA-256
`62f5c35620f993203a6842f1ac643c56b05067db881214d5298177ccfff65bdc`;
its `analysis.json` SHA-256 is
`65e66d4c5f3052eadc3dd3c23ca2df3d2203c67d10b08d2c7cba205936f59cde`.
The raw ZIP is archived on SSD2 under
`/Volumes/SSD990-2/n83-archive-20261015/raw/`. A separate local verifier
replayed all 16,000 answers and reproduced the analysis digest byte for byte.
`records_index.json` pins the compressed workload and measurements; the
current-head CI workflow regenerates both and requires a clean diff.

The two earlier 20,000-point orderings selected opposite bases. This
balanced run narrows the question: host or process-order effects remain large
enough to mask a K3300/K3350 timing difference in this allocation. The next
factor-base comparison should include a fresh same-point Tier 1 IC/rho
panel for a materially different K or root-index layout after its stage
candidate is selected.
