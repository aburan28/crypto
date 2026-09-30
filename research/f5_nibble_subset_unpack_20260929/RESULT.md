# Matrix-F5 nibble-subset unpack: rejected after qualified paired runs

The [protocol](PROTOCOL.md) was committed before candidate code or timing. The candidate stored all ordered subsets of each four-column group and copied selected monomial masks into exact output rows. The tested implementation is preserved as [measured_candidate.patch](measured_candidate.patch) and its [CI workflow](WORKFLOW_CANDIDATE.yml) is archived. The runtime candidate and one-off workflow were removed after this result. This is a matrix-F5 solver-stage experiment, not an IC one-target online DLP or rho comparison.

[CI run 36660474694](https://github.com/aburan28/crypto/actions/runs/36660474694) passed the release shared-eliminator and F5 tests, then produced a qualified same-binary paired result at one and two threads. Both jobs used release binary SHA-256 `13c23dcb45a3da44111c7b08ecfbb593fe0bef0663ab528a0b61282d4d5b00ba`. The one-thread host was AMD EPYC 7763; the two-thread host was Intel Xeon Platinum 8573C. Ratios compare only modes on the **same host and binary** within each job. Each selected seed block had 22/22 successful processes, seven cases per process, with matching raw and canonical row fingerprints, rank, output terms, column and row counts, criterion/build counts, and reduction word operations. Every selected reservation left zero eligible user threads and recorded zero contended samples at the predeclared two-second interval. The one-thread runner rejected five zero-call attempts for PSI preflight and then selected the first clean attempt per seed; the two-thread runner selected attempt one for every seed. All refusals and calls are archived.

| Frozen n24 degree-4 complete call | One thread | Two threads |
| --- | ---: | ---: |
| Paired reference/candidate median | **0.649×** | **0.649×** |
| Exact bootstrap 95% interval | 0.646–0.650× | 0.635–0.667× |
| A/A full-call range | 0.995–1.008× | 0.960–1.015× |
| Reference full-call marginal median | 110.40 ms | 88.15 ms |
| Candidate full-call marginal median | 169.98 ms | 135.08 ms |
| Reference output unpack marginal median | 55.26 ms | 34.14 ms |
| Candidate output unpack marginal median, table build included | 114.99 ms | 82.29 ms |

The candidate's actual table allocation was 1,049,112 bytes, below the frozen 2 MiB cap. On all four seeds, the primary complete-call medians ranged from 0.648–0.654× at one thread and 0.645–0.650× at two threads, all below their A/A minima. The n20 degree-4 case also regressed outside A/A on every seed: 0.690–0.698× at one thread and 0.638–0.684× at two. The unmodified route served the other five cases. The candidate misses the preset 1.05× median and 1.02× lower-bound gates and the smaller-case guard; it supplies no 2× result. Because this 1.05 MiB table was also substantially slower, the earlier byte-table regression cannot be attributed to table size alone. Extra per-chunk lookup and copy work is a plausible cause, but no cache or instruction attribution was measured here.

The deterministic [one-thread archive](runs/36660474694/MANIFEST-t1.json) preserves 37 raw files, including five PSI refusals; the [two-thread archive](runs/36660474694/MANIFEST-t2.json) preserves 22 files. Each manifest lists the SHA-256 of every raw file and the compressed tar, all selected attempts, exactness, ratios, source and binary hashes. [RUN_CANDIDATE.json](RUN_CANDIDATE.json) preserves GitHub run metadata. The unchanged-source [qualified baseline and binary archive](https://github.com/aburan28/crypto/pull/1013) are available separately; its absolute time was not used as a candidate denominator.
