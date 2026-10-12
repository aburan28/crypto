# Shared GF(2) zero-table selector: exact but no measurable full-call gain

The [protocol](PROTOCOL.md) was committed before candidate code or timing. The candidate compared each table offset with its known zero-entry base instead of dividing the offset by the runtime suffix length. The tested implementation is preserved as [measured_candidate.patch](measured_candidate.patch) and its [CI workflow](WORKFLOW_CANDIDATE.yml) is archived. The opt-in runtime path and one-off workflow were removed after this result.

[CI run 36662335651](https://github.com/aburan28/crypto/actions/runs/36662335651) passed the release shared-kernel and F5 tests under both explicit flag values. The same release binary, SHA-256 `b14c46434dc00df19ce90ecbeb53f07c80d90b277972d80001e9448601116d19`, ran each paired reference/candidate process. Both one- and two-thread jobs qualified all four seed blocks. Each selected block had 22/22 successful processes and all seven F5 cases matched raw and canonical row fingerprints, rank, output terms, rows/columns, criterion/build counts and reduction word operations. Every selected reservation left zero eligible user threads and recorded zero contended samples at the frozen two-second interval. Three one-thread and one two-thread preflights refused PSI and produced zero benchmark calls; all were preserved before the first clean attempt was selected. The one-thread host was AMD EPYC 7763 and the two-thread host AMD EPYC 9V45; only within-job paired ratios are compared.

| Frozen n24 degree-4 complete call | One thread | Two threads |
| --- | ---: | ---: |
| Paired reference/candidate full-call median | **1.001×** | **1.008×** |
| Exact bootstrap 95% interval | 0.957–1.040× | 0.957–1.032× |
| A/A full-call ratio range | 0.972–1.009× | 0.999–1.027× |
| Paired elimination median | 1.004× | 1.017× |
| Reference full-call marginal median | 104.307 ms | 58.674 ms |
| Candidate full-call marginal median | 104.345 ms | 58.486 ms |

The three one-thread holdout full-call medians were 1.005–1.006×; the two-thread holdouts were 0.998–1.015×. No smaller case had a complete-call median below its own A/A minimum on either thread count, but the frozen one-thread primary misses the preset 1.05× median and 1.02× lower-bound gate. The observed change is inside A/A noise; no complete-call speedup or further 2× result is established. The suspected division cost was therefore not a useful independent target under this implementation and workload. The experiment did not measure generated instructions or cache misses, so it does not determine whether the compiler eliminated the divisions.

The deterministic [one-thread archive](runs/36662335651/MANIFEST-t1.json) preserves 31 raw files, including all three PSI refusals; the [two-thread archive](runs/36662335651/MANIFEST-t2.json) preserves 25 files. Both manifests list every file SHA-256, the compressed tar SHA-256, selected attempts, exactness, source/binary hashes and ratio summaries. [RUN_CANDIDATE.json](RUN_CANDIDATE.json) preserves GitHub run metadata. This is a matrix-F5 solver-stage diagnostic, not an IC one-target online DLP or rho speedup.
