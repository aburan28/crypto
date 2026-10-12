# G7 native arithmetic comparison: no accepted speedup

The full run on one AWS `g7.4xlarge` completed all correctness gates and all
18 timed pairs, following eight excluded warmups. GPU: NVIDIA RTX PRO 4500
Blackwell Server Edition, sm_120, driver 595.91.07; compiler CUDA 13.3.73.
Each sample completed 201,863,462,912 scalar updates at the original geometry.
Three alternating pairs per candidate were measured for each workload.

| Variant (engineering) | Benchmark candidate / matched control (B/s) | Ratio | DP34 candidate / matched control (B/s) | Ratio | Correctness / acceptance |
|---|---|---|---|---|---|
| Native square | 4.688480 / 4.756701 | 0.985658 | 4.632716 / 4.547623 | 1.018712 | pass / not accepted |
| Three-limb Karatsuba | 4.188512 / 4.629912 | 0.904663 | 4.148491 / 4.547263 | 0.912305 | pass / not accepted |
| Both | 4.040647 / 4.579588 | 0.882317 | 4.007592 / 4.550207 | 0.880749 | pass / not accepted |

Rates are medians in billions of completed scalar updates per second. Each
candidate has its own interleaved control; ratios use that matched control.
The reference controls vary over the session, so candidate medians must not
be divided by a control from a different candidate's pairs.

**No candidate meets the predeclared acceptance rule:** at least 1% higher
median throughput in both workloads, with every pair faster. Native square
improves DP34 by 1.87% but slows the benchmark by 1.43%. Karatsuba slows the
benchmark by 9.53% and DP34 by 8.77%; both flags slow them by 11.77% and 11.93%.
Defaults remain off. This falsifies a broad speedup claim for these candidates
on this G7 configuration; it does not establish their performance on G7e.

All four arithmetic, storage, shared-sigma and integration gate groups passed,
as did cross-binary checkpoints and full sorted DP corpus comparisons. The
field-product count remains 5.3125 per scalar update, ratio 1.0 to the control.
The generic-group square-root boundary, walk rules and operation accounting
are unchanged. These are engineering practicality measurements, not a new
ECDLP algorithm or a lower generic-group exponent. The historical RTX PRO
6000 result of 14.637530 B/s is on a different GPU and is not a matched
reference for this run.

[summary.json](benchmarks/native-candidates/g7/summary.json) freezes all sample
rates, counts, corpus hashes, binary hashes and paired ratios.
[result.json.gz](benchmarks/native-candidates/g7/result.json.gz) retains the
complete runner JSON with raw compiler, correctness and timing output.
[source.tgz](benchmarks/native-candidates/g7/source.tgz) is the exact uploaded
source bundle; its SHA-256 is in the summary. The run used the launcher before
the concurrent cleanup-hardening commit a361973; those later changes are
preserved in this PR. Both benchmark instance IDs received termination calls.

The first G7 allocation failed before compilation because the archive did not
contain the empty `generated` directory. The runner now creates that directory
before `make generate`; generation passes from a fresh extracted source bundle.
The corrected allocation completed the full GPU comparison. Instance cleanup
status is retained in [cleanup.json](benchmarks/native-candidates/g7/cleanup.json).
