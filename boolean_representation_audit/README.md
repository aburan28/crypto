# Boolean row representation audit: exactness passes, performance is conditional

2026-09-29 UTC. Standalone synthetic algebra, 6 and 8 polynomial variables. These are not field degrees.

The packed backend preserves every certified row and algebraic-work counter of the recovered sparse backend under each matching schedule. The receipt audit reconciles 504 fresh-process samples, 1,416 independently verified computations, and 48 fresh output recomputations. Six new test methods include 60 random systems and the 12 frozen systems; the original nine tests also pass.

Classification: representation engineering / stage diagnostic. No novel algebraic algorithm, asymptotic improvement, blocked-workload result, ECC experiment or end-to-end cryptographic improvement is established. Full ideal closure still has an exponential output ceiling. The m=83 confidence gate is unperformed.

## Decision

The exhaustive schedule meets the predeclared numerical four-holdout criterion; frontier scheduling fails it because the chain and cycle controls stay near parity. This is not a blanket packed-backend performance claim. Baseline-only A/A controls expose substantial host noise (one tiny-case pair reaches 27.88x), and process/NUMA controls were incomplete. Treat timing improvements as provisional observations on this virtualized host, not validated hardware rates. No post-hoc tuning or reruns selected for favorable timing were performed.

The useful structural observation is that sparse input does not imply sparse intermediates: the new sparse quadratic system starts with at most six terms per equation and reaches 50 during frontier elimination. Both new random quadratic systems have a complete equation-level interaction graph (min-fill width 7), unlike the width-1 chain and width-2 cycle. These are diagnostic examples, not a statistical relationship or optimal-width proof for arbitrary systems.

## Paired representation comparison

Each ratio is packed time / sparse time under the same schedule; lower is less time. Intervals are exploratory 95% percentile bootstrap intervals on seven paired worker means, without multiple-comparison correction. Compute includes setup, conversion, products, queueing, elimination and certificate creation. The separate total adds independent certificate verification; neither includes process startup/imports, JSON or fingerprinting.

| Case | Schedule | Compute ratio [95% interval] | Compute + verify ratio | A/A min–max | Python allocation ratio |
| --- | --- | ---: | ---: | ---: | ---: |
| monomial-6 | exhaustive | 0.678 [0.608, 0.799] | 0.820 | 0.763–1.141 | 0.905 |
| monomial-6 | frontier | 0.937 [0.888, 1.095] | 0.965 | 0.763–1.141 | 0.910 |
| block-linear-6 | exhaustive | 0.973 [0.862, 1.260] | 0.918 | 0.524–1.775 | 0.919 |
| block-linear-6 | frontier | 0.963 [0.946, 1.203] | 0.952 | 0.524–1.775 | 0.913 |
| planted-quadratic-6 | exhaustive | 0.672 [0.562, 0.733] | 0.721 | 1.004–2.242 | 0.905 |
| planted-quadratic-6 | frontier | 0.678 [0.603, 1.301] | 0.725 | 1.004–2.242 | 0.895 |
| dependent-contradictory-6 | exhaustive | 0.973 [0.683, 1.753] | 0.982 | 0.834–27.878 | 0.905 |
| dependent-contradictory-6 | frontier | 0.989 [0.751, 1.235] | 0.981 | 0.834–27.878 | 0.916 |
| monomial-8 | exhaustive | 0.788 [0.720, 1.178] | 0.776 | 0.759–0.963 | 1.002 |
| monomial-8 | frontier | 1.033 [0.778, 1.140] | 0.999 | 0.759–0.963 | 0.999 |
| block-linear-8 | exhaustive | 0.667 [0.652, 0.774] | 0.728 | 0.327–1.093 | 0.982 |
| block-linear-8 | frontier | 0.974 [0.902, 1.021] | 0.973 | 0.327–1.093 | 0.981 |
| planted-quadratic-8 | exhaustive | 0.288 [0.253, 0.304] | 0.343 | 0.617–1.406 | 0.919 |
| planted-quadratic-8 | frontier | 0.301 [0.285, 0.309] | 0.328 | 0.617–1.406 | 0.887 |
| dependent-contradictory-8 | exhaustive | 0.897 [0.750, 0.996] | 0.950 | 0.688–1.383 | 0.978 |
| dependent-contradictory-8 | frontier | 0.978 [0.851, 1.099] | 0.980 | 0.688–1.383 | 0.982 |
| holdout-chain-linear-8 | exhaustive | 0.575 [0.534, 0.588] | 0.613 | 0.964–1.116 | 0.977 |
| holdout-chain-linear-8 | frontier | 0.974 [0.933, 1.080] | 0.998 | 0.964–1.116 | 0.977 |
| holdout-cycle-quadratic-8 | exhaustive | 0.735 [0.700, 0.760] | 0.758 | 0.392–1.096 | 1.000 |
| holdout-cycle-quadratic-8 | frontier | 0.994 [0.884, 1.051] | 0.962 | 0.392–1.096 | 0.994 |
| holdout-sparse-quadratic-8 | exhaustive | 0.179 [0.176, 0.219] | 0.223 | 0.883–1.635 | 0.902 |
| holdout-sparse-quadratic-8 | frontier | 0.291 [0.282, 0.312] | 0.321 | 0.883–1.635 | 0.878 |
| holdout-dense-quadratic-8 | exhaustive | 0.203 [0.174, 0.208] | 0.232 | 0.810–0.967 | 0.867 |
| holdout-dense-quadratic-8 | frontier | 0.203 [0.199, 0.210] | 0.228 | 0.810–0.967 | 0.791 |

A/A controls use the sparse frontier schedule only. Broad A/A ranges block fine-grained timing conclusions. Fresh untraced worker RSS values span 14,848–15,360 KiB and include interpreter/import/verifier memory; they do not establish a process-memory gain. Python allocations are separate traced compute-only runs. Packed payload bytes are minimal bit payload, while sparse payload counts 32-bit indices; neither is actual process memory.

## Absolute stage times and deterministic work

Every row below is verified. Counters are algebraic row counts, not calibrated machine operations. The exact output and work counts match between backends of the same schedule; representation cost alone changes.

| Case | Variant | Median compute ms | Median compute + verify ms | Submitted | XOR rows | Rank | Stored terms | Peak row terms |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| monomial-6 | sparse/exhaustive | 0.355 | 0.428 | 64 | 48 | 16 | 16 | 1 |
| monomial-6 | packed/exhaustive | 0.247 | 0.326 | 64 | 48 | 16 | 16 | 1 |
| monomial-6 | sparse/frontier | 0.467 | 0.542 | 33 | 17 | 16 | 16 | 1 |
| monomial-6 | packed/frontier | 0.466 | 0.542 | 33 | 17 | 16 | 16 | 1 |
| block-linear-6 | sparse/exhaustive | 3.191 | 3.641 | 192 | 280 | 56 | 112 | 2 |
| block-linear-6 | packed/exhaustive | 2.943 | 3.289 | 192 | 280 | 56 | 112 | 2 |
| block-linear-6 | sparse/frontier | 5.678 | 6.047 | 187 | 223 | 56 | 112 | 2 |
| block-linear-6 | packed/frontier | 5.790 | 6.124 | 187 | 223 | 56 | 112 | 2 |
| planted-quadratic-6 | sparse/exhaustive | 4.299 | 4.907 | 192 | 859 | 52 | 208 | 12 |
| planted-quadratic-6 | packed/exhaustive | 2.858 | 3.466 | 192 | 859 | 52 | 208 | 12 |
| planted-quadratic-6 | sparse/frontier | 6.337 | 6.917 | 229 | 1253 | 52 | 232 | 20 |
| planted-quadratic-6 | packed/frontier | 4.363 | 5.018 | 229 | 1253 | 52 | 232 | 20 |
| dependent-contradictory-6 | sparse/exhaustive | 1.456 | 1.835 | 192 | 128 | 64 | 64 | 1 |
| dependent-contradictory-6 | packed/exhaustive | 1.283 | 1.554 | 192 | 128 | 64 | 64 | 1 |
| dependent-contradictory-6 | sparse/frontier | 2.466 | 2.738 | 195 | 132 | 64 | 64 | 1 |
| dependent-contradictory-6 | packed/frontier | 2.328 | 2.712 | 195 | 132 | 64 | 64 | 1 |
| monomial-8 | sparse/exhaustive | 1.718 | 2.121 | 256 | 192 | 64 | 64 | 1 |
| monomial-8 | packed/exhaustive | 1.351 | 1.713 | 256 | 192 | 64 | 64 | 1 |
| monomial-8 | sparse/frontier | 2.705 | 3.137 | 193 | 129 | 64 | 64 | 1 |
| monomial-8 | packed/frontier | 2.903 | 3.239 | 193 | 129 | 64 | 64 | 1 |
| block-linear-8 | sparse/exhaustive | 9.060 | 10.825 | 1024 | 1936 | 240 | 480 | 2 |
| block-linear-8 | packed/exhaustive | 6.113 | 7.988 | 1024 | 1936 | 240 | 480 | 2 |
| block-linear-8 | sparse/frontier | 16.522 | 18.470 | 1052 | 1284 | 240 | 480 | 2 |
| block-linear-8 | packed/frontier | 16.102 | 17.983 | 1052 | 1284 | 240 | 480 | 2 |
| planted-quadratic-8 | sparse/exhaustive | 64.878 | 70.778 | 1024 | 14978 | 234 | 1670 | 40 |
| planted-quadratic-8 | packed/exhaustive | 18.377 | 23.656 | 1024 | 14978 | 234 | 1670 | 40 |
| planted-quadratic-8 | sparse/frontier | 147.205 | 153.163 | 1431 | 30469 | 234 | 1527 | 48 |
| planted-quadratic-8 | packed/frontier | 43.285 | 49.529 | 1431 | 30469 | 234 | 1527 | 48 |
| dependent-contradictory-8 | sparse/exhaustive | 4.246 | 5.565 | 768 | 512 | 256 | 256 | 1 |
| dependent-contradictory-8 | packed/exhaustive | 3.715 | 5.187 | 768 | 512 | 256 | 256 | 1 |
| dependent-contradictory-8 | sparse/frontier | 11.331 | 12.714 | 1027 | 772 | 256 | 256 | 1 |
| dependent-contradictory-8 | packed/frontier | 10.851 | 12.251 | 1027 | 772 | 256 | 256 | 1 |
| holdout-chain-linear-8 | sparse/exhaustive | 20.713 | 22.811 | 1792 | 9578 | 254 | 508 | 2 |
| holdout-chain-linear-8 | packed/exhaustive | 11.364 | 13.320 | 1792 | 9578 | 254 | 508 | 2 |
| holdout-chain-linear-8 | sparse/frontier | 18.169 | 20.192 | 1206 | 2372 | 254 | 508 | 2 |
| holdout-chain-linear-8 | packed/frontier | 17.963 | 20.161 | 1206 | 2372 | 254 | 508 | 2 |
| holdout-cycle-quadratic-8 | sparse/exhaustive | 9.670 | 10.720 | 2048 | 1839 | 209 | 209 | 1 |
| holdout-cycle-quadratic-8 | packed/exhaustive | 7.211 | 8.248 | 2048 | 1839 | 209 | 209 | 1 |
| holdout-cycle-quadratic-8 | sparse/frontier | 8.930 | 10.140 | 760 | 551 | 209 | 209 | 1 |
| holdout-cycle-quadratic-8 | packed/frontier | 8.639 | 9.857 | 760 | 551 | 209 | 209 | 1 |
| holdout-sparse-quadratic-8 | sparse/exhaustive | 120.734 | 127.314 | 1024 | 24903 | 239 | 2137 | 52 |
| holdout-sparse-quadratic-8 | packed/exhaustive | 21.744 | 28.373 | 1024 | 24903 | 239 | 2137 | 52 |
| holdout-sparse-quadratic-8 | sparse/frontier | 163.605 | 170.993 | 1488 | 30308 | 239 | 1893 | 50 |
| holdout-sparse-quadratic-8 | packed/frontier | 48.885 | 56.808 | 1488 | 30308 | 239 | 1893 | 50 |
| holdout-dense-quadratic-8 | sparse/exhaustive | 284.205 | 299.208 | 1024 | 36825 | 239 | 3424 | 68 |
| holdout-dense-quadratic-8 | packed/exhaustive | 54.454 | 64.697 | 1024 | 36825 | 239 | 3424 | 68 |
| holdout-dense-quadratic-8 | sparse/frontier | 893.176 | 919.376 | 1559 | 53677 | 239 | 3942 | 76 |
| holdout-dense-quadratic-8 | packed/frontier | 183.160 | 210.735 | 1559 | 53677 | 239 | 3942 | 76 |

The inherited peak-row counter observes post-XOR rows and stored pivots, not every incoming row; retain the maximum input-row size separately. It is sufficient to witness the reported 6-to-50 growth but must not be described as a universal all-transient peak.

## Host, provenance and deviations

- CPU: Intel Xeon Platinum 8370C, virtualized Linux x86-64, Python 3.12. CPU 0 affinity succeeded in all samples; one virtual NUMA node exposed. Affinity is not exclusive core reservation.
- Node-0 memory binding returned Invalid argument. Physical NUMA placement and DDR generation are unknown. The runner requested a process summary with ps -eo; that command failed with fatal library error, lookup self, leaving an empty list. No process-isolation claim is supported. Load averages are saved before/after.
- Five A/A pairs per case precede seven alternating AB/BA rounds per schedule. Each untraced worker performs three fresh computations. Four separate traced workers per case measure Python allocations. Raw failures/timeouts: none. Host noise and memory-binding failure are retained, not hidden.
- Source, inputs and protocol were committed before measurement: aa56cbb51ff53d8df64d65bdebd8aaefeb594d51. PRE_RUN_SHA256SUMS verifies every frozen source/input and the preserved previous artifact. The source run has not been retuned.
- A post-run receipt-auditor import-order error was fixed before its successful execution; measured sources and inputs were unchanged. This was audit plumbing, not a rerun or a numerical-data correction.
- The prior artifact remains unchanged in prior/. This follow-up adds a backend by loading the identical closure engine with different row/space/reducer classes. The packed backend is a bounded reference, not a scalable solver.

## Reproduce and audit

Run from boolean_representation_audit with Python 3.12 (standard library only):

```sh
sha256sum -c PRE_RUN_SHA256SUMS
python3 -m unittest -v test_backends
python3 audit_results.py
python3 measure.py --output runs/my-new-run
```

The last command refuses an existing directory. The delivered result is runs/20260929/results.json. The 504 individual receipts are compressed losslessly as samples.jsonl.gz; each JSON line has filename and data fields. Inspect with python3 -m gzip -d on a copy, or gzip.open in Python. SAMPLE_MANIFEST.json records bytes and SHA-256. Individual sample paths in results.json refer to the filename fields in that archive. AUDIT_RECEIPT.json and TEST_RECEIPT.txt retain verification output. The CI job replays correctness and audits evidence; it does not impose timing thresholds on shared runners.

## Next bounded question

- [x] Preserve and replay the original negative control.
- [x] Compare exact sparse and packed representations with identical scheduling and outputs.
- [x] Record intermediate fill-in, interaction graph, allocations, process RSS and raw timings.
- [x] Preserve failed performance gates and host-control limitations.
- [ ] Compare a compact support-sharing representation with these two references on a newly frozen, small synthetic suite. Include low-width and fully coupled inputs and charge construction/conversion. Do not raise the current variable limit or integrate an application.
- [ ] Test preservation of narrow variable interactions before considering a chordal method. A min-fill diagnostic is not chordal elimination or a completeness proof.
- [ ] Obtain a quieter, verifiably controlled host before treating the observed timing reductions as hardware results.

Prior art remains relevant: [PolyBoRi](https://polybori.sourceforge.net/features.html) uses shared ZDD representations for Boolean polynomials; [Cifuentes–Parrilo](https://arxiv.org/abs/1604.02618) develops chordal networks for structured polynomial ideals. Their features/abstract were reviewed in this follow-up. These experiments do not implement either system or establish novelty.
