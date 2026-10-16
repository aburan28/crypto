# `hwcost`: calibrated analytical time-and-cost model, validated on the rho iteration (2026-10-16)

A native model of how long a calculation takes on given hardware and what it
costs: counted operations × measured cycles per operation, a memory roofline
against the level that holds the working set, Amdahl with measured scaling,
latency + bytes/bandwidth for transfers, and machine-hours × price for cost
(tool docs: [`docs/hwcost/README.md`](../../docs/hwcost/README.md)). This
study calibrates it on one host and tests it against a Pollard rho iteration
whose operation counts are counted from the code.

- Protocol, frozen before the run at `a9a457ec1`: [`PROTOCOL.md`](PROTOCOL.md)
- Run: [`runs/final/`](runs/final/) (`measurements.json`, `samples.csv`,
  `medians.csv`, `profile.json`, `validation.{json,csv}`, `fit.txt`,
  `rho.csv`, `rho_negation.csv`)
- Run report (descriptive): [`reports/run-20261016T0459Z/run_report.pdf`](reports/run-20261016T0459Z/run_report.pdf)
- Results report: [`report.pdf`](report.pdf) (Typst source `report.typ`;
  tables read the run records directly; graphs from [`plot.py`](plot.py),
  which only lays out recorded values)

## Results (run `runs/final`, 7 repetitions per arm, isolation L0)

Host: Intel Xeon @ 2.10 GHz (family 6 model 207), 4 vCPUs under KVM, measured
clock 2.878 GHz (reg-reg `add` chain; `imul` chain 2.797 GHz). Build: x86-64
baseline (no BMI2/ADX/AVX2 code generation).

| claim | statement | value | threshold | verdict |
| --- | --- | --: | --: | :-: |
| C0 | add-chain and imul-chain clocks agree | 2.8% | ≤ 5% | PASS |
| V1.fp64 | calibrated model C, test batches (B ≥ 4096), median \|error\| | 2.7% | ≤ 30% | PASS |
| V1.w1 | 〃 | 5.5% | ≤ 30% | PASS |
| V1.w2 | 〃 | 3.7% | ≤ 30% | PASS |
| V1.w3 | 〃 | 1.8% | ≤ 30% | PASS |
| V1.w4 | 〃 | 2.8% | ≤ 30% | PASS |
| V2.fp64 | a priori model A (no fitted parameter), all batches | 54.2% | ≤ 30% | **FAIL** |
| V2.w1 | 〃 | 3.0% | ≤ 30% | PASS |
| V2.w2 | 〃 | 18.1% | ≤ 30% | PASS |
| V2.w3 | 〃 | 12.4% | ≤ 30% | PASS |
| V2.w4 | 〃 | 13.6% | ≤ 30% | PASS |
| V3.w4 | size-extrapolated model S (op costs fitted on 1–3 limbs), all batches | 7.4% | ≤ 30% | PASS |

- One fitted constant per field (the per-iteration residual, fitted on
  B ≤ 1024) predicts B = 4096 … 4,194,304 within 2–6% median error.
- Without it, counted operations × microbenchmarked costs explain 80–103% of
  the `FpWide` kernels' time but only about half of `fp64`'s. The `fp64`
  kernel is ~1.6× slower than `w1`'s although its operations are cheaper in
  isolation. Unverified lead: data-dependent branches in `Fp64::sub`/`redc`
  (no `perf` on this host).
- The roofline puts every arm on the compute side (the kernel draws ≤ 3.5 GB/s
  against ≥ 8.8 GB/s single-thread DRAM bandwidth); growing the working set
  past the last-level cache still adds 2–10%, which a `max()` roofline omits.
- Thread scaling at 4 threads: 0.83–0.94 by field (profile uses the median,
  0.903).

Application (descriptive, extrapolated; illustrative c7i.xlarge price):
parallel rho on a 2^112 prime-order group takes about 21.6 days and $92.5k on
1,000 four-core machines; on P-256 (2^256) about 2^69.4 years on the same
fleet. Every row lists the unmeasured constants it rests on (network, price,
cross-machine scaling).

Not an index-calculus result (no gap-to-rho row); no ECDLP methods are ranked,
so the §5a unpriced-work band does not apply.

## Reproduce

```sh
cargo build --release --bin hwcost
D=research/hwcost_model_20261016/runs/final
./target/release/hwcost measure --out $D 2> $D/stderr.txt
./target/release/hwcost fit --run $D --fit-max-batch 1024 > $D/fit.txt
./target/release/hwcost rho --hardware $D/profile.json \
  --group-log2 48,64,80,96,112,128,160,192,224,256 --machines 1,1000 --csv > $D/rho.csv
./target/release/hwcost rho --hardware $D/profile.json --negation-map \
  --group-log2 48,64,80,96,112,128,160,192,224,256 --machines 1,1000 --csv > $D/rho_negation.csv
python3 plot.py runs/final reports/run-20261016T0459Z   # matplotlib; layout only
# PDFs: typst compile --root . report.typ; typst compile --root . reports/run-20261016T0459Z/run_report.typ
```

## Next

1. Measure a `-C target-cpu=native` build (MULX/ADX, AVX-512 IFMA `Fp52`)
   and compare profiles.
2. Check the `fp64` branch-misprediction lead with hardware counters.
3. A GPU profile from `benchmarks/hardware-limits/` and the ECC2K-130
   per-pipe instruction counts, validated against the recorded 15.1 B it/s.
4. Price an index-calculus pipeline phase by phase from ecbench's counted
   units, so IC and rho are compared in time and dollars as well as in GAE.
