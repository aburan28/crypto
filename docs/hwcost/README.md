# `hwcost`: calibrated analytical time-and-cost model

A hardware profile and a workload go in; a time and cost breakdown comes out.
Source: [`src/cryptanalysis/hwcost/`](../../src/cryptanalysis/hwcost/),
CLI [`src/bin/hwcost.rs`](../../src/bin/hwcost.rs). Validation study:
[`research/hwcost_model_20261016/`](../../research/hwcost_model_20261016/)
(calibrated model within 2–6% median error on the rho kernel's held-out
batch sizes on this repository's cloud host; the a priori model without a
fitted constant is within 30% on 4 of 5 fields and 54% off on `fp64`).
Round 2: [`research/hwcost_round2_20261016/`](../../research/hwcost_round2_20261016/)
(a native build, the `fp64` branch cost, a GPU profile within 18% a priori
on 46 ECC2K-130 builds, and IC against rho in time and dollars).

## The model

Per phase of a workload, on `M` machines with `c` cores each:

| term | formula | where the constants come from |
| --- | --- | --- |
| compute | `iterations × Σ_op (count_op × cycles_op) / clock` | profile `ops` (`op` at reciprocal throughput, `op@latency` at dependent-chain latency) |
| memory | `iterations × bytes_per_iteration / bandwidth(level)`, level = first that holds `working_set_bytes_per_core` | profile `memory` |
| roofline | phase time = `max(compute, memory)` | |
| pipes | with `pipes_overlap`, cycles per iteration = the busiest pipe's `Σ count × cycles` (ops tagged with a `pipe`); otherwise the sum over ops | profile `ops[].pipe` |
| parallel | `T₁ · ((1 − f) + f / (M·c · efficiency(c)))` | profile `scaling` (measured; interpolated, flagged beyond the last point) |
| transfers | `messages × latency + bytes / bandwidth`; `overlaps_compute` → `max` with compute (summed per channel), otherwise added | profile `channels` |
| cost | `wall hours × M × $/h + egress + storage` | profile `price` |

Every constant has a `source`. Sources beginning `assumed` or
`illustrative` are listed under "rests on unmeasured constants" in every
prediction. An operation the workload names but the profile does not price
is an error.

## Commands

```sh
cargo build --release --bin hwcost

# Calibrate this host (≈ 3 min; --quick for a smoke test).
./target/release/hwcost measure --out runs/myhost
./target/release/hwcost fit --run runs/myhost      # profile.json, validation.{json,csv}

# Predict any workload.
./target/release/hwcost predict --hardware runs/myhost/profile.json \
    --workload docs/hwcost/examples/memory-bound-scan.json [--machines 10] [--cores 2] [--json]

# Pollard rho with distinguished points, from the kernel's counted ops.
./target/release/hwcost rho --hardware runs/myhost/profile.json \
    --group-log2 64,112,128 --machines 1,1000 [--negation-map] [--csv] [--write-workloads DIR]

# Two builds or variants measured in alternation (A1 B1 A2 B2): per-arm ratios,
# drift and a class; the kernel ratio of two fields inside the same runs.
./target/release/hwcost compare --base runs/a1,runs/a2 --cand runs/b1,runs/b2 --out runs/compare
./target/release/hwcost pair --runs runs/r1,runs/r2 --a fp64 --b fp64cm --min-batch 4096 \
    --threshold 0.85 --id F1prime --out runs/pair
./target/release/hwcost measure --fields fp64,fp64cm --out runs/r1   # a subset of fields
./target/release/hwcost measure --fields pf64 --prime-fast-rho-bits 44 --out runs/r2  # + prime_fast::rho_parallel
./target/release/hwcost fit --run runs/b1 --residual-from runs/a1   # carry fitted residuals across builds

# GPU profiles (SMs as cores, overlapping pipes) checked against recorded builds.
./target/release/hwcost gpu-check --calibration ecc2k130/benchmarks/roofline/calibration.json \
    --fit-arm base --reference-gpu rtx-pro-6000 --out runs/gpu

# IC and rho from frozen ecbench records: wall time, §5a band, dollars, projections.
./target/release/hwcost ic-price --band BAND.json --sessions SESSIONS_DIR --ic-arm ic-compact \
    --rho-arm rho-frob --tier 27,40 --targets 48,64,128 --out runs/ic
```

Fields `fp64bf` and `fp64cm` are `Fp64` with the reduction's and the
subtraction's selections written as masks and as
`std::hint::select_unpredictable`; results are bit for bit `Fp64`'s. Only
`fp64cm` changes the machine code (round 2, F3), and production `Fp64` has
used the same selections since
[`research/fp64_branchfree_20261016/`](../../research/fp64_branchfree_20261016/).
Field `pf64` times `prime_fast::Fp64` (`p < 2⁶²`, the production prime-field
rho's arithmetic) at 2⁶² − 57, and `--prime-fast-rho-bits` adds an arm
`pfrho.<bits>`: `prime_fast::rho_parallel` end to end, ns per walk step.
`hwcost rho` prices groups on `fp64` and `w1`–`w4` only.

## Profile (`hwcost.hardware.v1`)

```json
{
  "schema": "hwcost.hardware.v1",
  "name": "…", "kind": "cpu", "clock_ghz": 2.878, "clock_source": "measured: …",
  "cores": 4, "isa": ["avx2", "adx", "…"],
  "ops": { "w4.mul": { "cycles": 67.8, "latency_cycles": 82.2, "source": "measured: …" } },
  "memory": [ { "name": "L1d", "capacity_bytes": 49152, "shared": false,
                "read_gbps_one_core": 71.3, "read_gbps_all_cores": 269.4, "source": "measured: …" },
              { "name": "DRAM", "capacity_bytes": null, "shared": true, "…": "…" } ],
  "scaling": [ { "cores": 1, "efficiency": 1.0 }, { "cores": 4, "efficiency": 0.903 } ],
  "scaling_source": "measured: …",
  "channels": { "network": { "latency_s": 0.001, "bandwidth_bytes_per_s": 1.25e9, "source": "assumed: …" } },
  "price": { "usd_per_hour": 0.1785, "egress_usd_per_gb": 0.09, "storage_usd_per_gb_month": 0.023, "source": "illustrative: …" }
}
```

`hwcost fit` writes one from a measurement run. To price operations from
instruction tables instead (uops.info, Agner Fog, llvm-mca or uiCA on the
compiled hot loop), add entries to `ops` by hand and say so in `source`
(`"table: uops.info, MULX r64 r64 r64, Golden Cove"`). A GPU fits the same
schema: SMs as `cores`, instruction classes per SM-cycle as `ops`.

## Workload (`hwcost.workload.v1`)

Phases (`iterations` as a number or `{"log2": x}`, `ops_per_iteration`,
`bytes_per_iteration`, `working_set_bytes_per_core`, `parallel_fraction`),
`transfers` (`channel`, `bytes`, `messages`, `overlaps_compute`,
`per_machine_links`), `storage` (`bytes`, `months`, or `null` for the
predicted wall time) and `egress_bytes`. Examples:
[`examples/memory-bound-scan.json`](examples/memory-bound-scan.json)
(memory-bound phase, Amdahl serial fraction, blocking disk load) and
[`examples/rho-w2-n112-B1024-m1000.json`](examples/rho-w2-n112-B1024-m1000.json)
(written by `hwcost rho --write-workloads`).

```text
$ hwcost predict --hardware research/hwcost_model_20261016/runs/final/profile.json \
      --workload docs/hwcost/examples/memory-bound-scan.json
  phase addition scan: 2^33.00 iterations, 1.6 cycles/iter (0.54 ns compute, 0.89 ns DRAM memory per iter, memory-bound), efficiency 0.903
  transfer load from disk over disk: 2^36.0 bytes, 10.25 min
  wall time 10.29 min
  cost $1.61 = compute $0.03 (1.715e-1 machine-hours) + egress $0.00 + storage $1.58
  rests on unmeasured constants:
    - disk channel: assumed: block-storage baseline 125 MB/s, 1 ms per request (not measured)
    - price: illustrative: AWS us-east-1 on-demand c7i.xlarge (4 vCPU) $0.1785/h, …
```

## Limits

- A profile is one host, one build. The study's profile is an x86-64-baseline
  build (no BMI2/ADX/AVX2 code generation) on a shared 4-vCPU cloud guest
  (isolation L0); a `-C target-cpu=native` build needs its own `measure` run:
  round 2 measured its `fp64` kernel 1.6–2.6× faster, and fitted residuals do
  not carry across builds for `fp64`.
- The roofline is a `max()`: it misses the few-percent additive cost of a
  working set beyond the last-level cache that the study measured.
- Machines are assumed to scale perfectly with each other; coordination is
  priced only through the transfers you declare.
