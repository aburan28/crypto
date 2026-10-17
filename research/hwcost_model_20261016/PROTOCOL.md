# Protocol: calibrating and validating `hwcost` on the Pollard rho iteration

Frozen before the measured run. The commit that adds this file is the
protocol commit; the run report names it.

## Question

Does an analytical model of the form

```
time per iteration = max( Σ_op count_op × cost_op ,  bytes_per_iteration / bandwidth(level) )
```

with operation counts taken exactly from the code and per-operation costs
measured by microbenchmarks, predict the measured time of a Pollard rho
iteration on this host, and to what error? Is one fitted per-iteration
constant, fitted on small instances, enough to bring the error inside 30% on
larger ones? Does fitting per-operation costs against limb count on smaller
fields extrapolate to a larger one?

## Boundary

- **Measured:** one implementation of the rho iteration (`hwcost::kernel`:
  r-adding walk, 32 precomputed points, `B` affine walks per core, one
  Montgomery batch inversion per step) over the crate's `fp_backend` prime
  fields, on this host.
- **Not measured:** no discrete logarithm is solved; the walk is not
  restarted at distinguished points; no index-calculus pipeline is run.
  This is tooling. It is not an index-calculus improvement (AGENTS.md
  success criterion) and has no gap-to-rho row. No two ECDLP methods are
  ranked, so the §5a unpriced-work band does not apply.
- **Host:** the session's cloud container: `Intel Xeon Processor @ 2.10GHz`
  (family 6 model 207), 4 vCPUs under KVM, 1 socket, 15 GiB. Isolation L0
  (shared host, no core pinning, no frequency control). Cache sizes are read
  from sysfs as the guest reports them.
- **Build:** `cargo build --release --bin hwcost` with the default target
  (x86-64 baseline; the binary records the target features it was compiled
  with). An ISA-specific build (`-C target-cpu=native`: BMI2/ADX/AVX-512) is
  a different measurement, not this one.

## Arms (all interleaved across repetitions, order rotated per repetition)

Parameters: `reps = 7`, `target_ms = 100`, seed `20261016`.

| group | arms |
| --- | --- |
| M0 clock | dependent chain of 32 reg-reg `add` (1 cycle each); dependent chain of 32 `imul` (3 cycles each) as cross-check |
| M1 field ops | per field: `mul`, `sqr`, `add` (a subtraction), `inv` (`inv(x+1)` chain), each as throughput (8 independent chains; inv 4) and latency (1 chain) |
| M2 memory | sequential read of 16 KiB, 32 KiB, 128 KiB, 512 KiB, 1, 4, 16, 64, 256 MiB and 1 GiB per thread, on 1 thread and on all 4 |
| M3 kernel | per field × `B ∈ {16, 64, 256, 1024, 4096, 16384, 65536, 262144, 2^20, 2^22}`, one thread, ns per iteration |
| M4 scaling | per field, `B = 1024` per thread, on 1, 2 and 4 threads |
| M5 counts | one counted step of every M3 arm on the very walks that are timed |

Fields (`w1`–`w4` are one implementation, `FpWide<N>` Montgomery CIOS, with
primes filling every limb; `fp64` is the dedicated one-word backend):

| field | backend | prime |
| --- | --- | --- |
| fp64 | `Fp64` | 2^63 − 25 |
| w1 | `FpWide<1>` | 2^64 − 59 |
| w2 | `FpWide<2>` | 2^128 − 159 |
| w3 | `FpWide<3>` | P-192 |
| w4 | `FpWide<4>` | P-256 |

Every field's batched walk is checked against an unbatched reference (one
inversion per walk) and every point against the curve `y² = x³ + 7` before
any timing; the run aborts on a mismatch or on any degenerate step
(`x_P = x_R`).

## Models (all single-core, ns per iteration)

- **A, a priori, no fitted parameter.** Counted operations from M5 priced
  at M1 throughput, except the batch inversion and the two multiplications
  per iteration on its prefix-product chain (`Arith::mul_dep`), priced at
  M1 latency; under the M2 roofline: bytes moved per iteration
  `12 × element bytes` against the one-thread bandwidth of the first level
  holding the working set `4 × element bytes × B + table`.
- **L, all-latency (descriptive only).** Every operation at M1 latency.
- **C, calibrated, one parameter per field.** A plus a constant residual
  per iteration, the mean of `measured − A` over the **fit set
  `B ≤ 1024`**; evaluated on the **test set `B ≥ 4096`**.
- **S, size-extrapolated (w4 only).** Every M1 cost, throughput and latency,
  fitted against the limb count `n` on `w1, w2, w3` (`α + γn²` for `mul`,
  `sqr`, `inv`; `α + βn` for `add` and for the C residual), extrapolated to
  `n = 4`, then used as in C.

## Claims (thresholds frozen here)

| id | statement | threshold |
| --- | --- | --- |
| C0 | the M0 add-chain and imul-chain clock estimates agree, \|relative difference\| | ≤ 0.05 |
| V1.f | model C, test set, median \|error\|, for each field f | ≤ 0.30 |
| V2.f | model A, all batch sizes, median \|error\|, for each field f | ≤ 0.30 |
| V3.w4 | model S on w4, all batch sizes, median \|error\| | ≤ 0.30 |

Error is `predicted / measured − 1`; the measured value is the median over
the 7 repetitions, with quartiles reported. A FAIL is reported as a FAIL; it
is not refitted away. A V1 failure on a field falsifies "within about 30%
once calibrated" for that field on this kernel.

Reported whatever the verdicts: every row's measured value and quartiles,
all four predictions and their errors, the memory level and bound the
model assigns, the per-operation costs, the fitted residuals, the size fits
and the measured thread-scaling efficiency.

## Application (descriptive; no claim)

With the fitted profile, `hwcost rho` predicts wall time and cost of
parallel rho with distinguished points (expected `√(πn/2)` iterations, or
`√(πn/4)` with the negation map, plus the distinguished-point tail) for
group orders `2^48 … 2^256` on one 4-core machine and on 1,000 of them. The
price is illustrative (AWS us-east-1 on-demand c7i.xlarge, $0.1785/h, a
4-vCPU instance of the same CPU generation, not this container); network
and disk constants are assumed. The output names every unmeasured constant
it rests on.

## Commands

```sh
cargo build --release --bin hwcost
D=research/hwcost_model_20261016/runs/final
./target/release/hwcost measure --out $D 2> $D/stderr.txt
./target/release/hwcost fit --run $D --fit-max-batch 1024 > $D/fit.txt
./target/release/hwcost rho --hardware $D/profile.json \
  --group-log2 48,64,80,96,112,128,160,192,224,256 --machines 1,1000 --csv > $D/rho.csv
./target/release/hwcost rho --hardware $D/profile.json --negation-map \
  --group-log2 48,64,80,96,112,128,160,192,224,256 --machines 1,1000 --csv > $D/rho_negation.csv
```
