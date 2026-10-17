# `hwcost` round 2: ISA effect, the `fp64` branch lead, a GPU profile, IC in time and dollars (2026-10-16)

Round 1 ([`research/hwcost_model_20261016/`](../hwcost_model_20261016/))
built and validated the analytical time-and-cost model on the rho iteration.
This round follows its four next steps: a `-C target-cpu=native` build (A),
the `fp64` branch-misprediction lead (B), a GPU profile validated on the
recorded ECC2K-130 builds (C), and index calculus priced in time and dollars
next to rho (D).

- Protocol, frozen at `597502b4d` before any run; addendum 1 frozen at
  `4e5994fc9` before its runs: [`PROTOCOL.md`](PROTOCOL.md)
- Runs: [`runs/`](runs/) (`isa-base-{1,2}`, `isa-native-{1,2}`, `compare/`,
  `cm-{1,2}`, `cm-pair/`, `f3/`, `gpu/`, `ic-koblitz/`, `ic-prime/`)
- Run report (descriptive): [`reports/run-20261016T1654Z/run_report.pdf`](reports/run-20261016T1654Z/run_report.pdf)
- Results report: [`report.pdf`](report.pdf) (Typst source `report.typ`;
  tables read the run records; graphs from [`plot.py`](plot.py), layout only;
  F3 disassembly counts from [`f3_disasm.sh`](f3_disasm.sh))
- Benchmarks copy: [`benchmarks/20261016-hwcost-round2/`](../../benchmarks/20261016-hwcost-round2/)

**Not an index-calculus improvement and no gap-to-rho row**: nothing here
changes an IC pipeline. Part D ranks IC and rho only where the wall-time
order and the order at every §5a weight (0, 0.1, 1 GAE per unpriced native
operation) agree.

## Results

Host for A and B: Intel Xeon @ 2.10 GHz (family 6 model 207), 4 vCPUs under
KVM, measured clock 2.66–2.83 GHz, isolation L0, **no PMU** (the guest kernel
has no perf driver for this CPU, so branch misses could not be counted).

| claim | statement | value | threshold | verdict |
| --- | --- | --: | --: | :-: |
| I1 | native vs baseline, every arm, both pairs of runs | kernel: 22 faster, 1 slower, 37 no effect | descriptive | — |
| I2 | calibrated model on the native build, every field (V1) | ≤ 0.120 | ≤ 0.30 | PASS |
| I3 | baseline residuals carried to the native build | `fp64` 1.51 / 1.38, `fp64bf` 1.32 / 1.29; `w1`–`w4` ≤ 0.094 | ≤ 0.30 | **FAIL** (`fp64`, `fp64bf`) |
| F0 | `fp64bf`, `fp64cm` walks equal `Fp64`'s point for point | all runs | gate | PASS |
| F1 | `fp64bf` ≤ 0.85 × `fp64` at every B ≥ 4096 | 1.050 | ≤ 0.85 | **FAIL**, invalid test (F3) |
| F2 | a priori model on `fp64bf` (V2) | 0.530 / 0.551 | ≤ 0.30 | **FAIL** |
| F3 | baseline disassembly, `step_fp64` / `step_fp64bf` / `step_fp64cm` | 303 / 303 / 285 instructions; 32 / 32 / 20 jumps; 4 / 4 / 12 `cmov` | descriptive | — |
| F1′ | `fp64cm` ≤ 0.85 × `fp64` at every B ≥ 4096 (addendum 1) | 0.500 | ≤ 0.85 | PASS |
| F2′ | a priori model on `fp64cm` (V2) | 0.079 / 0.043 | ≤ 0.30 | PASS |
| G1 | GPU, a priori, 46 builds, median \|error\| | 0.179 (no build with η > 1) | ≤ 0.30 | PASS |
| G2 | one η per GPU, 44 other builds | 0.0145 | ≤ 0.30 (≤ 0.10) | PASS (PASS) |
| G3 | RTX PRO 6000 η carried to 25 B200 builds | 0.195 | ≤ 0.30 | PASS |
| D1 | measured IC/rho wall ratio inside the band | Koblitz 3 / 8, prime 0 / 6 | descriptive | — |
| D2 | projections to log2 r 48 … 128, §5a rank | IC > rho at every target | descriptive | — |

- **ISA.** The native build runs the `fp64` and `fp64bf` kernels 1.6–2.6×
  faster (native / baseline 0.38–0.63) and leaves the multi-limb `FpWide`
  kernels unchanged at 37 of 40 batch sizes; isolated field operations get no
  faster. Round 1's a priori failure on `fp64` disappears on the native build
  (V2 0.145 / 0.082).
- **Branches.** `Fp64Bf`'s masks compiled back to the same branchy code (F3),
  so F1 tested nothing. With `std::hint::select_unpredictable` (`Fp64Cm`, the
  same field and values bit for bit) the baseline kernel takes 0.40–0.50× of
  `fp64`'s time at B ≥ 4096, and the op-count model with no fitted parameter
  predicts it within 4–8%. Round 1's `fp64` anomaly is the data-dependent
  branches in `Fp64::sub`/`redc`, about half the kernel's time on random
  operands (causal evidence by substitution; branch misses not counted).
- **GPU.** SMs as cores, the busiest pipe binding: 17.9% a priori over 46
  builds, 1.45% with one efficiency per GPU, 19.5% across GPUs. The cost
  arithmetic gives $123,749 for ECC2K-130 at 14.64 B/s on the RTX PRO 6000
  ($3.031/h); the paper states $124k.
- **IC in time and dollars.** On every measured curve §5a ranks, IC costs
  more wall time than the matched rho (Koblitz 4 of 8 curves ranked, prime 5
  of 6; the rest are unranked because rho's ~3 ms per-solve floor puts IC's
  wall time below rho's on the smallest curves while every band weight ranks
  IC above). The prime band understates IC's wall cost at log2 r 32–40 (the
  measured ratio is 1.1–2.8× above the band's whole range). Projected to
  2^128 (extrapolated, post hoc GAE basis): Koblitz rho 1.4–2.4 × 10^3
  core-years, IC 4.5 × 10^5 – 5.3 × 10^7; prime rho 1.3–1.6 × 10^4, IC
  1.4–6.1 × 10^11.

Deviations, all in the reports: F1's arms compiled to the same code (addendum
1 re-tested with `fp64cm`); the pre-registered wall-time projection carries
rho's per-solve floor into its exponent (a post hoc GAE-based projection is
added and labelled); the per-curve §5a `rank` column was added after the
analysis.

## Reproduce

The protocol's two "Commands" blocks, then:

```sh
D=research/hwcost_round2_20261016/runs
research/hwcost_round2_20261016/f3_disasm.sh target/release/hwcost baseline $D/f3
research/hwcost_round2_20261016/f3_disasm.sh target-native/release/hwcost native $D/f3
cd research/hwcost_round2_20261016
python3 plot.py runs reports/run-20261016T1654Z       # matplotlib; layout only
typst compile --root . report.typ
typst compile --root . reports/run-20261016T1654Z/run_report.typ
```

## Next

1. Port the `select_unpredictable` selections to the production
   `fp_backend::Fp64` in its own PR, gated on the rho benchmarks against a
   frozen baseline (candidate: about 2× on the baseline-ISA `Fp64` kernel).
2. Re-run parts A and B on a host with a PMU to count branch misses and split
   the native speedup between `cmov` and `mulx`/vector code.
3. Price IC's memory-bound phases with the roofline, so the band's unpriced
   native work gets a measured weight per phase.
