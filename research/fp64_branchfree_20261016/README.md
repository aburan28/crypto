# Branch-free selections in production `Fp64` (`fp_backend` and `prime_fast`), 2026-10-16

`hwcost` round 2 ([`research/hwcost_round2_20261016/`](../hwcost_round2_20261016/),
addendum 1) found that `Fp64`'s data-dependent selections cost the rho
kernel about half its time on random operands. This study ports the fix to
both production one-word Montgomery fields:

- [`fp_backend::Fp64`](../../src/cryptanalysis/fp_backend.rs), with `p < 2^63`;
- [`prime_fast::Fp64`](../../src/cryptanalysis/prime_fast.rs), with
  `p < 2^62`, the arithmetic of the production prime-field rho
  `rho_parallel`.

In each, the Montgomery reduction's final subtraction, `add`'s final
subtraction and `sub`'s wrap-around become `std::hint::select_unpredictable`.
Results are bit for bit the previous ones: new unit tests compare them with
the old branchy formulas on edge values and random inputs.

- Protocol, frozen before the first run at `41a7a585c`
  (`runs/protocol_commit.txt`): [`PROTOCOL.md`](PROTOCOL.md).
- Runs: [`runs/`](runs/). These are `base-1`, `cand-1`, `base-2` and
  `cand-2` (A1 B1 A2 B2), plus `compare/`, `cm-pair/` and `disasm/`.
- Run report (descriptive):
  [`reports/run-20261016T2319Z/run_report.pdf`](reports/run-20261016T2319Z/run_report.pdf).
- Results report: [`report.pdf`](report.pdf), Typst source `report.typ`;
  its tables read the run records. Graphs come from [`plot.py`](plot.py),
  which only lays out recorded values.
- Benchmarks copy:
  [`benchmarks/20261016-fp64-branchfree/`](../../benchmarks/20261016-fp64-branchfree/).

This is not an index-calculus improvement and adds no gap-to-rho row.
Counted costs do not change; only wall time does.

## Results

Host: Intel Xeon @ 2.10 GHz (family 6 model 207), 4 vCPUs under KVM,
isolation L0, no PMU. Builds: baseline-ISA `--release`, the production
configuration. The baseline binary is the protocol commit with the two field
files restored to `main` (`f25d5e520`). Seven repetitions per arm.

| claim | statement | value | threshold | verdict |
| --- | --- | --: | --: | :-: |
| P1 | `fp_backend::Fp64` rho kernel, every B ≥ 4096, max(r1, r2) | 0.436 | ≤ 0.85 | PASS |
| P1pf | `prime_fast::Fp64` rho kernel, every B ≥ 4096, max(r1, r2) | 0.466 | ≤ 0.85 | PASS |
| P4 | production `prime_fast::rho_parallel` (44-bit curve, one thread), ns per step, max(r1, r2) | 0.721 | ≤ 0.95 | PASS |
| P2a–d | dependent mul and sqr chains, both fields, max(r1, r2) | ≤ 1.018 | ≤ 1.10 | PASS |
| P3 | `fp64 / Fp64Cm` kernel in the candidate runs, B ≥ 4096 | 1.027 | ≤ 1.10 | PASS |
| C1a, C1b | unchanged controls (`Fp64Cm`, `FpWide<1>` kernels), max \|r − 1\| | 0.025, 0.030 | ≤ 0.10 | PASS |

- **Rho kernel.** It runs 2.3–2.5× faster on `fp_backend::Fp64` (ratios
  0.40–0.44 at B ≥ 4096) and 2.1–2.2× faster on `prime_fast::Fp64`
  (0.45–0.47). The candidate's kernel step compiles to the study's `Fp64Cm`
  code: 303 → 285 instructions, 32 → 20 conditional jumps, 4 → 12 `cmov`.
- **Production rho.** `prime_fast::rho_parallel` drops from 98.9 / 98.6 to
  68.2 / 71.1 ns per step, 1.39–1.45× the steps per second. Its worker loop
  goes from 119 to 107 conditional jumps and from 18 to 30 `cmov`.
- **Cost model.** The a priori op-count model, which missed both baseline
  kernels by 51–55%, predicts both candidates within 6–9%.
- **Stated cost.** A dependent chain of subtractions with a fixed
  subtrahend runs 18–33% slower on `fp64` and 26–27% slower on `pf64`: a
  select replacing a branch that predicted well. The same arm moved −27%
  and −14% for the unchanged controls, so part of that move is code
  placement. Rho does not run such chains.

## Reproduce

Follow the protocol's "Commands" block, then:

```sh
cd research/fp64_branchfree_20261016
python3 plot.py runs reports/run-20261016T2319Z   # matplotlib; layout only
typst compile --root . report.typ
typst compile --root . reports/run-20261016T2319Z/run_report.typ
```

## Next

1. Time a latency-bound user of `fp_backend::Fp64`, such as `two_sylow_search`'s
   square-root chains or the F4 engine at an `Fp64` prime, to confirm that no
   caller relies on a predictable `sub` chain.
2. Repeat on a host with a PMU to count the removed branch misses directly.
