# `two_sylow_search` on `fp_backend::Fp64`, before and after the branch-free selections (2026-10-17)

PR #2034 ([`research/fp64_branchfree_20261016/`](../fp64_branchfree_20261016/))
made `fp_backend::Fp64`'s reduction, `add` and `sub` branch-free. That sped
up the rho kernel 2.3–2.5×, left dependent `mul`/`sqr` chains unchanged,
and made dependent `sub` chains with a predictable wrap 18–33% slower.
This study times a latency-bound production user of `fp_backend::Fp64`:
`two_sylow_search`'s halving search (`Halver::run_draw`).

- Protocol, frozen before the first run at `e71198313`:
  [`PROTOCOL.md`](PROTOCOL.md).
- Runs: [`runs/`](runs/). These are A1 B1 A2 B2, with `hwcost measure --fields fp64,w1 --two-sylow`.
- Run report:
  [`reports/run-20261017T0018Z/run_report.pdf`](reports/run-20261017T0018Z/run_report.pdf).
- Results report: [`report.pdf`](report.pdf), Typst source `report.typ`.
  Graphs come from [`plot.py`](plot.py), layout only.
- Benchmarks copy:
  [`benchmarks/20261017-two-sylow-fp64/`](../../benchmarks/20261017-two-sylow-fp64/).

## Result

Host: Intel Xeon @ 2.10 GHz (family 6 model 207), 4 vCPUs under KVM,
isolation L0, no PMU. The prime is 2^63 − 25, the largest prime ≡ 3 mod 4
below 2^63. The search covers both families and all three search arms;
each sample is wall ns per draw.

| claim | statement | value | threshold | verdict |
| --- | --- | --: | --: | :-: |
| S1 | no `two_sylow_search` arm slower, max(r1, r2) over 6 arms | 1.083 | ≤ 1.10 | PASS |
| K1 | positive control: `fp64` rho kernel, B ≥ 4096 | 0.442 | ≤ 0.85 | PASS |
| C1 | negative control: `FpWide<1>` kernel, B ≥ 4096, max \|r − 1\| | 0.045 | ≤ 0.10 | PASS |
| G | both binaries compute the same search (gate checksums identical in all 4 runs) | — | — | PASS |

**The change is neutral for this search** within this host's resolution:
- **First pair:** ratios 0.98–1.02.
- **Second pair:** ratios 1.02–1.08, a noisier pair in which unchanged
  controls also read high.
- **Classes:** five arms show no detected effect; one is classed slower
  (`z2xz2.fast-chi`: 1.023, 1.083).

A draw is dominated by `pow` (dependent `mul`/`sqr` chains, unchanged) and
`jacobi_u64` (untouched). Production P-192 runs use `FpP192`, which the
change does not touch.

## Reproduce

Follow the protocol's "Commands" block, then:

```sh
cd research/two_sylow_fp64_20261017
python3 plot.py runs reports/run-20261017T0018Z   # matplotlib; layout only
typst compile --root . report.typ
typst compile --root . reports/run-20261017T0018Z/run_report.typ
```
