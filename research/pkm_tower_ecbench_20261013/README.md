# Petit–Kosters–Messeng's prime-field index calculus in ecbench, with the root-finding oracles, against rho

2026-10-13.  Protocol: [`PROTOCOL.md`](PROTOCOL.md) (registered before the
runs).  Factor base: `prime-kummer:t=T` (`ic_framework::prime_tower`), PKM's
Kummer tower `V = μ_{2^t}`, the zeros of `x^{2^t} − 1`.  Oracles: the
root-finders of [`../pdp_root_finder_20261012/`](../pdp_root_finder_20261012/README.md)
as `ic.pipeline` plug-ins.  Companion binary ladder and the digest PDF:
[`../pdp_root_finder_ecbench_20261013/`](../pdp_root_finder_ecbench_20261013/README.md).
**Class: engineering.**

## Result

Index calculus on prime curves with PKM's tower base, run end to end in
ecbench (factor base, pair table, relations, linear algebra, verification,
every group operation charged), costs **29–141× rho with the negation map**
on five curves of 20 to 28 bits.  All 896 runs (672 measured, 224 warm-ups) verified and every
replay audit is `OK`.  The ratio grows with the curve, by
`0.26 ± 0.28` bits of `S_IC/S_rho` per bit of `r` (least squares over the
five sizes, ± 1.96 s.e.), consistent with the `p^{1/6}` the floors predict
(0.17 per bit) and with an interval too wide at these sizes to separate the
two from zero.  The tower base costs the same as a plain base of the same
size, within 1 %: PKM's structure helps an algebraic solver, not a table.

## Curves

From `examples/prime_tower_curves.rs --bits 20,22,24,26,28 --two-adicity 12
--seed 20261013 --max-cofactor 4` (`curves.jsonl`), registered in
`docs/curves/registry.json`:

| bits | `icv1` slug | `p` | cofactor | `log₂ r` | `t` | base points | columns |
|--:|:--|--:|--:|--:|--:|--:|--:|
| 20 | `icv1-fp20-tm728-d215a2a1` | 946,177 | 2 | 18.85 | 8 | 246 | 123 |
| 22 | `icv1-fp22-tm2947-e15fed69` | 2,387,969 | 1 | 21.19 | 9 | 502 | 251 |
| 24 | `icv1-fp24-t2705-27cbb01b` | 10,940,417 | 1 | 23.38 | 9 | 490 | 245 |
| 26 | `icv1-fp26-t2838-8e779ca8` | 40,857,601 | 4 | 23.28 | 10 | 980 | 490 |
| 28 | `icv1-fp28-t4155-50f492b1` | 231,657,473 | 3 | 26.20 | 11 | 2,012 | 1,006 |

## Tables

Mean `S` per arm (`ecbench table`; eight one-target workloads, three
measured rounds; floor `√(π/4)` = 0.886):

| bits | rho (negation) | IC + `mitm` | IC + `rootfind-mitm` | IC + `rootfind-mitm`, plain base | IC + `rootfind-pairs` |
|--:|--:|--:|--:|--:|--:|
| 20 | 1.154 | 32.40 | 39.26 | 39.02 | 1,522 |
| 22 | 1.210 | 44.59 | 50.74 | 50.70 | 1,388 |
| 24 | 1.145 | 31.13 | 33.64 | 33.25 | 3,885 |
| 26 | 0.997 | 91.29 | 96.32 | 95.81 | — |
| 28 | 0.936 | 128.3 | 131.6 | 131.7 | — |

Paired ratios (`ecbench compare`, 95 % cluster-bootstrap intervals; wall is
the harness's admitted median):

| bits | `rootfind` ÷ rho, `S` | ÷ rho, wall | tower ÷ plain base, `S` | `rootfind` ÷ `mitm`, `S` | ÷ `mitm`, wall |
|--:|--:|--:|--:|--:|--:|
| 20 | 34.0 [27.9, 42.6] | 15.5 | 1.006 [0.982, 1.032] | 1.212 [1.200, 1.225] | 1.31 |
| 22 | 41.9 [34.3, 53.4] | 26.4 | 1.001 [0.996, 1.006] | 1.138 [1.133, 1.142] | 1.13 |
| 24 | 29.4 [23.7, 37.3] | 22.8 | 1.012 [0.989, 1.035] | 1.081 [1.076, 1.086] | 1.17 |
| 26 | 96.6 [77.8, 124.1] | 83.4 | 1.005 [0.994, 1.020] | 1.055 [1.050, 1.065] | 1.06 |
| 28 | 140.6 [113.0, 182.3] | 178.6 | 0.999 [0.994, 1.004] | 1.026 [1.025, 1.026] | 1.09 |

The A/A control (`rho-neg-aa`) is 1.000 in `S` on every size.  Hit rates
are 0.82–1.00 at the chosen `t` (relations per trial), so the relation phase
is `trials × 2A` additions plus the table.  Median solve wall time grows from
3.6 ms to 292 ms for the root-finder against 0.25 ms to 1.45 ms for rho.

**On prime curves the new oracle is slightly slower than `mitm`** (3–31 % in
wall, 3–21 % in `S`): the prime group has no batched addition (each
addition inverts), so the block of 32 abscissae only adds the overshoot on
hits.  The binary group's shared inversion is what made it faster there.

## Decision

- P1 (every run verified, audits identical): **held** (896 runs; 144, 144,
  144, 120, 120 replays).
- P2 (`S_IC/S_rho > 1`, growing about as `p^{1/6}`): **held** in sign and
  size on every curve; the growth fit (0.26 ± 0.28 bits per bit of `r`) is
  consistent with 0.17 but its interval includes zero over this range.
- P3 (tower and plain base equal within 20 %): **held**, within 1.2 %.
- P4 (`rootfind` below `mitm` in wall): **failed** on prime curves (above,
  by 3–31 %); within 1.3× in `S`: held.
- Class: **engineering**.  PKM's question — is there anything better than
  Gröbner bases for System (4) — has a measured yes (10⁵–10¹² fewer
  operations per decomposition, `../pdp_root_finder_20261012/`), and end to
  end the index calculus it serves stays above rho, as the pair-table floor
  `p^{2/3}` against rho's `p^{1/2}` says it must.

No Gröbner arm runs here: `f4_fp_tower` is not an `ic.pipeline` oracle; its
per-decomposition cost on PKM's own cells is frozen in the PKM rounds and
compared in `../pdp_root_finder_20261012/`.

## Environment

Claude Code cloud container: Intel Xeon @ 2.10 GHz, 4 vCPU, 15 GiB RAM,
rustc 1.97.0; `ecbench run --cpus auto`, one thread.  Binary: built from
commit `0ad5bc38b` (the tower base and the registry entries) in a separate
target directory, so the binary ladder's sessions kept the binary that
measured them.  Conductor's control plane was unreachable: no task card.

## Files

`PROTOCOL.md`; `curves.jsonl` and `curves.log`; `specs/SPEC-p{20..28}.json`;
`run_sessions.sh`; `sessions/p*/` (sealed); `audits/p*.audit.json`;
`tables/`; `figures/` (`plot.py`, `s_by_size`, `phase_split`,
`wall_by_size` as SVG and PNG).
