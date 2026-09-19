# AutoLab: ECC2K-130 Crossbred α

Local control plane for the Crossbred thread in
[`RESEARCH_ECC2K130_ROUTE_TARGETS.md`](../../RESEARCH_ECC2K130_ROUTE_TARGETS.md).
Harbor is not required. Public synthetic Koblitz systems only.

**Incumbent, frozen before this harness.** X1 and X3 each have two usable
`m = 3` rungs and no `α` fit. GPU search stays off (`filters = 0`). The
product-law floor is untouched.

## Quick start

```bash
python3 research/ecc2k130_crossbred_autolab_20260918/crossbred_autolab.py plan
python3 research/ecc2k130_crossbred_autolab_20260918/crossbred_autolab.py preflight
python3 research/ecc2k130_crossbred_autolab_20260918/test_crossbred_autolab.py

python3 research/ecc2k130_crossbred_autolab_20260918/crossbred_autolab.py \
  launch --beat smoke.x1_n5
python3 research/ecc2k130_crossbred_autolab_20260918/crossbred_autolab.py \
  launch --beat replay.x3_k1_n7
python3 research/ecc2k130_crossbred_autolab_20260918/crossbred_autolab.py \
  launch --beat fit.alpha
python3 research/ecc2k130_crossbred_autolab_20260918/crossbred_autolab.py \
  launch --beat x5.ffd_chained_sym_m4
python3 research/ecc2k130_crossbred_autolab_20260918/crossbred_autolab.py promote --all
```

## Beats

| id | what | promotion |
|---|---|---|
| `smoke.x1_n5` | Replay chained `K_0` `n=5` `m=3` | matches frozen `Q/C` |
| `replay.x3_k1_n7` | Replay symmetrised `K_1` `n=7` `m=3` | matches frozen `Q/C` |
| `fit.alpha` | OLS **per frame** on frozen usable rungs | ≥4 rungs in one frame and `α ≤ 1.5` |
| `x5.ffd_chained_m4` | Chained-`x` FFD at `m=3,4`, `n=9,15`, 4 draws | FFD that grows with `n` at `m=4` |
| `x5.ffd_chained_m4_16` | Same, 16 draws | **done** 2026-09-18: FFD max=3, does not grow |
| `x5.ffd_chained_sym_m4_smoke` | Chained symmetrised `S₃` at `m=4`, `n=7,9,15`, 4 draws, `d_max=4` | FFD that grows with `n` on `ell>1` |
| `x5.ffd_chained_sym_m4` | Same, 16 draws | **done** 2026-09-18: FFD max=4 at `n=9,15`; `n=7` `ell=1` is a collapse, not H1 |

A two-point sketch is not a fit. Changing T4, `Q_enum`, or the divisor
after seeing a cell is inadmissible.

## Layout

```
protocol.json              # beats, boundaries, frozen pointers
crossbred_autolab.py       # plan / preflight / launch / status / verify / promote
test_crossbred_autolab.py
instruction.md             # Harbor-shaped agent brief
task.toml
runs/                      # gitignored scratch
evidence/                  # promoted receipts
```
