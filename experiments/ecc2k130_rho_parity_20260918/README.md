# ECC2K-130 toy-suite rho parity receipt

Cited from this directory. Protocol:
[`RESEARCH_ECC2K130_RHO_PARITY.md`](../../RESEARCH_ECC2K130_RHO_PARITY.md).

Host `ip-172-31-19-103`, 2026-09-18. Unit: exclusive group operations.
Class: engineering. Not n=131.

| iteration | what | all-cases gate |
|---|---|---|
| 0 | folded-row rule + walked collection, two-pass table build | unmet: mean α 1.05–1.24, max 1.31. `iteration-0/` |
| 1 | same, one-pass table build | **met**: mean α 0.76–0.91, max 0.918. `iteration-1/` |

`table.txt` is the runner's per-pair print. `summary.json` is the
machine receipt. `run.log` is the raw cargo session including compiler
warnings; it is not a cost claim.
