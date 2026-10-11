# Direct polynomial delta-3/4 result on G7

The candidate is a decisive regression and does not advance. The selected
eight-jump executable remains unchanged, and the 12 B/s goal remains open.

| Variant | Median B updates/s | Paired raw ratio | Raw 95% CI | Collision-adjusted ratio | Adjusted 95% CI | Rate / 12 B/s | Decision |
|---|---:|---:|---:|---:|---:|---:|---|
| selected eight-jump | 6.272357 | 1.000000 | [1.000000, 1.000000] | 1.000000 | [1.000000, 1.000000] | 0.522696 | reference |
| indexed two-jump | 6.289420 | 1.001295 | [0.992584, 1.010083] | 0.988367 | [0.979768, 0.997041] | 0.524118 | algorithmic regression |
| direct polynomial delta-3/4 | 4.363499 | 0.689846 | [0.676097, 0.703874] | 0.680939 | [0.667367, 0.694786] | 0.363625 | engineering regression |

Each of the nine retained timing samples performs 34,359,738,368 complete
scalar updates on the same AWS g7.2xlarge / RTX PRO 4500 at 165 W. Process
ownership observations found no foreign GPU process. Every run reports the
expected update count and zero drops. Paired intervals use Student t on three
log ratios. The two-jump rows divide raw rates by the frozen 1.0130803718
collision-work ratio from 6,000 matched planted DLPs.

The direct candidate preserves the validated two-jump map but forms
`x + x^(2^j)` and `y + y^(2^j)` through three or four polynomial squarings.
It passes 16,652 independent delta comparisons, exact full and ragged GPU
state/DP comparisons, CPU replay, bidirectional resume, and memcheck,
initcheck and synccheck. The timing kernel uses 80 registers, 16 local bytes
per thread and 28,032 shared bytes per block. The result shows that repeated
polynomial squaring costs substantially more than the removed normal-basis
conversion and routing networks on this GPU.

Generic work remains `sqrt(n/262)`. The selected work ratio is 1, full-DLP S
is null, and no total-DLP speedup is claimed.

[Frozen comparison](comparison.json) · [Validation audit](validation-audit.json) · [Protocol](README.md)
