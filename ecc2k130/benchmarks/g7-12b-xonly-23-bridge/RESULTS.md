# Result: full-generator x-only bridge rejected

The sparse bridge repairs the scalar-generation failure of the [2]/[3]
x-only walk, but both valid GPU implementations are decisively slower than
the selected G7 binary. The selected binary remains
`ecc2k130/build/ecc2k130-local-packed`, SHA-256
`8c4ed76a152cc135280a68fd3559c3bf71bef726ab73eae0e04014c2e27c6d02`.
The 12 billion complete scalar updates/s target remains unmet.

## Boundary and decision table

The algorithmic boundary stays `sqrt(n/262)` and every measured row has
generic-work ratio 1. Full-DLP `S` is unmeasured. Each timing sample performs
34,359,738,368 complete updates on the same AWS g7.2xlarge / RTX PRO 4500 at
165 W. Intervals are Student-t 95% intervals over three paired log ratios.

| Variant | B complete updates/s | Paired speedup | 95% CI | Rate / 12 B/s | Generic-work ratio | Correctness | Class / decision |
|:--|--:|--:|:--|--:|--:|:--|:--|
| Selected packed runtime, initial screen | 6.311324 | 1.000000 | [1.000000, 1.000000] | 0.525944 | 1 | passed | reference |
| Sparse bridge, initial | 5.248565 | 0.827384 | [0.814735, 0.840229] | 0.437380 | 1 | passed | engineering regression |
| Selected packed runtime, recompute screen | 6.310157 | 1.000000 | [1.000000, 1.000000] | 0.525846 | 1 | passed | reference |
| Sparse bridge, low-live recompute | 5.264523 | 0.829866 | [0.816504, 0.843447] | 0.438710 | 1 | passed | engineering regression |

The initial bridge used 80 registers, a 48-byte stack frame, 72 bytes of spill
stores, 60 bytes of spill loads and 28,032 shared bytes. Outlining the rare
bridge worsened spills to 80/64 bytes and was rejected before timing. The
preregistered low-live rewrite reduced the stack to 32 bytes and spills to
64/40 bytes, but its entire confidence interval remained below one.

## Correctness and persistence

The transition set is complete: `1+s^3` is a nonsquare, and its order together
with the order of 2 has lcm `ell-1`. Independent arithmetic established

```
x(P + sigma^3(P)) = A(x)^2 / (x B(x)^2)
A = x^7 + x^6 + x^4 + x^3 + x + 1
B = x^6 + x^5 + x^4 + x^3 + x^2 + x + 1.
```

Both timed implementations passed all 117 built-in checks. The initial build
CPU-replayed 8,896 dense distinguished endpoints and 272 bounded partial
endpoints with zero drops; the recompute build replayed another 272 partial
endpoints with zero drops. Both passed memcheck, initcheck and synccheck.
Checkpoint version 4 resumed bit-for-bit: uninterrupted and split runs produced
the same SHA-256
`837babe22813e58184b4381c90df4b9e590fee899631078e039bd321fa9f89e0`.
Versions 3 and 4 rejected one another in both directions with exit status 6.

## Attribution diagnostic

The preregistered non-promotable diagnostic measured the zero-spill [2]/[3]
x-only core at 6.290649 B/s against 6.334035 B/s selected, ratio 0.993150.
That binary cannot be selected because its transition multipliers generate
only an index-two scalar subgroup. The diagnostic shows that the x-only base
map offers no material speed headroom and attributes the valid candidates'
roughly 17% regression to the sparse bridge path and its compiler costs.

Frozen comparisons are [comparison.json](comparison.json),
[comparison-recompute.json](comparison-recompute.json), and
[diagnostic-xonly23.json](diagnostic-xonly23.json). All raw timing, telemetry,
replay, persistence, compiler and sanitizer logs remain in this directory.
