## C3: epsilon at kappa_R2 (R1 primitives, W-walk collection), ratio to os = 0

| arm | sizes | gm ratio over sizes |
|---|---|---|
| lpsgebw2 | 4 | 1.0000 |
| lpsgebw2o100c@16 | 4 | 0.9401 |
| lpsgebw2o35c@16 | 4 | 0.9447 |
| lpsgebw2o70c@16 | 4 | 0.9368 |
| lpsgebw2o70c@32 | 4 | 0.9331 |

C3 rule: epsilon* = `lpsgebw2o70c@32` (ratio 0.9331)

## C4: C1 + C2 + C3, kappa re-tuned: quadratic in ln kappa per size, per epsilon

| arm | log2 q | grid points | kappa* | predicted min cost | at kappa_R2 |
|---|---|---|---|---|---|
| lpsgebw2+C1C2 | 16 | 4 | 1.6051 | 2.4783e8 | 2.4834e8 |
| lpsgebw2+C1C2 | 18 | 4 | 1.2611 | 1.5420e9 | 1.6081e9 |
| lpsgebw2+C1C2 | 20 | 4 | 1.0447 | 1.0494e10 | 1.0852e10 |
| lpsgebw2+C1C2 | 21 | 4 | 0.9914 | 2.7244e10 | 2.7795e10 |
| lpsgebw2o100c@16+C1C2 | 16 | 4 | 3.2347 | 1.9583e8 | 2.1978e8 |
| lpsgebw2o100c@16+C1C2 | 18 | 4 | 1.8951 | 1.3546e9 | 1.3720e9 |
| lpsgebw2o100c@16+C1C2 | 20 | 4 | 1.3818 | 9.7285e9 | 9.7524e9 |
| lpsgebw2o100c@16+C1C2 | 21 | 4 | 1.2368 | 2.5709e10 | 2.5711e10 |
| lpsgebw2o35c@16+C1C2 | 16 | 4 | 3.0030 | 2.1738e8 | 2.2377e8 |
| lpsgebw2o35c@16+C1C2 | 18 | 4 | 1.4532 | 1.4379e9 | 1.4416e9 |
| lpsgebw2o35c@16+C1C2 | 20 | 4 | 1.1827 | 1.0029e10 | 1.0062e10 |
| lpsgebw2o35c@16+C1C2 | 21 | 4 | 1.1017 | 2.6235e10 | 2.6287e10 |
| lpsgebw2o70c@16+C1C2 | 16 | 4 | 3.2347 | 2.0163e8 | 2.1876e8 |
| lpsgebw2o70c@16+C1C2 | 18 | 4 | 1.6955 | 1.3904e9 | 1.3902e9 |
| lpsgebw2o70c@16+C1C2 | 20 | 4 | 1.2924 | 9.8615e9 | 9.8446e9 |
| lpsgebw2o70c@16+C1C2 | 21 | 4 | 1.1818 | 2.5879e10 | 2.5885e10 |
| lpsgebw2o70c@32+C1C2 | 16 | 4 | 3.2347 | 1.9959e8 | 2.1743e8 |
| lpsgebw2o70c@32+C1C2 | 18 | 4 | 1.7163 | 1.3811e9 | 1.3818e9 |
| lpsgebw2o70c@32+C1C2 | 20 | 4 | 1.3032 | 9.8168e9 | 9.8029e9 |
| lpsgebw2o70c@32+C1C2 | 21 | 4 | 1.1885 | 2.5774e10 | 2.5774e10 |

C4 rule: arm `lpsgebw2o100c@16+C1C2`; ln kappa_C4(b) = 0.3478 + (-0.1916) (b - 20)

| b | 12 | 14 | 16 | 18 | 20 | 21 | 22 |
|---|---|---|---|---|---|---|---|
| kappa_C4(b) | 6.555 | 4.469 | 3.047 | 2.077 | 1.416 | 1.169 | 0.965 |
| kappa_R2(b) | 2.772 | 2.297 | 1.903 | 1.576 | 1.306 | 1.189 | 1.082 |
