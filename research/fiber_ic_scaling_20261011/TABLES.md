## `scaling_h4`

| curve | bits | h | log₂ #E | blind-b128 | closed-b128 | closed-b16 | closed-b256 | closed-b32 | closed-b64 | rho-neg | envelope arm | envelope ÷ rho | blind ÷ envelope |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|:--|--:|--:|
| `icv1-fp18-t334-97871af7` | 18 | 4 | 17.58 | 80.4 | 303.7 | 16.1 | 1196.0 | 25.6 | 80.3 | 1.9 | closed-b16 | 8.43 | 4.98 |
| `icv1-fp18-tm204-54b2fadf` | 18 | 4 | 17.83 | 73.7 | 278.2 | 15.5 | 1095.7 | 25.1 | 73.6 | 1.8 | closed-b16 | 8.83 | 4.76 |
| `icv1-fp22-tm2826-0afa05b6` | 22 | 4 | 21.08 | 28.7 | 91.9 | 21.6 | 355.8 | 15.0 | 28.5 | 1.4 | closed-b32 | 11.04 | 1.91 |
| `icv1-fp22-tm1214-827be133` | 22 | 4 | 21.34 | 27.3 | 84.3 | 21.5 | 325.4 | 16.4 | 25.4 | 1.5 | closed-b32 | 11.26 | 1.67 |
| `icv1-fp26-t374-04f2b312` | 26 | 4 | 25.37 | 26.9 | 27.9 | 65.9 | 83.6 | 27.7 | 19.3 | 0.9 | closed-b64 | 20.60 | 1.40 |
| `icv1-fp26-tm7664-a6fbaf53` | 26 | 4 | 25.80 | 36.3 | 25.8 | 79.3 | 72.7 | 39.4 | 25.8 | 1.4 | closed-b64 | 17.85 | 1.40 |
| `icv1-fp30-t8770-57966dac` | 30 | 4 | 29.34 | 90.9 | 33.5 | 273.8 | 32.9 | 131.6 | 57.8 | 0.9 | closed-b256 | 35.10 | 2.76 |
| `icv1-fp30-tm29402-75bd2650` | 30 | 4 | 29.47 | 79.2 | 37.6 | 386.2 | 34.7 | 162.0 | 69.6 | 0.9 | closed-b256 | 36.99 | 2.29 |

Envelope ÷ rho against size, h = 4: slope of log₂(ratio) on log₂ #E = **0.177** (95 % interval 0.145 to 0.209; prediction 1/6 = 0.167), over 8 curves. Intercept -0.14: the fitted ratio would reach 1 at log₂ #E = 0.8.

## `scaling_h72`

| curve | bits | h | log₂ #E | blind-b128 | closed-b128 | closed-b16 | closed-b256 | closed-b32 | closed-b64 | rho-neg | envelope arm | envelope ÷ rho | blind ÷ envelope |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|:--|--:|--:|
| `icv1-fp18-tm22-a5ba8c3f` | 18 | 72 | 17.56 | 369.9 | 22967.6 | 393.4 | 91468.1 | 1482.4 | 5797.3 | 2.7 | closed-b16 | 146.94 | 0.94 |
| `icv1-fp18-t872-49a6dd87` | 18 | 72 | 17.94 | 296.3 | 20134.3 | 344.1 | 80186.5 | 1299.0 | 5081.6 | 3.8 | closed-b16 | 91.04 | 0.86 |
| `icv1-fp22-tm2412-15c27732` | 22 | 72 | 21.43 | 101.7 | 6001.5 | 111.0 | 23897.2 | 390.3 | 1516.3 | 2.3 | closed-b16 | 48.10 | 0.92 |
| `icv1-fp22-tm1504-967cbbde` | 22 | 72 | 21.48 | 104.6 | 5894.6 | 108.8 | 23469.8 | 384.0 | 1489.8 | 2.2 | closed-b16 | 49.74 | 0.96 |
| `icv1-fp26-t8250-b50c04e1` | 26 | 72 | 25.47 | 66.6 | 1482.6 | 47.7 | 5895.5 | 107.3 | 377.7 | 1.1 | closed-b16 | 42.54 | 1.40 |
| `icv1-fp26-tm1038-caffe525` | 26 | 72 | 25.88 | 58.1 | 1287.1 | 49.6 | 5116.5 | 94.2 | 330.0 | 1.3 | closed-b16 | 39.07 | 1.17 |
| `icv1-fp30-tm48354-87c391fd` | 30 | 72 | 29.61 | 141.3 | 359.8 | 82.9 | 1407.4 | 64.9 | 106.0 | 1.1 | closed-b32 | 60.36 | 2.18 |
| `icv1-fp30-tm45748-31aed113` | 30 | 72 | 29.92 | 181.2 | 325.8 | 105.9 | 1264.1 | 60.4 | 100.9 | 1.0 | closed-b32 | 63.06 | 3.00 |

Envelope ÷ rho against size, h = 72: slope of log₂(ratio) on log₂ #E = **-0.073** (95 % interval -0.183 to 0.038; prediction 1/6 = 0.167), over 8 curves. Intercept 7.66: the fitted ratio would reach 1 at log₂ #E = NaN.

## `scaling_h72_small_b`

| curve | bits | h | log₂ #E | closed-b4 | closed-b8 | rho-neg | envelope arm | envelope ÷ rho | blind ÷ envelope |
|:--|--:|--:|--:|--:|--:|--:|:--|--:|--:|
| `icv1-fp18-tm22-a5ba8c3f` | 18 | 72 | 17.56 | 53.9 | 118.9 | 2.7 | closed-b4 | 20.12 | – |
| `icv1-fp18-t872-49a6dd87` | 18 | 72 | 17.94 | 50.8 | 104.9 | 3.8 | closed-b4 | 13.44 | – |
| `icv1-fp22-tm2412-15c27732` | 22 | 72 | 21.43 | 40.8 | 45.6 | 2.3 | closed-b4 | 17.67 | – |
| `icv1-fp22-tm1504-967cbbde` | 22 | 72 | 21.48 | 46.9 | 46.3 | 2.2 | closed-b8 | 21.18 | – |
| `icv1-fp26-t8250-b50c04e1` | 26 | 72 | 25.47 | 122.8 | 57.2 | 1.1 | closed-b8 | 51.02 | – |
| `icv1-fp26-tm1038-caffe525` | 26 | 72 | 25.88 | 152.6 | 47.8 | 1.3 | closed-b8 | 37.63 | – |
| `icv1-fp30-tm48354-87c391fd` | 30 | 72 | 29.61 | 412.0 | 191.9 | 1.1 | closed-b8 | 178.44 | – |
| `icv1-fp30-tm45748-31aed113` | 30 | 72 | 29.92 | 585.6 | 238.8 | 1.0 | closed-b8 | 249.46 | – |

Envelope ÷ rho against size, h = 72: slope of log₂(ratio) on log₂ #E = **0.305** (95 % interval 0.184 to 0.425; prediction 1/6 = 0.167), over 8 curves. Intercept -1.83: the fitted ratio would reach 1 at log₂ #E = 6.0.

## `scaling_h72_combined`

| curve | bits | h | log₂ #E | blind-b128 | closed-b128 | closed-b16 | closed-b256 | closed-b32 | closed-b4 | closed-b64 | closed-b8 | rho-neg | envelope arm | envelope ÷ rho | blind ÷ envelope |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|:--|--:|--:|
| `icv1-fp18-tm22-a5ba8c3f` | 18 | 72 | 17.56 | 369.9 | 22967.6 | 393.4 | 91468.1 | 1482.4 | 53.9 | 5797.3 | 118.9 | 2.7 | closed-b4 | 20.12 | 6.87 |
| `icv1-fp18-t872-49a6dd87` | 18 | 72 | 17.94 | 296.3 | 20134.3 | 344.1 | 80186.5 | 1299.0 | 50.8 | 5081.6 | 104.9 | 3.8 | closed-b4 | 13.44 | 5.83 |
| `icv1-fp22-tm2412-15c27732` | 22 | 72 | 21.43 | 101.7 | 6001.5 | 111.0 | 23897.2 | 390.3 | 40.8 | 1516.3 | 45.6 | 2.3 | closed-b4 | 17.67 | 2.50 |
| `icv1-fp22-tm1504-967cbbde` | 22 | 72 | 21.48 | 104.6 | 5894.6 | 108.8 | 23469.8 | 384.0 | 46.9 | 1489.8 | 46.3 | 2.2 | closed-b8 | 21.18 | 2.26 |
| `icv1-fp26-t8250-b50c04e1` | 26 | 72 | 25.47 | 66.6 | 1482.6 | 47.7 | 5895.5 | 107.3 | 122.8 | 377.7 | 57.2 | 1.1 | closed-b16 | 42.54 | 1.40 |
| `icv1-fp26-tm1038-caffe525` | 26 | 72 | 25.88 | 58.1 | 1287.1 | 49.6 | 5116.5 | 94.2 | 152.6 | 330.0 | 47.8 | 1.3 | closed-b8 | 37.63 | 1.21 |
| `icv1-fp30-tm48354-87c391fd` | 30 | 72 | 29.61 | 141.3 | 359.8 | 82.9 | 1407.4 | 64.9 | 412.0 | 106.0 | 191.9 | 1.1 | closed-b32 | 60.36 | 2.18 |
| `icv1-fp30-tm45748-31aed113` | 30 | 72 | 29.92 | 181.2 | 325.8 | 105.9 | 1264.1 | 60.4 | 585.6 | 100.9 | 238.8 | 1.0 | closed-b32 | 63.06 | 3.00 |

Envelope ÷ rho against size, h = 72: slope of log₂(ratio) on log₂ #E = **0.168** (95 % interval 0.115 to 0.221; prediction 1/6 = 0.167), over 8 curves. Intercept 0.92: the fitted ratio would reach 1 at log₂ #E = -5.5.

## Cofactor 72 (`scaling_h72`) against cofactor 4 at matched size

| bits | mean envelope ÷ rho, h = 4 | h = 72 | ratio |
|--:|--:|--:|--:|
| 18 | 8.63 | 118.99 | 13.79 |
| 22 | 11.15 | 48.92 | 4.39 |
| 26 | 19.22 | 40.81 | 2.12 |
| 30 | 36.04 | 61.71 | 1.71 |

## Cofactor 72 (`scaling_h72_combined`) against cofactor 4 at matched size

| bits | mean envelope ÷ rho, h = 4 | h = 72 | ratio |
|--:|--:|--:|--:|
| 18 | 8.63 | 16.78 | 1.94 |
| 22 | 11.15 | 19.42 | 1.74 |
| 26 | 19.22 | 40.09 | 2.09 |
| 30 | 36.04 | 61.71 | 1.71 |

