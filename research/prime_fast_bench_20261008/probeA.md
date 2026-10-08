# Prime-field ECDLP breakthrough probe

Host threads: 14.

## A. Preprocessing rho on prime-field curves (native engine)

`S = 2^s` table entries, `2^dp ≈ √(n/S)`, online walks single-threaded; plain rho = `rho_parallel` on the same targets.

| bits | S (2^s) | dp | precompute steps | precompute wall s | online steps mean (targets) | 2√(n/S) | plain rho steps √(πn/4) | per-target speed-up (steps) | online wall ms mean | plain rho wall s | break-even targets |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 32 | 2^11 | 10 | 2.225e6 | 0.0 | 5.372e3 (16/16) | 2.896e3 | 5.808e4 (measured 6.029e4) | 11× | 7.80 | 0.061 | 1 |
| 36 | 2^12 | 12 | 2.005e7 | 0.3 | 6.140e3 (16/16) | 8.192e3 | 2.323e5 (measured 1.702e5) | 38× | 8.58 | 0.160 | 2 |
| 40 | 2^13 | 14 | 3.046e8 | 10.3 | 2.195e4 (16/16) | 2.317e4 | 9.293e5 (measured 7.567e5) | 42× | 38.71 | 0.617 | 18 |
| 44 | 2^15 | 15 | 3.188e9 | 64.5 | 6.188e4 (16/16) | 4.634e4 | 3.717e6 (measured 5.994e6) | 60× | 161.56 | 0.719 | 116 |

JSON written to /private/tmp/claude-501/-Volumes-SSD990-crypto/b4dec92e-7670-423d-ac27-7fc98f39d636/scratchpad/probeA2.json
