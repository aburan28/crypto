# Prime-field ECDLP speed-up bench

Host: "macos" "aarch64" with 14 hardware threads.

## 1. Field and point micro-benchmarks (ns per operation)

| p bits | op | BigUint `FieldElement`/`Point` | `Fp64`/`FastCurve` | speed-up |
|---|---|---:|---:|---:|
| 28 | mul | 133.7 | 10.23 | 13× |
| 28 | inv | 19613 | 1921.8 | 10× |
| 28 | inv, batched W=1024 | — | 223.96 | 88× vs BigUint inv |
| 61 | mul | 568.1 | 9.96 | 57× |
| 61 | inv | 177098 | 626.9 | 283× |
| 61 | inv, batched W=1024 | — | 16.88 | 10494× vs BigUint inv |
| 28-bit ladder | affine add (own inverse) | 25923 | 464 | 56× |
| 28-bit ladder | affine add, batched W=1024 | — | 132.7 | 195× vs BigUint add |
| 56-bit generated | affine add (own inverse) | 99863 | 1004 | 100× |
| 56-bit generated | affine add, batched W=1024 | — | 64.1 | 1558× vs BigUint add |

## 2. Pollard rho: BigUint Floyd (`pollard_rho_ecdlp`) vs native parallel DP (`rho_parallel`)

Known-answer targets `Q = [k]G`, `k` from seed; medians over 3 seed(s).  Baseline uses the ladder examples' step cap `2^(bits/2+8)`.

| bits | curve | baseline wall s | fast wall s | speed-up | fast steps/s (total) | threads × walkers | dp bits | fast steps (median) |
|---|---|---:|---:|---:|---:|---|---|---:|
| 16 | a3-bench-16bit-p256class | 0.073 | 0.0289 | 3× | 6.66e3 | 1 × 32 | 0 | 1.920e2 |
| 20 | a3-bench-20bit-p256class | 0.369 | 0.0341 | 11× | 2.63e4 | 1 × 32 | 1 | 8.640e2 |
| 24 | a3-bench-24bit-p256class | 6.045 | 0.1007 | 60× | 4.89e4 | 1 × 56 | 3 | 5.096e3 |
| 28 | a3-bench-28bit-p256class | 22.659 | 0.0486 | 466× | 1.60e5 | 3 × 75 | 3 | 2.318e4 |
<!-- generated a3-fast-32bit-p256class: p=4294967291 b=31 n=4294959973 G=(2, 2120051858) in 0.07s -->
| 32 | a3-fast-32bit-p256class | — | 0.0968 | — | 4.07e5 | 14 × 64 | 3 | 5.005e4 |
<!-- generated a3-fast-36bit-p256class: p=68719476731 b=33 n=68719022539 G=(2, 51741851039) in 0.08s -->
| 36 | a3-fast-36bit-p256class | — | 0.1400 | — | 1.96e6 | 14 × 259 | 3 | 3.553e5 |
<!-- generated a3-fast-40bit-p256class: p=1099511627563 b=56 n=1099512551701 G=(1, 15267301451) in 0.45s -->
| 40 | a3-fast-40bit-p256class | — | 0.2443 | — | 1.36e6 | 14 × 1024 | 3 | 3.328e5 |
<!-- generated a3-fast-44bit-p256class: p=17592186044399 b=95 n=17592192476281 G=(1, 8379642404149) in 1.88s -->
| 44 | a3-fast-44bit-p256class | — | 1.3449 | — | 3.16e6 | 14 × 1024 | 5 | 4.243e6 |
<!-- generated a3-fast-48bit-p256class: p=281474976710591 b=9 n=281474991897067 G=(1, 111984877719828) in 0.38s -->
| 48 | a3-fast-48bit-p256class | — | 1.7238 | — | 1.16e7 | 14 × 1024 | 7 | 2.046e7 |

JSON written to /private/tmp/claude-501/-Volumes-SSD990-crypto/b4dec92e-7670-423d-ac27-7fc98f39d636/scratchpad/bench2.json
