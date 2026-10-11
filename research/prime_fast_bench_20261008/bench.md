# Prime-field ECDLP speed-up bench

Host: "macos" "aarch64" with 14 hardware threads.

## 1. Field and point micro-benchmarks (ns per operation)

| p bits | op | BigUint `FieldElement`/`Point` | `Fp64`/`FastCurve` | speed-up |
|---|---|---:|---:|---:|
| 28 | mul | 1198.2 | 7.85 | 153× |
| 28 | inv | 123428 | 23569.2 | 5× |
| 28 | inv, batched W=1024 | — | 75.79 | 1629× vs BigUint inv |
| 61 | mul | 2817.1 | 7.85 | 359× |
| 61 | inv | 918145 | 525.7 | 1747× |
| 61 | inv, batched W=1024 | — | 26.62 | 34496× vs BigUint inv |
| 28-bit ladder | affine add (own inverse) | 207110 | 341 | 608× |
| 28-bit ladder | affine add, batched W=1024 | — | 899.3 | 230× vs BigUint add |
| 56-bit generated | affine add (own inverse) | 875698 | 14349 | 61× |
| 56-bit generated | affine add, batched W=1024 | — | 106.0 | 8260× vs BigUint add |

## 2. Pollard rho: BigUint Floyd (`pollard_rho_ecdlp`) vs native parallel DP (`rho_parallel`)

Known-answer targets `Q = [k]G`, `k` from seed; medians over 3 seed(s).  Baseline uses the ladder examples' step cap `2^(bits/2+8)`.

| bits | curve | baseline wall s | fast wall s | speed-up | fast steps/s (total) | threads × walkers | dp bits | fast steps (median) |
|---|---|---:|---:|---:|---:|---|---|---:|
| 16 | a3-bench-16bit-p256class | 0.424 | 0.0373 | 11× | 5.11e3 | 1 × 32 | 0 | 1.920e2 |
| 20 | a3-bench-20bit-p256class | 3.214 | 0.3685 | 9× | 2.26e3 | 1 × 32 | 1 | 8.640e2 |
| 24 | a3-bench-24bit-p256class | 22.135 | 0.4280 | 52× | 1.47e4 | 1 × 56 | 3 | 5.096e3 |
| 28 | a3-bench-28bit-p256class | 147.866 | 0.4144 | 357× | 4.15e4 | 3 × 75 | 3 | 1.230e4 |
<!-- generated a3-fast-32bit-p256class: p=4294967291 b=31 n=4294959973 G=(2, 2120051858) in 0.02s -->
| 32 | a3-fast-32bit-p256class | — | 0.3989 | — | 8.52e4 | 14 × 64 | 3 | 2.611e4 |
<!-- generated a3-fast-36bit-p256class: p=68719476731 b=33 n=68719022539 G=(2, 51741851039) in 0.38s -->
| 36 | a3-fast-36bit-p256class | — | 0.6901 | — | 5.60e5 | 14 × 259 | 3 | 4.634e5 |
<!-- generated a3-fast-40bit-p256class: p=1099511627563 b=56 n=1099512551701 G=(1, 15267301451) in 1.50s -->
| 40 | a3-fast-40bit-p256class | — | 0.9139 | — | 1.34e6 | 14 × 1024 | 3 | 9.411e5 |
<!-- generated a3-fast-44bit-p256class: p=17592186044399 b=95 n=17592192476281 G=(1, 8379642404149) in 5.81s -->
| 44 | a3-fast-44bit-p256class | — | 1.5978 | — | 2.87e6 | 14 × 1024 | 5 | 4.689e6 |
<!-- generated a3-fast-48bit-p256class: p=281474976710591 b=9 n=281474991897067 G=(1, 111984877719828) in 1.17s -->
| 48 | a3-fast-48bit-p256class | — | 2.7366 | — | 7.19e6 | 14 × 1024 | 7 | 2.536e7 |
<!-- generated a3-fast-52bit-p256class: p=4503599627370323 b=124 n=4503599513826641 G=(2, 1220001265989541) in 43.20s -->
| 52 | a3-fast-52bit-p256class | — | 9.0461 | — | 1.07e7 | 14 × 1024 | 9 | 9.652e7 |
<!-- generated a3-fast-56bit-p256class: p=72057594037927931 b=60 n=72057594202405003 G=(2, 795920743654909) in 54.20s -->
| 56 | a3-fast-56bit-p256class | — | 19.9999 | — | 1.28e7 | 14 × 1024 | 11 | 2.570e8 |

## 3. Index calculus, 2-decomposition: Semaev S₃ sweep (BigUint) vs direct subtraction (native)

Factor base policy `m = 120` smallest-x points, `+6` extra relations (the a=−3 ladder policy).  Both drivers verify `[x]G = Q`.

| bits | baseline total s (fb / rel / LA) | fast total s (fb / rel / LA) | speed-up | fast targets | fast rows (1/2) | rho fast s |
|---|---|---|---:|---:|---|---:|
| 16 | 14.50 (0.006 / 14.48 / 0.011) | 0.004 (0.0001 / 0.004 / 0.0001) | 3481× | 298 | 1/125 | 0.0174 |
| 20 | 351.69 (0.006 / 351.50 / 0.137) | 0.499 (0.0001 / 0.499 / 0.0002) | 704× | 4788 | 3/123 | 0.0193 |
| 24 | — | 6.899 (0.0001 / 6.816 / 0.0825) | — | 64988 | 3/123 | 0.1806 |
| 28 | — | 137.345 (0.0001 / 137.344 / 0.0001) | — | 1230190 | 3/123 | 0.5407 |
| 32 | — | 1744.111 (0.0022 / 1744.109 / 0.0002) | — | 17689251 | 1/125 | 0.5102 |

## 4. Index calculus, 3-decomposition: S₄ quartic roots (BigUint) vs pair table (native)

Baseline cost is measured per relation with `find_one_relation_s4_counted` (fb = 120); the fast driver runs end to end with the same base and `+6` extra relations.

| bits | baseline s/relation (relations timed) | fast total s (table / rel / LA) | fast s/relation | per-relation speed-up | pair table entries | fast targets | rows (1/2/3) |
|---|---:|---|---:|---:|---:|---:|---|
| 16 | 42.153 (2) | 0.002 (0.001 / 0.001 / 0.0003) | 0.000005 | 9055901× | 14400 | 4 | 0/4/122 |
| 20 | — | 0.009 (0.007 / 0.002 / 0.0003) | 0.000017 | — | 14400 | 54 | 0/4/122 |
| 24 | — | 0.033 (0.001 / 0.032 / 0.0002) | 0.000250 | — | 14400 | 979 | 0/5/121 |
| 28 | — | 1.252 (0.001 / 1.251 / 0.0002) | 0.009925 | — | 14400 | 14800 | 0/4/122 |
| 32 | — | 20.519 (0.001 / 20.518 / 0.0005) | 0.162840 | — | 14400 | 223193 | 0/2/124 |

JSON written to /private/tmp/claude-501/-Volumes-SSD990-crypto/b4dec92e-7670-423d-ac27-7fc98f39d636/scratchpad/bench.json
