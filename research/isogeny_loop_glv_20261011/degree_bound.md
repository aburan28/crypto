# Class-group-free lower bound on loop-endomorphism cost

KAPPA (M per bit of degree): projective 14.196, affine 5.047, xonly 3.786; saving = 8 * ceil(log2 r / 2) M.

| curve | family | log2 r | log2 d_min(D_K) | Delta = f^2 D_K verified | saving M | bound projective M (ratio) | bound affine M (ratio) | bound xonly M (ratio) | cheapest loop found, M (frozen; projective) |
|:--|:--|--:|--:|:--|--:|--:|--:|--:|--:|
| ECC2K-130 | koblitz | 130 | 1.0 | yes | 520 | 14 (0.03) | 5 (0.01) | 4 (0.01) | 2.0 |
| icv1-f2m83-tm6151469093347-debefd74 | koblitz | 82 | 1.0 | yes | 328 | 14 (0.04) | 5 (0.02) | 4 (0.01) | 2.0 |
| icv1-f2m83-t6151469093347-cdcc5432 | koblitz | 53 | 1.0 | yes | 216 | 14 (0.07) | 5 (0.02) | 4 (0.02) | 2.0 |
| ECC2K-95 | koblitz | 95 | 1.0 | yes | 384 | 14 (0.04) | 5 (0.01) | 4 (0.01) | 2.0 |
| ECC2K-108 | koblitz | 108 | 1.0 | yes | 432 | 14 (0.03) | 5 (0.01) | 4 (0.01) | 2.0 |
| secp256k1 | prime | 256 | 0.0 | yes | 1024 | 0 (0.00) | 0 (0.00) | 0 (0.00) | 52.5 |
| id-GostR3410-2001-CryptoPro-B-ParamSet | prime | 256 | 7.3 | yes | 1024 | 103 (0.10) | 37 (0.04) | 28 (0.03) | 127.5 |
| id-GostR3410-2001-CryptoPro-A-ParamSet | prime | 256 | 255.9 | yes | 1024 | 3632 (3.55) | 1292 (1.26) | 969 (0.95) | not computed |
| sect113r1 | binary | 113 | 112.4 | yes | 456 | 1595 (3.50) | 567 (1.24) | 425 (0.93) | 9594.5 |
| ECC2-131 | binary | 131 | 129.1 | yes | 528 | 1833 (3.47) | 652 (1.23) | 489 (0.93) | 21449.5 |
| sect131r1 | binary | 131 | 129.7 | yes | 528 | 1841 (3.49) | 655 (1.24) | 491 (0.93) | 82224.5 |
| sect131r2 | binary | 131 | 130.6 | yes | 528 | 1854 (3.51) | 659 (1.25) | 494 (0.94) | 19806.0 |
| ECCp-131 | prime | 131 | 130.0 | yes | 528 | 1845 (3.49) | 656 (1.24) | 492 (0.93) | 49072.5 |
| ECC2-109 | binary | 109 | 103.4 | yes | 440 | 1468 (3.34) | 522 (1.19) | 391 (0.89) | 15836.5 |
| ECCp-109 | prime | 109 | 108.4 | yes | 440 | 1538 (3.50) | 547 (1.24) | 410 (0.93) | 21075.0 |
| ECC2-97 | binary | 97 | 91.2 | yes | 392 | 1295 (3.30) | 460 (1.17) | 345 (0.88) | 21601.0 |
| ECCp-97 | prime | 97 | 96.5 | yes | 392 | 1370 (3.50) | 487 (1.24) | 365 (0.93) | 8422.5 |
| ECC2-89 | binary | 89 | 83.0 | yes | 360 | 1178 (3.27) | 419 (1.16) | 314 (0.87) | 7345.5 |
| ECCp-89 | prime | 89 | 87.0 | yes | 360 | 1235 (3.43) | 439 (1.22) | 329 (0.92) | 17790.0 |
| ECC2-79 | binary | 79 | 78.8 | yes | 320 | 1118 (3.49) | 398 (1.24) | 298 (0.93) | 7500.0 |
| ECCp-79 | prime | 79 | 78.1 | yes | 320 | 1109 (3.46) | 394 (1.23) | 296 (0.92) | 6517.5 |

Ratio >= 1: no loop endomorphism of that curve can cost less than the doublings a 2-dimensional split saves, under that per-step model.
