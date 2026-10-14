# IC next-10 item 05: multi-polynomial / shared-target Coppersmith reach

## LLL validation: p = 1099511627563 (40-bit), a = -3, b = 31, δ_LLL = 0.99

Success = at least n reduced rows vanish over Z at the planted root (n = 2 single, 3 pair).

| system | ℓ | M (closure-optimal) | dim | LP δ_ℓ | δ | B = ⌊p^δ⌋ | success | mean vanishing rows | mean LLL ms |
|:--|--:|:--|--:|--:|--:|--:|:--|--:|--:|
| C1 S3 | 1 | max [4] | 25 | 0.09000 | 0.110 | 21 | 24/24 | 5.96 | 78 |
| C1 S3 | 1 | max [4] | 25 | 0.09000 | 0.120 | 27 | 19/24 | 3.50 | 54 |
| C1 S3 | 1 | max [4] | 25 | 0.09000 | 0.130 | 36 | 11/24 | 1.88 | 47 |
| C1 S3 | 1 | max [4] | 25 | 0.09000 | 0.140 | 48 | 7/24 | 0.92 | 52 |
| C1 S3 | 1 | max [4] | 25 | 0.09000 | 0.150 | 63 | 6/24 | 0.79 | 55 |
| C1 S3 | 1 | max [4] | 25 | 0.09000 | 0.160 | 84 | 1/24 | 0.29 | 54 |
| V3 S3 pair s=1 | 1 | max [4, 4] | 45 | 0.09474 | 0.110 | 21 | 24/24 | 12.29 | 4005 |
