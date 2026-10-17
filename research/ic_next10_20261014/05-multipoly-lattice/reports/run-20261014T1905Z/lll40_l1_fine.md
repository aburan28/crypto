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
| V3 S3 pair s=1 | 1 | max [4, 4] | 45 | 0.09474 | 0.120 | 27 | 24/24 | 7.62 | 3501 |
| V3 S3 pair s=1 | 1 | max [4, 4] | 45 | 0.09474 | 0.130 | 36 | 17/24 | 4.33 | 2847 |
| V3 S3 pair s=1 | 1 | max [4, 4] | 45 | 0.09474 | 0.140 | 48 | 9/24 | 2.46 | 2172 |
| V3 S3 pair s=1 | 1 | max [4, 4] | 45 | 0.09474 | 0.150 | 63 | 4/24 | 1.00 | 1912 |
| V3 S3 pair s=1 | 1 | max [4, 4] | 45 | 0.09474 | 0.160 | 84 | 2/24 | 0.96 | 1747 |

total wall 396.8 s
wrote lll40_l1_fine.json
