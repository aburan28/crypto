### Per family (SGGM bound, group-operation units)

| family | rows | rows with δ_upper < 1 | min measured / min(√(r/A), 1/δ) | min measured / min(√(r/A), 1/δ), δ_upper < 1 only | min measured / √(r/A) | rows below any bound |
|:--|--:|--:|--:|--:|--:|--:|
| char2 | 35 | 23 | 56.777 | 56.777 | 30.129 | 0 |
| koblitz | 6 | 1 | 2.838e6 | 2.838e6 | 92.172 | 0 |
| koblitz_batch | 9 | 1 | 79898 | 1.958e6 | 69.776 | 0 |
| koblitz_compact_orbit (probes, not group ops) | 14 | 0 | 1077 | — | 0.002 | 0 |
| koblitz_round5 | 72 | 13 | 355 | 616 | 32.511 | 0 |
| prime | 72 | 48 | 38.868 | 38.868 | 11.899 | 0 |

### Rows with δ_upper < 1 (the bound is not vacuous)

| configuration | log2 r | A | m | abs F | δ measured | δ upper (count) | 1/δ | √(r/A) | bound | Thm 3.2, ε=½ | measured (GAE) | measured/bound | flags |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|:--|
| icv1-fp10-t5-192cb216/mitm_m2_negfold_walk_balanced | 9.68 | 2 | 2 | 10 | 0.062 | 0.079 | 12.662 | 20.285 | 12.662 | 0.386 | 655 | 51.771 |  |
| icv1-fp10-t5-192cb216/mitm_m2_negfold_walk | 9.68 | 2 | 2 | 32 | 0.417 | 0.680 | 1.470 | 20.285 | 1.470 | 0 | 744 | 506 |  |
| icv1-fp10-t5-192cb216/mitm_m2_negfold | 9.68 | 2 | 2 | 32 | 0.467 | 0.680 | 1.470 | 20.285 | 1.470 | 0 | 1054 | 717 |  |
| icv1-fp10-t5-192cb216/semaev_s3_roots_m2 | 9.68 | 2 | 2 | 32 | 0.467 | 0.680 | 1.470 | 20.285 | 1.470 | 0 | 1095 | 745 |  |
| icv1-fp10-t5-192cb216/mitm_m2 | 9.68 | 2 | 2 | 32 | 0.467 | 0.680 | 1.470 | 20.285 | 1.470 | 0 | 1310 | 891 |  |
| icv1-fp10-t5-192cb216/direct_subtraction_m2 | 9.68 | 2 | 2 | 32 | 0.467 | 0.680 | 1.470 | 20.285 | 1.470 | 0 | 1451 | 987 |  |
| icv1-fp12-t19-d4a315df/mitm_m2_negfold_walk | 11.93 | 2 | 2 | 32 | 0.147 | 0.144 | 6.786 | 44.096 | 6.786 | 0 | 974 | 144 | measured_rate_above_counting_ceiling |
| icv1-fp12-t19-d4a315df/mitm_m2_negfold_walk_balanced | 11.93 | 2 | 2 | 16 | 0.030 | 0.039 | 25.586 | 44.096 | 25.586 | 1.464 | 1033 | 40.392 |  |
| icv1-fp12-t19-d4a315df/mitm_m2_negfold | 11.93 | 2 | 2 | 32 | 0.127 | 0.144 | 6.945 | 44.096 | 6.945 | 0 | 3217 | 463 |  |
| icv1-fp12-t19-d4a315df/mitm_m2 | 11.93 | 2 | 2 | 32 | 0.127 | 0.144 | 6.945 | 44.096 | 6.945 | 0 | 3473 | 500 |  |
| icv1-fp12-t19-d4a315df/semaev_s3_roots_m2 | 11.93 | 2 | 2 | 32 | 0.127 | 0.144 | 6.945 | 44.096 | 6.945 | 0 | 4676 | 673 |  |
| icv1-fp12-t19-d4a315df/direct_subtraction_m2 | 11.93 | 2 | 2 | 32 | 0.127 | 0.144 | 6.945 | 44.096 | 6.945 | 0 | 5701 | 821 |  |
| icv1-fp14-tm43-158b6914/mitm_m2_negfold_walk_balanced | 13.99 | 2 | 2 | 26 | 0.023 | 0.023 | 42.857 | 90.086 | 42.857 | 2.904 | 1666 | 38.868 | measured_rate_above_counting_ceiling |
| icv1-fp14-tm43-158b6914/mitm_m2_negfold_walk | 13.99 | 2 | 2 | 64 | 0.132 | 0.132 | 7.570 | 90.086 | 7.570 | 0 | 1998 | 264 |  |
| icv1-fp14-tm43-158b6914/mitm_m2_negfold | 13.99 | 2 | 2 | 64 | 0.149 | 0.132 | 6.722 | 90.086 | 6.722 | 0 | 7150 | 1064 | measured_rate_above_counting_ceiling |
| icv1-fp14-tm43-158b6914/mitm_m2 | 13.99 | 2 | 2 | 64 | 0.149 | 0.132 | 6.722 | 90.086 | 6.722 | 0 | 8174 | 1216 | measured_rate_above_counting_ceiling |
| icv1-fp14-tm43-158b6914/semaev_s3_roots_m2 | 13.99 | 2 | 2 | 64 | 0.149 | 0.132 | 6.722 | 90.086 | 6.722 | 0 | 12687 | 1887 | measured_rate_above_counting_ceiling |
| icv1-fp14-tm43-158b6914/direct_subtraction_m2 | 13.99 | 2 | 2 | 64 | 0.149 | 0.132 | 6.722 | 90.086 | 6.722 | 0 | 15424 | 2295 | measured_rate_above_counting_ceiling |
| icv1-fp16-t351-d433f24e/mitm_m2_negfold_walk_balanced | 15.99 | 2 | 2 | 40 | 0.013 | 0.013 | 74.102 | 180 | 74.102 | 5.508 | 2979 | 40.198 | measured_rate_above_counting_ceiling |
| icv1-fp16-t351-d433f24e/mitm_m2_negfold_walk | 15.99 | 2 | 2 | 128 | 0.142 | 0.129 | 7.037 | 180 | 7.037 | 0 | 5736 | 815 | measured_rate_above_counting_ceiling |
| icv1-fp16-t351-d433f24e/mitm_m2_negfold | 15.99 | 2 | 2 | 128 | 0.117 | 0.129 | 7.753 | 180 | 7.753 | 0 | 20120 | 2595 |  |
| icv1-fp16-t351-d433f24e/mitm_m2 | 15.99 | 2 | 2 | 128 | 0.117 | 0.129 | 7.753 | 180 | 7.753 | 0 | 24216 | 3123 |  |
| icv1-fp16-t351-d433f24e/direct_subtraction_m2 | 15.99 | 2 | 2 | 128 | 0.117 | 0.129 | 7.753 | 180 | 7.753 | 0 | 58808 | 7585 |  |
| icv1-fp16-t351-d433f24e/semaev_s3_roots_m2 | 15.99 | 2 | 2 | 128 | 0.117 | 0.129 | 7.753 | 180 | 7.753 | 0 | 154746 | 19959 |  |
| icv1-fp18-t175-45770ebc/mitm_m2_negfold_walk_balanced | 18.00 | 2 | 2 | 64 | 0.010 | 0.008 | 104 | 362 | 104 | 8.032 | 4304 | 41.233 | measured_rate_above_counting_ceiling |
| icv1-fp18-t175-45770ebc/mitm_m2_negfold_walk | 18.00 | 2 | 2 | 128 | 0.030 | 0.032 | 31.211 | 362 | 31.211 | 1.934 | 6601 | 212 |  |
| icv1-fp18-t175-45770ebc/mitm_m2_negfold | 18.00 | 2 | 2 | 128 | 0.031 | 0.032 | 31.211 | 362 | 31.211 | 1.934 | 56718 | 1817 |  |
| icv1-fp18-t175-45770ebc/mitm_m2 | 18.00 | 2 | 2 | 128 | 0.031 | 0.032 | 31.211 | 362 | 31.211 | 1.934 | 60814 | 1948 |  |
| icv1-fp18-t175-45770ebc/semaev_s3_roots_m2 | 18.00 | 2 | 2 | 128 | 0.031 | 0.032 | 31.211 | 362 | 31.211 | 1.934 | 154578 | 4953 |  |
| icv1-fp18-t175-45770ebc/direct_subtraction_m2 | 18.00 | 2 | 2 | 128 | 0.031 | 0.032 | 31.211 | 362 | 31.211 | 1.934 | 187797 | 6017 |  |
| icv1-fp20-t1293-3c2999ba/mitm_m2_negfold_walk_balanced | 20.00 | 2 | 2 | 102 | 0.005 | 0.005 | 196 | 724 | 196 | 15.626 | 10972 | 56.117 |  |
| icv1-fp20-t1293-3c2999ba/mitm_m2_negfold_walk | 20.00 | 2 | 2 | 256 | 0.031 | 0.032 | 31.582 | 724 | 31.582 | 1.965 | 20259 | 641 |  |
| icv1-fp20-t1293-3c2999ba/mitm_m2_negfold | 20.00 | 2 | 2 | 256 | 0.038 | 0.032 | 26.589 | 724 | 26.589 | 1.549 | 112963 | 4249 | measured_rate_above_counting_ceiling |
| icv1-fp20-t1293-3c2999ba/mitm_m2 | 20.00 | 2 | 2 | 256 | 0.038 | 0.032 | 26.589 | 724 | 26.589 | 1.549 | 129347 | 4865 | measured_rate_above_counting_ceiling |
| icv1-fp20-t1293-3c2999ba/semaev_s3_roots_m2 | 20.00 | 2 | 2 | 256 | 0.038 | 0.032 | 26.589 | 724 | 26.589 | 1.549 | 434803 | 16353 | measured_rate_above_counting_ceiling |
| icv1-fp20-t1293-3c2999ba/direct_subtraction_m2 | 20.00 | 2 | 2 | 256 | 0.038 | 0.032 | 26.589 | 724 | 26.589 | 1.549 | 536937 | 20194 | measured_rate_above_counting_ceiling |
| icv1-fp22-tm1267-a5acaf6d/mitm_m2_negfold_walk_balanced | 21.65 | 2 | 2 | 148 | 0.003 | 0.003 | 295 | 1283 | 295 | 23.882 | 23169 | 78.650 |  |
| icv1-fp22-tm1267-a5acaf6d/mitm_m2_negfold_walk | 21.65 | 2 | 2 | 512 | 0.037 | 0.040 | 24.967 | 1283 | 24.967 | 1.414 | 71270 | 2855 |  |
| icv1-fp22-tm1267-a5acaf6d/mitm_m2_negfold | 21.65 | 2 | 2 | 512 | 0.040 | 0.040 | 24.967 | 1283 | 24.967 | 1.414 | 253057 | 10136 |  |
| icv1-fp22-tm1267-a5acaf6d/mitm_m2 | 21.65 | 2 | 2 | 512 | 0.040 | 0.040 | 24.967 | 1283 | 24.967 | 1.414 | 318593 | 12760 |  |
| icv1-fp22-tm1267-a5acaf6d/semaev_s3_roots_m2 | 21.65 | 2 | 2 | 512 | 0.040 | 0.040 | 24.967 | 1283 | 24.967 | 1.414 | 1.392e6 | 55772 |  |
| icv1-fp22-tm1267-a5acaf6d/direct_subtraction_m2 | 21.65 | 2 | 2 | 512 | 0.040 | 0.040 | 24.967 | 1283 | 24.967 | 1.414 | 1.761e6 | 70513 |  |
| icv1-fp24-t1577-77336155/mitm_m2_negfold_walk_balanced | 23.38 | 2 | 2 | 222 | 0.002 | 0.002 | 438 | 2338 | 438 | 35.816 | 46179 | 105 |  |
| icv1-fp24-t1577-77336155/mitm_m2_negfold_walk | 23.38 | 2 | 2 | 512 | 0.012 | 0.012 | 82.932 | 2338 | 82.932 | 6.244 | 82157 | 991 |  |
| icv1-fp24-t1577-77336155/mitm_m2_negfold | 23.38 | 2 | 2 | 512 | 0.012 | 0.012 | 82.932 | 2338 | 82.932 | 6.244 | 821584 | 9907 |  |
| icv1-fp24-t1577-77336155/mitm_m2 | 23.38 | 2 | 2 | 512 | 0.012 | 0.012 | 82.932 | 2338 | 82.932 | 6.244 | 887120 | 10697 |  |
| icv1-fp24-t1577-77336155/direct_subtraction_m2 | 23.38 | 2 | 2 | 512 | 0.012 | 0.012 | 82.932 | 2338 | 82.932 | 6.244 | 6.734e6 | 81201 |  |
| icv1-fp24-t1577-77336155/semaev_s3_roots_m2 | 23.38 | 2 | 2 | 512 | 0.012 | 0.012 | 82.932 | 2338 | 82.932 | 6.244 | 1.689e7 | 203701 |  |
| icv1-f2m15-t131-ed7256f0/mitm_m2_negfold_walk | 13.99 | 2 | 2 | 30 | 0.018 | 0.061 | 16.484 | 90.330 | 16.484 | 0.707 | 3361 | 204 |  |
| icv1-f2m15-t131-ed7256f0/mitm_m2_negfold_walk_balanced | 13.99 | 2 | 2 | 14 | 0.003 | 0.015 | 68.567 | 90.330 | 68.567 | 5.047 | 3893 | 56.777 |  |
| icv1-f2m15-t131-ed7256f0/mitm_m3_negfold_walk | 13.99 | 2 | 3 | 30 | 0.100 | 0.669 | 1.496 | 90.330 | 1.496 | 0 | 6785 | 4536 |  |
| icv1-f2m15-t131-ed7256f0/mitm_m3_negfold | 13.99 | 2 | 3 | 30 | 0.144 | 0.669 | 1.496 | 90.330 | 1.496 | 0 | 8587 | 5741 |  |
| icv1-f2m15-t131-ed7256f0/mitm_m3 | 13.99 | 2 | 3 | 30 | 0.144 | 0.669 | 1.496 | 90.330 | 1.496 | 0 | 8812 | 5891 |  |
| icv1-f2m15-t131-ed7256f0/mitm_m2_negfold | 13.99 | 2 | 2 | 30 | 0.013 | 0.061 | 16.484 | 90.330 | 16.484 | 0.707 | 25071 | 1521 |  |
| icv1-f2m15-t131-ed7256f0/semaev_s4_pairs_and_solve_m3 | 13.99 | 2 | 3 | 30 | 0.134 | 0.669 | 1.496 | 90.330 | 1.496 | 0 | 174243 | 116490 |  |
| icv1-f2m18-t393-d2c04152/mitm_m2_negfold_walk_balanced | 15.00 | 2 | 2 | 30 | 0.004 | 0.121 | 8.262 | 128 | 8.262 | 0.022 | 3854 | 466 |  |
| icv1-f2m18-t393-d2c04152/mitm_m2_negfold_walk | 15.00 | 2 | 2 | 64 | 0.017 | 0.524 | 1.908 | 128 | 1.908 | 0 | 6308 | 3307 |  |
| icv1-f2m18-t393-d2c04152/mitm_m2_negfold | 15.00 | 2 | 2 | 64 | 0.016 | 0.524 | 1.908 | 128 | 1.908 | 0 | 40875 | 21428 |  |
| icv1-f2m21-t2591-de54d72e/mitm_m2_negfold_walk_balanced | 20.00 | 2 | 2 | 62 | 0.001 | 0.004 | 260 | 724 | 260 | 20.989 | 24961 | 96.050 |  |
| icv1-f2m21-t2591-de54d72e/mitm_m2_negfold_walk | 20.00 | 2 | 2 | 122 | 0.003 | 0.015 | 68.674 | 724 | 68.674 | 5.056 | 25530 | 372 |  |
| icv1-f2m21-t2591-de54d72e/mitm_m3_negfold_walk | 20.00 | 2 | 3 | 122 | 0.122 | 0.607 | 1.648 | 724 | 1.648 | 0 | 66352 | 40263 |  |
| icv1-f2m21-t2591-de54d72e/mitm_m3_negfold | 20.00 | 2 | 3 | 122 | 0.142 | 0.607 | 1.648 | 724 | 1.648 | 0 | 79156 | 48032 |  |
| icv1-f2m21-t2591-de54d72e/mitm_m3 | 20.00 | 2 | 3 | 122 | 0.142 | 0.607 | 1.648 | 724 | 1.648 | 0 | 82877 | 50290 |  |
| icv1-f2m21-t2591-de54d72e/mitm_m2_negfold | 20.00 | 2 | 2 | 122 | 0.004 | 0.015 | 68.674 | 724 | 68.674 | 5.056 | 570242 | 8304 |  |
| icv1-f2m21-t2591-de54d72e/semaev_s4_pairs_and_solve_m3 | 20.00 | 2 | 3 | 122 | 0.141 | 0.607 | 1.648 | 724 | 1.648 | 0 | 8.720e6 | 5.291e6 |  |
| icv1-f2m24-tm3047-2d0c0e8c/mitm_m2_negfold_walk_balanced | 21.00 | 2 | 2 | 134 | 0.001 | 0.035 | 28.564 | 1024 | 28.564 | 1.714 | 39452 | 1381 |  |
| icv1-f2m24-tm3047-2d0c0e8c/mitm_m2_negfold_walk | 21.00 | 2 | 2 | 274 | 0.005 | 0.145 | 6.909 | 1024 | 6.909 | 0 | 55081 | 7972 |  |
| icv1-f2m24-tm3047-2d0c0e8c/mitm_m2_negfold | 21.00 | 2 | 2 | 274 | 0.004 | 0.145 | 6.909 | 1024 | 6.909 | 0 | 993140 | 143745 |  |
| icv1-f2m27-t12543-569dca8b/mitm_m2_negfold_walk_balanced | 24.41 | 2 | 2 | 254 | 2.578e-4 | 0.009 | 114 | 3344 | 114 | 8.851 | 174655 | 1529 |  |
| icv1-f2m27-t12543-569dca8b/mitm_m2_negfold_walk | 24.41 | 2 | 2 | 526 | 0.001 | 0.037 | 26.795 | 3344 | 26.795 | 1.566 | 197196 | 7359 |  |
| icv1-f2m27-t12543-569dca8b/mitm_m2_negfold | 24.41 | 2 | 2 | 526 | 0.001 | 0.037 | 26.795 | 3344 | 26.795 | 1.566 | 7.778e6 | 290286 |  |
| icv1-f2m41-tm2308219-7f48b14a/koblitz_one_target | 39.00 | 82 | 3 | 6560 | — | 0.343 | 2.918 | 81880 | 2.918 | 0 | 8.281e6 | 2.838e6 |  |
| icv1-f2m41-tm2308219-7f48b14a/koblitz_collection_m3_aimed | 39.00 | 82 | 3 | 6560 | 0.009 | 0.343 | 2.918 | 81880 | 2.918 | 0 | 5.713e6 | 1.958e6 |  |
| icv1-f2m17-tm101-00378d4e/mitm_m2_signed_orbit_columns | 16.00 | 34 | 2 | 239 | — | 0.883 | 1.133 | 43.904 | 1.133 | 0 | 29594 | 26117 |  |
| icv1-f2m17-tm101-00378d4e/mitm_m2_abscissa_columns_control | 16.00 | 34 | 2 | 239 | — | 0.883 | 1.133 | 43.904 | 1.133 | 0 | 38733 | 34183 |  |
| icv1-f2m17-tm101-00378d4e/mitm_m2_signed_orbit_columns_negfold | 16.00 | 34 | 2 | 239 | — | 0.883 | 1.133 | 43.904 | 1.133 | 0 | 15444 | 13630 |  |
| icv1-f2m17-tm101-00378d4e/mitm_m2_signed_orbit_columns_frobfold | 16.00 | 34 | 2 | 239 | — | 0.883 | 1.133 | 43.904 | 1.133 | 0 | 2524 | 2227 |  |
| icv1-f2m17-tm101-00378d4e/mitm_m2_signed_orbit_columns_frobfold_walk | 16.00 | 34 | 2 | 239 | — | 0.883 | 1.133 | 43.904 | 1.133 | 0 | 2323 | 2050 |  |
| icv1-f2m19-tm797-9c54981b/mitm_m2_signed_orbit_columns_frobfold_walk | 18.00 | 38 | 2 | 305 | — | 0.358 | 2.791 | 83.057 | 2.791 | 0 | 3063 | 1098 |  |
| icv1-f2m23-tm5197-1f85e9e1/mitm_m2_signed_orbit_columns_frobfold_walk | 22.00 | 46 | 2 | 875 | — | 0.183 | 5.460 | 302 | 5.460 | 0 | 12685 | 2324 |  |
| icv1-f2m41-tm2308219-7f48b14a/mitm_m3_signed_orbit_columns | 39.00 | 82 | 3 | 5003 | — | 0.152 | 6.577 | 81880 | 6.577 | 0 | 4.362e7 | 6.632e6 |  |
| icv1-f2m41-tm2308219-7f48b14a/mitm_m3_signed_orbit_columns_negfold | 39.00 | 82 | 3 | 5003 | — | 0.152 | 6.577 | 81880 | 6.577 | 0 | 3.736e7 | 5.680e6 |  |
| icv1-f2m41-tm2308219-7f48b14a/mitm_m3_signed_orbit_columns_frobfold | 39.00 | 82 | 3 | 5003 | — | 0.152 | 6.577 | 81880 | 6.577 | 0 | 3.348e7 | 5.091e6 |  |
| icv1-f2m41-tm2308219-7f48b14a/mitm_m3_signed_orbit_columns_frobfold_walk | 39.00 | 82 | 3 | 5003 | — | 0.152 | 6.577 | 81880 | 6.577 | 0 | 3.044e7 | 4.629e6 |  |
| icv1-f2m41-tm2308219-7f48b14a/mitm_m2_signed_orbit_columns_frobfold_walk | 39.00 | 82 | 2 | 5003 | — | 9.111e-5 | 10975 | 81880 | 10975 | 914 | 6.764e6 | 616 |  |
| icv1-f2m41-tm2308219-7f48b14a/mitm_m2_signed_orbit_columns_frobfold_walk_balanced | 39.00 | 82 | 2 | 20501 | — | 0.002 | 654 | 81880 | 654 | 53.827 | 3.793e6 | 5800 |  |

### Compact-orbit ladder (4-sum relations; probes are S3 solve + root-table lookup, not group operations)

| rung | log2 r | A | B = abs F | K | δ upper | measured | probes | probes / √(r/A) | S = K²n | log2(S·T²) − log2(r/A) | S·T/(r/A) | frozen / rank mean | P(T ≤ frozen) if exponential |
|:--|--:|--:|--:|--:|--:|:--|--:|--:|--:|--:|--:|--:|--:|
| compact-n37 | 27.78 | 74 | 518 | 7 | 1.000 | target_probes_online | 1077 | 0.610 | 1813 | 9.4 | 0.627 | 0.573 | 0.436 |
| compact-n37 | 27.78 | 74 | 518 | 7 | 1.000 | rank_plus_target_probes | 14237 | 8.065 | 1813 | 16.8 | 8.283 | — | — |
| compact-n41 | 39.00 | 82 | 6970 | 85 | 1.000 | target_probes_online | 21248 | 0.260 | 296225 | 14.3 | 0.939 | 0.860 | 0.577 |
| compact-n41 | 39.00 | 82 | 6970 | 85 | 1.000 | rank_plus_target_probes | 2.121e6 | 25.910 | 296225 | 27.6 | 93.736 | — | — |
| compact-n53 | 44.26 | 106 | 23320 | 220 | 1.000 | target_probes_online | 26485 | 0.059 | 2.565e6 | 13.1 | 0.342 | 0.316 | 0.271 |
| compact-n53 | 44.26 | 106 | 23320 | 220 | 1.000 | rank_plus_target_probes | 1.844e7 | 41.385 | 2.565e6 | 32.0 | 238 | — | — |
| compact-n61 | 47.21 | 122 | 73200 | 600 | 1.000 | target_probes_online | 227437 | 0.197 | 2.196e7 | 19.7 | 3.741 | 3.554 | 0.971 |
| compact-n61 | 47.21 | 122 | 73200 | 600 | 1.000 | rank_plus_target_probes | 3.862e7 | 33.426 | 2.196e7 | 34.5 | 635 | — | — |
| compact-n71 | 52.29 | 142 | 85200 | 600 | 1.000 | target_probes_online | 14554 | 0.002 | 2.556e7 | 7.1 | 0.010 | 0.006 | 0.006 |
| compact-n71 | 52.29 | 142 | 85200 | 600 | 1.000 | rank_plus_target_probes | 1.347e9 | 216 | 2.556e7 | 40.1 | 887 | — | — |
| compact-n73 | 56.26 | 146 | 87600 | 600 | 1.000 | target_probes_online | 457561 | 0.019 | 2.628e7 | 13.2 | 0.020 | 0.007 | 0.007 |
| compact-n73 | 56.26 | 146 | 87600 | 600 | 1.000 | rank_plus_target_probes | 3.678e10 | 1515 | 2.628e7 | 45.8 | 1641 | — | — |
| compact-n83 | 52.93 | 166 | 99600 | 600 | 1.000 | target_probes_online | 8.845e6 | 1.231 | 2.988e7 | 25.4 | 5.120 | 3.613 | 0.973 |
| compact-n83 | 52.93 | 166 | 99600 | 600 | 1.000 | rank_plus_target_probes | 1.478e9 | 206 | 2.988e7 | 40.2 | 855 | — | — |

At the rank-stage mean probes per relation, S·T/(r/A) is: compact-n37 1.094, compact-n41 1.092, compact-n53 1.081, compact-n61 1.053, compact-n71 1.478, compact-n73 2.734, compact-n83 1.417.

### Multiplier closure (non-adaptive table arm)

| arm | r | A | pairs | canonical entries | mean first-hit probes T | pairs·T/r | canonical·T/(r/A) | logs | target |
|:--|--:|--:|--:|--:|--:|--:|--:|:--|:--|
| multiplier-closure/c23/frob-K12-c1 | 2095853 | 46 | 151800 | 3163 | 15.094 | 1.093 | 1.048 | true | true |
| multiplier-closure/c23/tau-K6-c2 | 2095853 | 46 | 151800 | 3117 | 18.766 | 1.359 | 1.284 | true | true |
| multiplier-closure/c23/int-K6-c2 | 2095853 | 46 | 151800 | 3095 | 13.734 | 0.995 | 0.933 | true | true |
| multiplier-closure/c23/tau-K4-c3 | 2095853 | 46 | 151800 | 3034 | 14.562 | 1.055 | 0.970 | true | true |
| multiplier-closure/c23/int-K4-c3 | 2095853 | 46 | 151800 | 3015 | 16.219 | 1.175 | 1.073 | true | true |
| multiplier-closure/c23/tau-K2-c6 | 2095853 | 46 | 151800 | 2666 | 19.609 | 1.420 | 1.147 | true | true |
| multiplier-closure/c23/int-K2-c6 | 2095853 | 46 | 151800 | 2672 | 12.609 | 0.913 | 0.739 | true | true |
| multiplier-closure/c23/tau-K1-c12 | 2095853 | 46 | 151800 | 1834 | 26.062 | 1.888 | 1.049 | true | true |
| multiplier-closure/c23/int-K1-c12 | 2095853 | 46 | 151800 | 1961 | 26.172 | 1.896 | 1.126 | true | true |
| multiplier-closure/c37/frob-K36-c1 | 230603167 | 74 | 3545784 | 47454 | 63.859 | 0.982 | 0.972 | false | false |
| multiplier-closure/c37/tau-K18-c2 | 230603167 | 74 | 3545784 | 47160 | 78.141 | 1.202 | 1.183 | true | true |
| multiplier-closure/c37/int-K18-c2 | 230603167 | 74 | 3545784 | 47270 | 68.969 | 1.060 | 1.046 | true | true |
| multiplier-closure/c37/tau-K9-c4 | 230603167 | 74 | 3545784 | 46580 | 74.438 | 1.145 | 1.113 | true | true |
| multiplier-closure/c37/int-K9-c4 | 230603167 | 74 | 3545784 | 46625 | 74.125 | 1.140 | 1.109 | true | true |
| multiplier-closure/c37/tau-K3-c12 | 230603167 | 74 | 3545784 | 43201 | 80.406 | 1.236 | 1.115 | true | true |
| multiplier-closure/c37/int-K3-c12 | 230603167 | 74 | 3545784 | 43457 | 66.984 | 1.030 | 0.934 | true | true |
| multiplier-closure/c37/tau-K1-c36 | 230603167 | 74 | 3545784 | 30471 | 109 | 1.678 | 1.067 | true | true |
| multiplier-closure/c37/int-K1-c36 | 230603167 | 74 | 3545784 | 33907 | 95.375 | 1.466 | 1.038 | true | true |
| multiplier-closure/c37b/frob-K36-c1 | 230603167 | 74 | 3545784 | 47473 | 68.039 | 1.046 | 1.037 | true | true |
| multiplier-closure/c37b/tau-K18-c2 | 230603167 | 74 | 3545784 | 47255 | 75.312 | 1.158 | 1.142 | true | true |
| multiplier-closure/c37b/int-K18-c2 | 230603167 | 74 | 3545784 | 47263 | 70.719 | 1.087 | 1.073 | true | true |
| multiplier-closure/c37b/tau-K9-c4 | 230603167 | 74 | 3545784 | 46595 | 59.953 | 0.922 | 0.896 | true | true |
| multiplier-closure/c37b/int-K9-c4 | 230603167 | 74 | 3545784 | 46676 | 68.070 | 1.047 | 1.020 | true | true |
| multiplier-closure/c37b/tau-K3-c12 | 230603167 | 74 | 3545784 | 42792 | 80.297 | 1.235 | 1.103 | true | true |
| multiplier-closure/c37b/int-K3-c12 | 230603167 | 74 | 3545784 | 43594 | 75.070 | 1.154 | 1.050 | true | true |
| multiplier-closure/c37b/tau-K1-c36 | 230603167 | 74 | 3545784 | 30471 | 146 | 2.251 | 1.432 | true | true |
| multiplier-closure/c37b/int-K1-c36 | 230603167 | 74 | 3545784 | 33907 | 101 | 1.559 | 1.103 | true | true |

### Record checks

- `docs/ic/boundary_targets.json /regimes/koblitz/records/vs_rho/history/2/metrics/subgroup_log2_order`: recorded log2 r = 47.5407, log2 of the exact r = 47.2109
- `docs/ic/boundary_targets.json /regimes/koblitz/records/vs_rho/history/2/metrics/ic_precomputation_ms_excluded_from_online` (n = 61): recorded excluded precompute 13098 ms; the rung's own components are 653 / 12441 / 11053 / 0.781 (base load / index build / rank stage / linear algebra, ms); includes rank stage: false; includes index build: true
- `docs/ic/boundary_targets.json /regimes/koblitz/records/vs_rho/history/1/metrics/ic_precomputation_ms_excluded_from_online` (n = 71): recorded excluded precompute 39113 ms; the rung's own components are 290 / 38807 / 1.155e6 / — (base load / index build / rank stage / linear algebra, ms); includes rank stage: false; includes index build: true

### Multi-target slope d ln(total gae) / d ln k

| source | curve | log2 r | A | k | slope | 95% bootstrap | Kuhn–Struik law slope (same k) | all verified |
|:--|:--|--:|--:|:--|--:|:--|--:|:--|
| committed records | icv1-f2m41-tm2308219-7f48b14a | 39.00 | 82 | [1, 4, 16, 32] | 0.547 | [0.487, 0.622] | 0.533 | true (ledger §19.3) |
| committed records | icv1-f2m53-tm56619371-dac20a85 | 44.26 | 106 | [1, 4, 16, 32] | 0.585 | [0.504, 0.671] | 0.533 | true (ledger §19.3) |
| committed records | icv1-f2m61-t158598901-ab42b6c5 | 47.21 | 122 | [1, 32] | 0.521 | [0.416, 0.681] | 0.534 | true (ledger §19.3) |
| fresh measurement | prime-floor q=68723001073 N=68723001023 (ell=120283, t=51) | 36.00 | 2 | [1, 2, 4, 8, 16] | 0.526 | [0.453, 0.608] | 0.539 | true |
| fresh measurement | icv1-f2m41-tm2308219-7f48b14a | 39.00 | 82 | [1, 2, 4, 8, 16] | 0.561 | [0.500, 0.637] | 0.539 | true |

| source | curve | k | batches | mean total gae | total / √(k r/A) | per target over k=1 | batch law |
|:--|:--|--:|--:|--:|--:|--:|--:|
| records | icv1-f2m41-tm2308219-7f48b14a | 1 | 16 | 105937 | 1.294 | 1.000 | 1.000 |
| records | icv1-f2m41-tm2308219-7f48b14a | 4 | 16 | 219965 | 1.343 | 0.519 | 0.547 |
| records | icv1-f2m41-tm2308219-7f48b14a | 16 | 16 | 469266 | 1.433 | 0.277 | 0.280 |
| records | icv1-f2m41-tm2308219-7f48b14a | 32 | 16 | 710725 | 1.534 | 0.210 | 0.199 |
| records | icv1-f2m53-tm56619371-dac20a85 | 1 | 16 | 482740 | 1.083 | 1.000 | 1.000 |
| records | icv1-f2m53-tm56619371-dac20a85 | 4 | 16 | 1.149e6 | 1.290 | 0.595 | 0.547 |
| records | icv1-f2m53-tm56619371-dac20a85 | 16 | 16 | 2.496e6 | 1.400 | 0.323 | 0.280 |
| records | icv1-f2m53-tm56619371-dac20a85 | 32 | 16 | 3.693e6 | 1.465 | 0.239 | 0.199 |
| records | icv1-f2m61-t158598901-ab42b6c5 | 1 | 8 | 1.469e6 | 1.271 | 1.000 | 1.000 |
| records | icv1-f2m61-t158598901-ab42b6c5 | 32 | 8 | 8.941e6 | 1.368 | 0.190 | 0.199 |
| fresh | prime-floor q=68723001073 N=68723001023 (ell=120283, t=51) | 1 | 24 | 255072 | 1.376 | 1.000 | 1.000 |
| fresh | prime-floor q=68723001073 N=68723001023 (ell=120283, t=51) | 2 | 24 | 383172 | 1.462 | 0.751 | 0.750 |
| fresh | prime-floor q=68723001073 N=68723001023 (ell=120283, t=51) | 4 | 24 | 554337 | 1.495 | 0.543 | 0.547 |
| fresh | prime-floor q=68723001073 N=68723001023 (ell=120283, t=51) | 8 | 24 | 727279 | 1.387 | 0.356 | 0.393 |
| fresh | prime-floor q=68723001073 N=68723001023 (ell=120283, t=51) | 16 | 24 | 1.146e6 | 1.546 | 0.281 | 0.280 |
| fresh | icv1-f2m41-tm2308219-7f48b14a | 1 | 24 | 96511 | 1.179 | 1.000 | 1.000 |
| fresh | icv1-f2m41-tm2308219-7f48b14a | 2 | 24 | 173602 | 1.499 | 0.899 | 0.750 |
| fresh | icv1-f2m41-tm2308219-7f48b14a | 4 | 24 | 246785 | 1.507 | 0.639 | 0.547 |
| fresh | icv1-f2m41-tm2308219-7f48b14a | 8 | 24 | 339394 | 1.465 | 0.440 | 0.393 |
| fresh | icv1-f2m41-tm2308219-7f48b14a | 16 | 24 | 481947 | 1.472 | 0.312 | 0.280 |
