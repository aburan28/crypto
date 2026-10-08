
### Kernel->isogeny, 32-bit p, task `codomain_from_point` (median)

| l | velu | kohel(poly+formulas) | sqrt_velu |
|---|---|---|---|
| 3 | 340 ns | 814 ns | 441 ns |
| 5 | 430 ns | 1.1 us | 580 ns |
| 7 | 785 ns | 2.1 us | 1.8 us |
| 11 | 1.5 us | 2.7 us | 2.7 us |
| 13 | 1.6 us | 4.0 us | 3.4 us |
| 31 | 4.4 us | 12.3 us | 4.1 us |
| 61 | 8.2 us | 26.3 us | 6.8 us |
| 101 | 14.0 us | 60.9 us | 10.7 us |
| 211 | 29.6 us | 226.8 us | 20.7 us |
| 401 | 56.3 us | 776.8 us | 34.6 us |
| 1009 | 142.6 us | 4.51 ms | 80.3 us |

### Kernel->isogeny, 32-bit p, task `codomain_given_h` (median)

| l | kohel(h given) |
|---|---|
| 3 | 491 ns |
| 5 | 658 ns |
| 7 | 799 ns |
| 11 | 1.2 us |
| 13 | 2.2 us |
| 31 | 7.1 us |
| 61 | 16.3 us |
| 101 | 41.6 us |
| 211 | 172.6 us |
| 401 | 644.7 us |
| 1009 | 3.87 ms |

### Kernel->isogeny, 32-bit p, task `eval_point` (median)

| l | velu | kohel | sqrt_velu |
|---|---|---|---|
| 3 | 251 ns | 589 ns | 455 ns |
| 5 | 593 ns | 604 ns | 538 ns |
| 7 | 718 ns | 605 ns | 1.5 us |
| 11 | 1.4 us | 908 ns | 1.7 us |
| 13 | 1.7 us | 1.0 us | 2.7 us |
| 31 | 3.6 us | 1.7 us | 2.6 us |
| 61 | 7.2 us | 3.2 us | 4.6 us |
| 101 | 11.8 us | 5.0 us | 7.1 us |
| 211 | 25.0 us | 10.2 us | 14.6 us |
| 401 | 48.6 us | 19.2 us | 26.1 us |
| 1009 | 123.7 us | 47.7 us | 60.3 us |

### Kernel->isogeny, 61-bit p, task `codomain_from_point` (median)

| l | velu | kohel(poly+formulas) | sqrt_velu |
|---|---|---|---|
| 3 | 417 ns | 900 ns | 522 ns |
| 5 | 1.0 us | 1.7 us | 1.2 us |
| 7 | 1.6 us | 2.4 us | 2.9 us |
| 11 | 2.2 us | 3.5 us | 3.9 us |
| 13 | 2.7 us | 4.3 us | 3.6 us |
| 31 | 7.0 us | 12.5 us | 5.9 us |
| 61 | 14.1 us | 32.0 us | 9.2 us |
| 101 | 23.0 us | 70.5 us | 13.3 us |
| 211 | 49.3 us | 255.7 us | 26.2 us |
| 401 | 93.1 us | 801.0 us | 42.8 us |

### Kernel->isogeny, 61-bit p, task `codomain_given_h` (median)

| l | kohel(h given) |
|---|---|
| 3 | 487 ns |
| 5 | 646 ns |
| 7 | 813 ns |
| 11 | 1.5 us |
| 13 | 1.5 us |
| 31 | 5.1 us |
| 61 | 16.4 us |
| 101 | 41.4 us |
| 211 | 178.2 us |
| 401 | 623.5 us |

### Kernel->isogeny, 61-bit p, task `eval_point` (median)

| l | velu | kohel | sqrt_velu |
|---|---|---|---|
| 3 | 527 ns | 867 ns | 746 ns |
| 5 | 953 ns | 1.1 us | 814 ns |
| 7 | 1.3 us | 1.2 us | 1.9 us |
| 11 | 2.2 us | 1.4 us | 2.0 us |
| 13 | 2.9 us | 1.4 us | 2.2 us |
| 31 | 6.3 us | 2.3 us | 3.1 us |
| 61 | 13.0 us | 3.7 us | 5.1 us |
| 101 | 21.7 us | 5.8 us | 7.8 us |
| 211 | 45.6 us | 11.6 us | 15.1 us |
| 401 | 90.4 us | 20.5 us | 26.4 us |

### Kernel finding, 30-bit p (median; Phi setup is a single run)

| l | #kernels | divpoly factor | Phi roots+Elkies+BMSS | Phi roots only | Phi_l setup (once) |
|---|---|---|---|---|---|
| 3 | 1 | 11.2 us | 26.3 us | 7.9 us | 34.8 us |
| 5 | 2 | 223.3 us | 163.3 us | 31.5 us | 138.9 us |
| 7 | 1 | 2.62 ms | 211.9 us | 20.6 us | 347.8 us |
| 11 | 2 | 18.51 ms | 1.63 ms | 52.3 us | 2.43 ms |
| 13 | 2 | 27.16 ms | 2.84 ms | 53.1 us | 5.01 ms |
| 17 | 2 | 113.24 ms | 7.58 ms | 121.0 us | 19.50 ms |
| 19 | 2 | 825.48 ms | 11.35 ms | 96.1 us | 44.62 ms |
| 23 | 1 | 1.23 s | 11.92 ms | 136.2 us | 95.24 ms |

### Kernel finding, 61-bit p (median; Phi setup is a single run)

| l | #kernels | divpoly factor | Phi roots+Elkies+BMSS | Phi roots only | Phi_l setup (once) |
|---|---|---|---|---|---|
| 3 | 1 | 20.5 us | 36.0 us | 15.4 us | 43.2 us |
| 5 | 2 | 457.2 us | 177.5 us | 45.9 us | 114.8 us |
| 7 | 2 | 2.20 ms | 427.4 us | 55.3 us | 328.0 us |
| 11 | 1 | 49.23 ms | 873.9 us | 81.6 us | 2.27 ms |
| 13 | 2 | 38.18 ms | 2.99 ms | 120.4 us | 5.01 ms |
| 17 | 2 | 234.85 ms | 7.99 ms | 182.6 us | 19.67 ms |
| 19 | 2 | 815.62 ms | 11.75 ms | 370.3 us | 34.20 ms |
| 23 | 2 | 2.10 s | 23.01 ms | 334.0 us | 150.49 ms |

### Path finding (median over verified instances)

| algorithm | p bits | verified | median time | path len | work counters (median) |
|---|---|---|---|---|---|
| couveignes_hhs_mitm | 20 | 5/5 | 216.92 ms | 4 | nodes=28 |
| couveignes_hhs_mitm | 24 | 5/5 | 566.93 ms | 6 | nodes=143 |
| couveignes_hhs_mitm | 28 | 4/4 | 482.14 ms | 6 | nodes=205 |
| delfs_galbraith_supersingular | 22 | 5/5 | 97.22 ms | 1231 | walk_steps=1075 bfs_nodes=182 |
| delfs_galbraith_supersingular | 26 | 5/5 | 386.11 ms | 4174 | walk_steps=3767 bfs_nodes=947 |
| delfs_galbraith_supersingular | 30 | 5/5 | 2.80 s | 23503 | walk_steps=23470 bfs_nodes=5829 |
| delfs_galbraith_supersingular | 34 | 5/5 | 7.06 s | 30950 | walk_steps=30870 bfs_nodes=53227 |
| delfs_galbraith_supersingular | 38 | 3/5 | 37.81 s | 282752 | walk_steps=253781 bfs_nodes=107099 |
| explicit_chain_from_path | 20 | 5/5 | 299.5 us | - | steps=6 |
| explicit_chain_from_path | 24 | 5/5 | 883.4 us | - | steps=12 |
| explicit_chain_from_path | 28 | 5/5 | 1.26 ms | - | steps=14 |
| explicit_chain_from_path | 32 | 5/5 | 986.6 us | - | steps=6 |
| galbraith_bidirectional_bfs | 20 | 5/5 | 792.7 us | 6 | nodes_expanded=19 |
| galbraith_bidirectional_bfs | 24 | 5/5 | 4.95 ms | 12 | nodes_expanded=105 |
| galbraith_bidirectional_bfs | 28 | 5/5 | 3.59 ms | 14 | nodes_expanded=61 |
| galbraith_bidirectional_bfs | 32 | 5/5 | 862.7 us | 6 | nodes_expanded=14 |
| ghs_random_walk_collision | 20 | 5/5 | 2.62 ms | 35 | steps=63 |
| ghs_random_walk_collision | 24 | 5/5 | 10.05 ms | 102 | steps=200 |
| ghs_random_walk_collision | 28 | 5/5 | 9.95 ms | 112 | steps=184 |
| ghs_random_walk_collision | 32 | 5/5 | 27.03 ms | 238 | steps=475 |
| ghs_with_kohel_volcano | 17 | 4/5 | 1.47 ms | 10 | steps=12 |
| ghs_with_kohel_volcano | 19 | 4/5 | 425.9 us | 4 | steps=3 |
| ghs_with_kohel_volcano | 20 | 5/5 | 2.56 ms | 35 | steps=63 |
| ghs_with_kohel_volcano | 21 | 4/5 | 4.22 ms | 18 | steps=20 |
| ghs_with_kohel_volcano | 24 | 5/5 | 9.91 ms | 102 | steps=200 |
| ghs_with_kohel_volcano | 28 | 5/5 | 10.22 ms | 112 | steps=184 |
| ghs_with_kohel_volcano | 32 | 5/5 | 27.82 ms | 238 | steps=475 |
| kohel_volcano_crater_walk | 17 | 5/5 | 393.4 us | 2 |  |
| kohel_volcano_crater_walk | 19 | 5/5 | 225.8 us | 2 |  |
| kohel_volcano_crater_walk | 21 | 5/5 | 544.7 us | 2 |  |

Unverified / failed instances (kept in the raw data):

- ghs_with_kohel_volcano p_bits=17 instance=1 ns=36657562 note: ascend l=3, walk with l=5,7,11,13; fails when the split primes do not generate the class group orbit
- ghs_with_kohel_volcano p_bits=19 instance=0 ns=168506957 note: ascend l=3, walk with l=5,7,11,13; fails when the split primes do not generate the class group orbit
- ghs_with_kohel_volcano p_bits=21 instance=3 ns=132004351 note: ascend l=3, walk with l=5,7,11,13; fails when the split primes do not generate the class group orbit
- delfs_galbraith_supersingular p_bits=38 instance=3 ns=45346931958 note: j-path only; Fp2 arithmetic
- delfs_galbraith_supersingular p_bits=38 instance=4 ns=95265923987 note: j-path only; Fp2 arithmetic
