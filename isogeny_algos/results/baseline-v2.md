
### V2 kernel -> isogeny, odd degree, 32-bit p (median): general Velu from P, x-only Velu from x(P), Kohel with h given

| l | velu_general_from_point | velu_xonly_from_x | kohel_given_h |
|---|---|---|---|
| 3 | 287 ns | 189 ns | 1.4 us |
| 5 | 566 ns | 380 ns | 2.0 us |
| 7 | 831 ns | 672 ns | 2.6 us |
| 11 | 1.5 us | 1.2 us | 3.0 us |
| 13 | 1.7 us | 1.4 us | 3.3 us |
| 31 | 4.4 us | 4.0 us | 7.8 us |
| 61 | 8.8 us | 8.5 us | 20.2 us |
| 101 | 14.1 us | 14.0 us | 47.1 us |
| 211 | 29.4 us | 29.6 us | 185.2 us |
| 401 | 56.4 us | 56.1 us | 650.1 us |
| 1009 | 144.1 us | 145.2 us | 4.06 ms |

### V2 kernel -> isogeny, even / composite cyclic degree n, 32-bit p (median)

| n | velu_general_from_point | kohel_given_h |
|---|---|---|
| 2 | 90 ns | 856 ns |
| 4 | 617 ns | 1.7 us |
| 6 | 1.1 us | 4.3 us |
| 8 | 1.3 us | 2.5 us |
| 10 | 1.4 us | 3.0 us |
| 12 | 1.6 us | 3.3 us |
| 16 | 2.6 us | 4.1 us |
| 20 | 2.7 us | 4.9 us |

### V2 Montgomery x-only Velu vs Weierstrass Velu for the same CSIDH kernel (median)

| l | montgomery_xonly_velu_codomain | weierstrass_velu_same_kernel |
|---|---|---|
| 3 | 249 ns | 398 ns |
| 5 | 517 ns | 692 ns |
| 7 | 464 ns | 1.0 us |
| 11 | 504 ns | 1.6 us |
| 13 | 643 ns | 2.0 us |
| 17 | 737 ns | 2.8 us |
| 19 | 686 ns | 2.9 us |
| 23 | 846 ns | 4.0 us |
| 29 | 1.0 us | 4.8 us |
| 43 | 1.3 us | 7.1 us |

### V2 l^e-kernel chains over F_p^2 (median); identical codomain for all strategies

| l | e | p bits | strategy | time | l-mults | point evals | isogeny builds |
|---|---|---|---|---|---|---|---|
| 2 | 8 | 16 | balanced | 9.4 us | 12 | 12 | 8 |
| 2 | 8 | 16 | naive | 14.4 us | 28 | 7 | 8 |
| 2 | 8 | 16 | optimal_calibrated | 10.7 us | 10 | 15 | 8 |
| 2 | 12 | 28 | balanced | 20.5 us | 24 | 20 | 12 |
| 2 | 12 | 28 | naive | 42.3 us | 66 | 11 | 12 |
| 2 | 12 | 28 | optimal_calibrated | 23.3 us | 19 | 25 | 12 |
| 2 | 16 | 35 | balanced | 36.7 us | 32 | 32 | 16 |
| 2 | 16 | 35 | naive | 89.0 us | 120 | 15 | 16 |
| 2 | 16 | 35 | optimal_calibrated | 41.5 us | 29 | 36 | 16 |
| 2 | 20 | 43 | balanced | 58.4 us | 48 | 40 | 20 |
| 2 | 20 | 43 | naive | 162.6 us | 190 | 19 | 20 |
| 2 | 20 | 43 | optimal_calibrated | 66.6 us | 37 | 52 | 20 |
| 2 | 24 | 49 | balanced | 80.9 us | 60 | 52 | 24 |
| 2 | 24 | 49 | naive | 263.1 us | 276 | 23 | 24 |
| 2 | 24 | 49 | optimal_calibrated | 96.1 us | 43 | 75 | 24 |
| 3 | 5 | 16 | balanced | 7.9 us | 7 | 5 | 5 |
| 3 | 5 | 16 | naive | 9.7 us | 10 | 4 | 5 |
| 3 | 5 | 16 | optimal_calibrated | 7.7 us | 5 | 7 | 5 |
| 3 | 8 | 28 | balanced | 19.4 us | 12 | 12 | 8 |
| 3 | 8 | 28 | naive | 31.2 us | 28 | 7 | 8 |
| 3 | 8 | 28 | optimal_calibrated | 19.4 us | 10 | 15 | 8 |
| 3 | 10 | 35 | balanced | 31.5 us | 19 | 15 | 10 |
| 3 | 10 | 35 | naive | 56.2 us | 45 | 9 | 10 |
| 3 | 10 | 35 | optimal_calibrated | 30.7 us | 13 | 23 | 10 |
| 3 | 12 | 43 | balanced | 46.2 us | 24 | 20 | 12 |
| 3 | 12 | 43 | naive | 92.8 us | 66 | 11 | 12 |
| 3 | 12 | 43 | optimal_calibrated | 44.6 us | 17 | 29 | 12 |
| 3 | 14 | 49 | balanced | 61.1 us | 29 | 25 | 14 |
| 3 | 14 | 49 | naive | 137.7 us | 91 | 13 | 14 |
| 3 | 14 | 49 | optimal_calibrated | 60.8 us | 21 | 36 | 14 |

### V2 (E, E~) -> isogeny, 40-bit p (median; sigma supplied where required)

| l | linear_algebra | stark1972 | atkin1992 | atkin_modcomp | elkies1992 | elkies1998 | fast_elkies | fast_elkies_prime | reference_kohel_given_g | v1_pade_weierstrass_series |
|---|---|---|---|---|---|---|---|---|---|---|
| 5 | 15.9 us | 46.9 us | 4.2 us | 5.9 us | 2.8 us | 2.4 us | 8.1 us | 69.7 us | 2.1 us | 63.0 us |
| 7 | 51.5 us | 86.4 us | 7.5 us | 8.7 us | 3.9 us | 3.3 us | 12.9 us | 140.0 us | 2.8 us | 223.7 us |
| 11 | 255.1 us | 224.0 us | 10.9 us | 20.1 us | 5.6 us | 4.5 us | 21.5 us | 382.5 us | 4.2 us | 799.4 us |
| 13 | 522.3 us | 305.5 us | 13.7 us | 19.5 us | 6.3 us | 4.8 us | 24.8 us | 546.7 us | 3.5 us | 1.45 ms |
| 17 | 1.33 ms | 586.7 us | 23.4 us | 32.5 us | 8.8 us | 6.4 us | 39.3 us | 1.03 ms | 4.4 us | 3.72 ms |
| 19 | 2.09 ms | 747.5 us | 29.1 us | 40.3 us | 10.0 us | 7.1 us | 50.0 us | 1.37 ms | 4.8 us | 5.60 ms |
| 23 | 4.53 ms | 1.23 ms | 43.7 us | 59.9 us | 13.0 us | 8.8 us | 72.4 us | 2.49 ms | 5.8 us | 11.90 ms |
| 29 | 11.40 ms | 2.30 ms | 75.3 us | 99.6 us | 18.3 us | 11.7 us | 113.0 us | 4.25 ms | 7.5 us | 27.34 ms |
| 31 | 14.89 ms | 2.64 ms | 93.4 us | 122.7 us | 19.9 us | 12.6 us | 139.1 us | 4.97 ms | 7.9 us | 36.57 ms |
| 37 | 29.89 ms | 4.33 ms | 153.6 us | 186.7 us | 26.6 us | 16.3 us | 206.1 us | 8.17 ms | 10.1 us | 66.89 ms |
| 41 | 44.50 ms | 5.91 ms | 183.5 us | 231.5 us | 30.5 us | 18.5 us | 255.4 us | 10.91 ms | 11.4 us | 126.16 ms |
| 53 | 119.34 ms | 11.39 ms | 370.1 us | 473.2 us | 45.7 us | 26.9 us | 488.3 us | 22.27 ms | 16.7 us | 275.01 ms |
| 61 | 211.67 ms | 26.02 ms | 532.8 us | 643.6 us | 57.2 us | 32.8 us | 677.2 us | 33.27 ms | 20.5 us | 496.22 ms |
| 101 | 1.65 s | 73.08 ms | 2.18 ms | 2.51 ms | 147.1 us | 74.3 us | 2.59 ms | 145.73 ms | 47.5 us | 3.62 s |

### V2 sigma from Phi and dual isogeny (median)

| l | sigma_from_phi | dual_isogeny |
|---|---|---|
| 5 | 4.8 us | 64.5 us |
| 7 | 8.9 us | 234.3 us |
| 11 | 13.8 us | 815.4 us |
| 13 | 17.8 us | 1.47 ms |
| 17 | 31.1 us | 3.81 ms |
| 19 | 33.9 us | 5.67 ms |
| 23 | 48.0 us | 11.38 ms |

### V2 end-to-end from (E, l) with Phi_l given (median)

| l | e2e_elkies1998 | e2e_fast_elkies_prime |
|---|---|---|
| 5 | 28.2 us | 98.4 us |
| 7 | 78.1 us | 375.5 us |
| 11 | 121.0 us | 838.9 us |
| 13 | 144.0 us | 1.17 ms |
| 17 | 279.1 us | 2.32 ms |
| 19 | 275.0 us | 2.97 ms |
| 23 | 367.0 us | 4.74 ms |

### V2 CSIDH-style action

| what | parameters / counters | time | verified |
|---|---|---|---|
| class_group_action | {"p_bits": "20", "n_primes": "6", "exp_bound": "3", "isogeny_steps": "14", "ns_per_step": 2007.2857142857142} | 28.1 us | yes |
| mitm_group_action_inversion | {"p_bits": "20", "n_primes": "6", "box_m": "2", "ns": 3294758, "nodes": 1458} | 3.29 ms | yes |
| class_group_action | {"p_bits": "30", "n_primes": "8", "exp_bound": "3", "isogeny_steps": "6", "ns_per_step": 2973.1666666666665} | 17.8 us | yes |
| mitm_group_action_inversion | {"p_bits": "30", "n_primes": "8", "box_m": "2", "ns": 44173903, "nodes": 13122} | 44.17 ms | yes |
| class_group_action | {"p_bits": "40", "n_primes": "10", "exp_bound": "3", "isogeny_steps": "13", "ns_per_step": 4037.769230769231} | 52.5 us | yes |
| mitm_group_action_inversion | {"p_bits": "40", "n_primes": "10", "box_m": "1", "ns": 8897349, "nodes": 2048} | 8.90 ms | yes |
| class_group_action | {"p_bits": "50", "n_primes": "12", "exp_bound": "3", "isogeny_steps": "17", "ns_per_step": 5017.941176470588} | 85.3 us | yes |
| ideal_order_by_cycle | {"p_bits": "20", "ell": "3", "class_number": "1905", "ns": 5431041, "order": 1905} | 5.43 ms | yes |
| ideal_order_by_cycle | {"p_bits": "20", "ell": "5", "class_number": "1905", "ns": 3761358, "order": 1905} | 3.76 ms | yes |
| ideal_order_by_cycle | {"p_bits": "20", "ell": "7", "class_number": "1905", "ns": 3675585, "order": 1905} | 3.68 ms | yes |

### V2 Kohel End(E) conductor (median)

| p bits | instance | height of 3-volcano | level of E | time | verified |
|---|---|---|---|---|---|
| 17 | 0 | 2 | 2 | 4.8 us | yes |
| 17 | 1 | 1 | 1 | 4.9 us | yes |
| 17 | 2 | 2 | 2 | 4.8 us | yes |
| 19 | 0 | 1 | 1 | 5.4 us | yes |
| 19 | 1 | 2 | 2 | 5.3 us | yes |
| 19 | 2 | 1 | 1 | 5.3 us | yes |
| 21 | 0 | 1 | 0 | 65.2 us | yes |
| 21 | 1 | 1 | 1 | 6.5 us | yes |
| 21 | 2 | 2 | 1 | 84.6 us | yes |

### V2 walks on identical instances (median over verified instances)

| algorithm | p bits | verified | median time | median steps / nodes |
|---|---|---|---|---|
| galbraith_bfs | 24 | 5/5 | 10.95 ms | 83 |
| galbraith_stolbunov_w_16_8_4_2_1 | 24 | 5/5 | 2.61 ms | 115 |
| galbraith_stolbunov_w_1_1_1_1_1 | 24 | 5/5 | 4.35 ms | 93 |
| ghs_uniform_over_all_neighbours | 24 | 5/5 | 16.37 ms | 92 |
| galbraith_bfs | 28 | 5/5 | 11.83 ms | 71 |
| galbraith_stolbunov_w_16_8_4_2_1 | 28 | 5/5 | 4.56 ms | 211 |
| galbraith_stolbunov_w_1_1_1_1_1 | 28 | 5/5 | 2.81 ms | 66 |
| ghs_uniform_over_all_neighbours | 28 | 5/5 | 28.04 ms | 188 |
| galbraith_bfs | 32 | 5/5 | 17.50 ms | 93 |
| galbraith_stolbunov_w_16_8_4_2_1 | 32 | 5/5 | 16.95 ms | 732 |
| galbraith_stolbunov_w_1_1_1_1_1 | 32 | 5/5 | 2.35 ms | 48 |
| ghs_uniform_over_all_neighbours | 32 | 5/5 | 16.97 ms | 95 |
