# Frozen toy relation-geometry tables

Structural diagnostics only. All S, rho ratios, DLP-floor ratios and full-DLP costs are null. Rank is the dimension of span(V*V), not a measured solving degree. Pairs use distinct affine points and exclude infinity sums from covered targets. All families in each cell use the same curve and coordinate dimension.

G = generic_linear; U = subfield_linear; scaled U = scaled_subfield_linear.

## Initial frozen run

| n/k | Curve case | l | Family | Base points | Product rank | Distinct targets / curve order | Pair attempts | Coverage / counting ceiling |
|---|---|---:|---|---:|---:|---:|---:|---:|
| 6/2 | subfield_control | 2 | G | 3 | 3 | 2 / 56 | 3 | 0.666667 |
| 6/2 | subfield_control | 2 | U | 7 | 2 | 7 / 56 | 21 | 0.333333 |
| 6/2 | subfield_control | 2 | scaled U | 3 | 2 | 2 / 56 | 3 | 0.666667 |
| 6/2 | subfield_control | 4 | G | 15 | 6 | 55 / 56 | 105 | 1.000000 |
| 6/2 | subfield_control | 4 | U | 15 | 6 | 45 / 56 | 105 | 0.818182 |
| 6/2 | subfield_control | 4 | scaled U | 15 | 6 | 53 / 56 | 105 | 0.963636 |
| 6/2 | full_degree_main | 2 | G | 7 | 3 | 18 / 64 | 21 | 0.857143 |
| 6/2 | full_degree_main | 2 | U | 7 | 2 | 18 / 64 | 21 | 0.857143 |
| 6/2 | full_degree_main | 2 | scaled U | 3 | 2 | 2 / 64 | 3 | 0.666667 |
| 6/2 | full_degree_main | 4 | G | 13 | 6 | 38 / 64 | 78 | 0.603175 |
| 6/2 | full_degree_main | 4 | U | 15 | 6 | 52 / 64 | 105 | 0.825397 |
| 6/2 | full_degree_main | 4 | scaled U | 15 | 6 | 51 / 64 | 105 | 0.809524 |
| 6/2 | full_degree_holdout | 2 | G | 5 | 3 | 8 / 64 | 10 | 0.800000 |
| 6/2 | full_degree_holdout | 2 | U | 7 | 2 | 18 / 64 | 21 | 0.857143 |
| 6/2 | full_degree_holdout | 2 | scaled U | 7 | 2 | 18 / 64 | 21 | 0.857143 |
| 6/2 | full_degree_holdout | 4 | G | 13 | 6 | 52 / 64 | 78 | 0.825397 |
| 6/2 | full_degree_holdout | 4 | U | 19 | 6 | 59 / 64 | 171 | 0.936508 |
| 6/2 | full_degree_holdout | 4 | scaled U | 15 | 6 | 54 / 64 | 105 | 0.857143 |
| 9/3 | subfield_control | 3 | G | 9 | 6 | 32 / 508 | 36 | 0.888889 |
| 9/3 | subfield_control | 3 | U | 3 | 3 | 2 / 508 | 3 | 0.666667 |
| 9/3 | subfield_control | 3 | scaled U | 11 | 3 | 50 / 508 | 55 | 0.909091 |
| 9/3 | subfield_control | 6 | G | 55 | 9 | 483 / 508 | 1485 | 0.952663 |
| 9/3 | subfield_control | 6 | U | 59 | 9 | 498 / 508 | 1711 | 0.982249 |
| 9/3 | subfield_control | 6 | scaled U | 67 | 9 | 502 / 508 | 2211 | 0.990138 |
| 9/3 | full_degree_main | 3 | G | 7 | 6 | 18 / 512 | 21 | 0.857143 |
| 9/3 | full_degree_main | 3 | U | 7 | 3 | 18 / 512 | 21 | 0.857143 |
| 9/3 | full_degree_main | 3 | scaled U | 11 | 3 | 50 / 512 | 55 | 0.909091 |
| 9/3 | full_degree_main | 6 | G | 71 | 9 | 511 / 512 | 2485 | 1.000000 |
| 9/3 | full_degree_main | 6 | U | 59 | 9 | 502 / 512 | 1711 | 0.982387 |
| 9/3 | full_degree_main | 6 | scaled U | 71 | 9 | 511 / 512 | 2485 | 1.000000 |
| 9/3 | full_degree_holdout | 3 | G | 9 | 6 | 32 / 480 | 36 | 0.888889 |
| 9/3 | full_degree_holdout | 3 | U | 7 | 3 | 18 / 480 | 21 | 0.857143 |
| 9/3 | full_degree_holdout | 3 | scaled U | 11 | 3 | 50 / 480 | 55 | 0.909091 |
| 9/3 | full_degree_holdout | 6 | G | 45 | 9 | 429 / 480 | 990 | 0.895616 |
| 9/3 | full_degree_holdout | 6 | U | 59 | 9 | 471 / 480 | 1711 | 0.983299 |
| 9/3 | full_degree_holdout | 6 | scaled U | 63 | 9 | 477 / 480 | 1953 | 0.995825 |
| 10/2 | subfield_control | 2 | G | 5 | 3 | 8 / 968 | 10 | 0.800000 |
| 10/2 | subfield_control | 2 | U | 7 | 2 | 7 / 968 | 21 | 0.333333 |
| 10/2 | subfield_control | 2 | scaled U | 3 | 2 | 2 / 968 | 3 | 0.666667 |
| 10/2 | subfield_control | 4 | G | 13 | 10 | 72 / 968 | 78 | 0.923077 |
| 10/2 | subfield_control | 4 | U | 19 | 6 | 149 / 968 | 171 | 0.871345 |
| 10/2 | subfield_control | 4 | scaled U | 11 | 6 | 50 / 968 | 55 | 0.909091 |
| 10/2 | full_degree_main | 2 | G | 5 | 3 | 8 / 1024 | 10 | 0.800000 |
| 10/2 | full_degree_main | 2 | U | 7 | 2 | 18 / 1024 | 21 | 0.857143 |
| 10/2 | full_degree_main | 2 | scaled U | 3 | 2 | 2 / 1024 | 3 | 0.666667 |
| 10/2 | full_degree_main | 4 | G | 21 | 8 | 184 / 1024 | 210 | 0.876190 |
| 10/2 | full_degree_main | 4 | U | 19 | 6 | 150 / 1024 | 171 | 0.877193 |
| 10/2 | full_degree_main | 4 | scaled U | 11 | 6 | 50 / 1024 | 55 | 0.909091 |
| 10/2 | full_degree_holdout | 2 | G | 3 | 3 | 2 / 1024 | 3 | 0.666667 |
| 10/2 | full_degree_holdout | 2 | U | 7 | 2 | 18 / 1024 | 21 | 0.857143 |
| 10/2 | full_degree_holdout | 2 | scaled U | 3 | 2 | 2 / 1024 | 3 | 0.666667 |
| 10/2 | full_degree_holdout | 4 | G | 5 | 9 | 8 / 1024 | 10 | 0.800000 |
| 10/2 | full_degree_holdout | 4 | U | 15 | 6 | 98 / 1024 | 105 | 0.933333 |
| 10/2 | full_degree_holdout | 4 | scaled U | 15 | 6 | 94 / 1024 | 105 | 0.895238 |
| 12/4 | subfield_control | 4 | G | 17 | 9 | 128 / 4144 | 136 | 0.941176 |
| 12/4 | subfield_control | 4 | U | 15 | 4 | 15 / 4144 | 105 | 0.142857 |
| 12/4 | subfield_control | 4 | scaled U | 11 | 4 | 50 / 4144 | 55 | 0.909091 |
| 12/4 | subfield_control | 8 | G | 293 | 12 | 4143 / 4144 | 42778 | 1.000000 |
| 12/4 | subfield_control | 8 | U | 243 | 12 | 4139 / 4144 | 29403 | 0.999035 |
| 12/4 | subfield_control | 8 | scaled U | 247 | 12 | 4141 / 4144 | 30381 | 0.999517 |
| 12/4 | full_degree_main | 4 | G | 19 | 10 | 162 / 4128 | 171 | 0.947368 |
| 12/4 | full_degree_main | 4 | U | 23 | 4 | 242 / 4128 | 253 | 0.956522 |
| 12/4 | full_degree_main | 4 | scaled U | 15 | 4 | 96 / 4128 | 105 | 0.914286 |
| 12/4 | full_degree_main | 8 | G | 267 | 12 | 4127 / 4128 | 35511 | 1.000000 |
| 12/4 | full_degree_main | 8 | U | 243 | 12 | 4126 / 4128 | 29403 | 0.999758 |
| 12/4 | full_degree_main | 8 | scaled U | 247 | 12 | 4124 / 4128 | 30381 | 0.999273 |
| 12/4 | full_degree_holdout | 4 | G | 15 | 10 | 98 / 4096 | 105 | 0.933333 |
| 12/4 | full_degree_holdout | 4 | U | 15 | 4 | 98 / 4096 | 105 | 0.933333 |
| 12/4 | full_degree_holdout | 4 | scaled U | 15 | 4 | 98 / 4096 | 105 | 0.933333 |
| 12/4 | full_degree_holdout | 8 | G | 259 | 12 | 4095 / 4096 | 33411 | 1.000000 |
| 12/4 | full_degree_holdout | 8 | U | 247 | 12 | 4093 / 4096 | 30381 | 0.999512 |
| 12/4 | full_degree_holdout | 8 | scaled U | 275 | 12 | 4095 / 4096 | 37675 | 1.000000 |

## Independent odd-prime holdout

| n/k | Curve case | l | Family | Base points | Product rank | Distinct targets / curve order | Pair attempts | Coverage / counting ceiling |
|---|---|---:|---|---:|---:|---:|---:|---:|
| 6/2 | odd_prime_holdout | 2 | G | 3 | 3 | 2 / 52 | 3 | 0.666667 |
| 6/2 | odd_prime_holdout | 2 | U | 3 | 2 | 2 / 52 | 3 | 0.666667 |
| 6/2 | odd_prime_holdout | 2 | scaled U | 3 | 2 | 2 / 52 | 3 | 0.666667 |
| 6/2 | odd_prime_holdout | 4 | G | 11 | 6 | 37 / 52 | 55 | 0.725490 |
| 6/2 | odd_prime_holdout | 4 | U | 11 | 6 | 37 / 52 | 55 | 0.725490 |
| 6/2 | odd_prime_holdout | 4 | scaled U | 15 | 6 | 46 / 52 | 105 | 0.901961 |
| 9/3 | odd_prime_holdout | 3 | G | 7 | 5 | 18 / 548 | 21 | 0.857143 |
| 9/3 | odd_prime_holdout | 3 | U | 11 | 3 | 50 / 548 | 55 | 0.909091 |
| 9/3 | odd_prime_holdout | 3 | scaled U | 3 | 3 | 2 / 548 | 3 | 0.666667 |
| 9/3 | odd_prime_holdout | 6 | G | 71 | 9 | 541 / 548 | 2485 | 0.989031 |
| 9/3 | odd_prime_holdout | 6 | U | 71 | 9 | 543 / 548 | 2485 | 0.992687 |
| 9/3 | odd_prime_holdout | 6 | scaled U | 67 | 9 | 543 / 548 | 2211 | 0.992687 |
| 10/2 | odd_prime_holdout | 2 | G | 1 | 3 | 0 / 1048 | 0 | undefined (zero pairs) |
| 10/2 | odd_prime_holdout | 2 | U | 7 | 2 | 18 / 1048 | 21 | 0.857143 |
| 10/2 | odd_prime_holdout | 2 | scaled U | 7 | 2 | 18 / 1048 | 21 | 0.857143 |
| 10/2 | odd_prime_holdout | 4 | G | 15 | 10 | 86 / 1048 | 105 | 0.819048 |
| 10/2 | odd_prime_holdout | 4 | U | 23 | 6 | 187 / 1048 | 253 | 0.739130 |
| 10/2 | odd_prime_holdout | 4 | scaled U | 23 | 6 | 204 / 1048 | 253 | 0.806324 |
| 12/4 | odd_prime_holdout | 4 | G | 13 | 10 | 72 / 4184 | 78 | 0.923077 |
| 12/4 | odd_prime_holdout | 4 | U | 23 | 4 | 242 / 4184 | 253 | 0.956522 |
| 12/4 | odd_prime_holdout | 4 | scaled U | 19 | 4 | 162 / 4184 | 171 | 0.947368 |
| 12/4 | odd_prime_holdout | 8 | G | 241 | 12 | 4181 / 4184 | 28920 | 0.999522 |
| 12/4 | odd_prime_holdout | 8 | U | 291 | 12 | 4183 / 4184 | 42195 | 1.000000 |
| 12/4 | odd_prime_holdout | 8 | scaled U | 263 | 12 | 4182 / 4184 | 34453 | 0.999761 |
