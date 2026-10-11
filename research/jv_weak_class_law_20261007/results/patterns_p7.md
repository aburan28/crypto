## patterns, p = 7

| height | pattern | nodes | weak | weak fraction × q |
|--:|--:|--:|--:|--:|
| 1 | 0 | 4958 | 213 | 2.11 |
| 1 | 1 | 9748 | 375 | 1.89 |
| 2 | 3 | 3655 | 294 | 3.94 |
| 3 | 3 | 1064 | 243 | 11.19 |
| 4 | 3 | 154 | 60 | 19.09 |
| 5 | 3 | 23 | 12 | 25.57 |
| 6 | 3 | 6 | 3 | 24.50 |

M1 — unique halvable point's isogeny, (below crater, direction: +1 up, −1 down, 0 level, −99 floor) → count:
- below_crater=false dir=0 : 4869
- below_crater=true dir=1 : 4879

M2 — weak nodes at height 1 with pattern 0 in depth ≥ 2 classes: 213

Crater test — (D_K mod 8, at crater) → nodes, weak, weak fraction × q:
- D_K mod 8 = 0, at crater = false : 2568 nodes, 396 weak, 7.56
- D_K mod 8 = 0, at crater = true : 8416 nodes, 24 weak, 0.14
- D_K mod 8 = 1, at crater = false : 2304 nodes, 156 weak, 3.32
- D_K mod 8 = 1, at crater = true : 552 nodes, 126 weak, 11.18
- D_K mod 8 = 4, at crater = false : 2568 nodes, 114 weak, 2.18
- D_K mod 8 = 4, at crater = true : 1020 nodes, 186 weak, 8.94
- D_K mod 8 = 5, at crater = false : 1995 nodes, 126 weak, 3.09
- D_K mod 8 = 5, at crater = true : 185 nodes, 72 weak, 19.07

Height 1 by class depth — (depth = 1, pattern) → nodes, weak:
- depth1 = false, pattern 0 : 2431 nodes, 213 weak
- depth1 = false, pattern 1 : 4879 nodes, 375 weak
- depth1 = true, pattern 0 : 2527 nodes, 0 weak
- depth1 = true, pattern 1 : 4869 nodes, 0 weak

Pattern-3 ascent rule — (below crater, outcome) → count:
- below_crater = false : max-depth point level : 1083
- below_crater = false : tie -> a_est fallback level : 1008
- below_crater = false : tie -> unresolved : 686
- below_crater = true : max-depth point ascends : 1035
- below_crater = true : tie -> a_est fallback ascends : 1051
- below_crater = true : tie -> unresolved : 39

Pattern-3 ascent rule by height — (height, outcome) → count:
- height 2 : max-depth point ascends : 894
- height 2 : max-depth point level : 912
- height 2 : tie -> a_est fallback ascends : 937
- height 2 : tie -> a_est fallback level : 912
- height 3 : max-depth point ascends : 141
- height 3 : max-depth point level : 96
- height 3 : tie -> a_est fallback ascends : 114
- height 3 : tie -> a_est fallback level : 96
- height 3 : tie -> unresolved : 617
- height 4 : max-depth point level : 75
- height 4 : tie -> unresolved : 79
- height 5 : tie -> unresolved : 23
- height 6 : tie -> unresolved : 6
M3 — pattern-3 nodes: 4902, minimum height 2

Intrinsic estimate a_est = 1 + min_i depth2(T_i) (cap 3) against height (both capped at 4): equal 18905, different 703
| height | a_est | nodes |
|--:|--:|--:|
| 1 | 1 | 14706 |
| 2 | 2 | 3655 |
| 3 | 2 | 617 |
| 3 | 3 | 447 |
| 4 | 3 | 86 |
| 4 | 4 | 97 |