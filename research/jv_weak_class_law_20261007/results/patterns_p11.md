## patterns, p = 11

| height | pattern | nodes | weak | weak fraction × q |
|--:|--:|--:|--:|--:|
| 1 | 0 | 73657 | 1088 | 1.79 |
| 1 | 1 | 147788 | 2542 | 2.08 |
| 2 | 3 | 55278 | 1837 | 4.02 |
| 3 | 3 | 16183 | 1842 | 13.77 |
| 4 | 3 | 2053 | 0 | 0.00 |
| 5 | 3 | 269 | 0 | 0.00 |
| 6 | 3 | 28 | 0 | 0.00 |
| 7 | 3 | 3 | 0 | 0.00 |

M1 — unique halvable point's isogeny, (below crater, direction: +1 up, −1 down, 0 level, −99 floor) → count:
- below_crater=false dir=0 : 73865
- below_crater=true dir=1 : 73923

M2 — weak nodes at height 1 with pattern 0 in depth ≥ 2 classes: 1088

Crater test — (D_K mod 8, at crater) → nodes, weak, weak fraction × q:
- D_K mod 8 = 0, at crater = false : 41608 nodes, 48 weak, 0.14
- D_K mod 8 = 0, at crater = true : 126509 nodes, 1116 weak, 1.07
- D_K mod 8 = 1, at crater = false : 28356 nodes, 1537 weak, 6.56
- D_K mod 8 = 1, at crater = true : 7876 nodes, 960 weak, 14.75
- D_K mod 8 = 4, at crater = false : 41112 nodes, 1836 weak, 5.40
- D_K mod 8 = 4, at crater = true : 15732 nodes, 252 weak, 1.94
- D_K mod 8 = 5, at crater = false : 31371 nodes, 1560 weak, 6.02
- D_K mod 8 = 5, at crater = true : 2695 nodes, 0 weak, 0.00

Height 1 by class depth — (depth = 1, pattern) → nodes, weak:
- depth1 = false, pattern 0 : 36633 nodes, 1088 weak
- depth1 = false, pattern 1 : 73923 nodes, 2542 weak
- depth1 = true, pattern 0 : 37024 nodes, 0 weak
- depth1 = true, pattern 1 : 73865 nodes, 0 weak

Pattern-3 ascent rule — (below crater, outcome) → count:
- below_crater = false : max-depth point level : 20746
- below_crater = false : tie -> a_est fallback level : 13688
- below_crater = false : tie -> unresolved : 7489
- below_crater = true : max-depth point ascends : 13888
- below_crater = true : tie -> a_est fallback ascends : 13878
- below_crater = true : tie -> unresolved : 4125

Pattern-3 ascent rule by height — (height, outcome) → count:
- height 2 : max-depth point ascends : 13888
- height 2 : max-depth point level : 13824
- height 2 : tie -> a_est fallback ascends : 13878
- height 2 : tie -> a_est fallback level : 13688
- height 3 : max-depth point level : 6922
- height 3 : tie -> unresolved : 9261
- height 4 : tie -> unresolved : 2053
- height 5 : tie -> unresolved : 269
- height 6 : tie -> unresolved : 28
- height 7 : tie -> unresolved : 3
M3 — pattern-3 nodes: 73814, minimum height 2

Intrinsic estimate a_est = 1 + min_i depth2(T_i) (cap 3) against height (both capped at 4): equal 285967, different 9292
| height | a_est | nodes |
|--:|--:|--:|
| 1 | 1 | 221445 |
| 2 | 2 | 55278 |
| 3 | 2 | 6939 |
| 3 | 3 | 9244 |
| 4 | 2 | 2353 |