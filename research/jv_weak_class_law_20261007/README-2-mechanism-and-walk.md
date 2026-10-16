# Protocol 2 results: the mechanism behind the depth-1 exclusion, and the walk policies on ground truth

Run 2026-10-08 against [`PROTOCOL-2-mechanism-and-walk.md`](PROTOCOL-2-mechanism-and-walk.md),
on the exhaustive tables of PR #1556.  **Stage diagnostic, toy sizes.**  `S`,
end-to-end cost and speedup are **unset**.  Nothing here concerns a
prime-field or deployed curve.  Instrument: `census.rs` modes `patterns`,
`walk`, `walkselftest` (standalone Rust, `rustc -O census.rs -o census`).

## 1. Answer

1. **The depth-1 exclusion is a congruence on the public trace.**  With
   `t` the trace over `F_{q³}`, `q = p²`, a class has 2-volcano depth 1
   exactly when `t ≡ ±2 (mod 16)` for `p ≡ ±3 (mod 8)` and exactly when
   `t ≡ ±6 (mod 16)` for `p ≡ ±1 (mod 8)`; the census splits cleanly on
   this at every size (§3).  So A1 of PR #1556 reads: **a weak curve's
   trace never takes those residues**, equivalently `64 | t² − 4q³`,
   equivalently `#E(F_{q³}) ≡ 0 or 4 (mod 16)` on the curve and its
   twist.  The earlier sampled census looked at `t mod 8` and so could not
   see it.
2. **The ascent rule is exact (M1).**  For a curve with exactly one
   halvable 2-torsion point below the crater, the 2-isogeny with that
   kernel ascends: 4,879 of 4,879 such curves at `p = 7`, 73,923 of
   73,923 at `p = 11`.  At the crater it is horizontal, every time.  For
   curves with three halvable points at height 2, the kernel with the
   unique largest 2-adic depth ascends, and when the depths tie the
   codomain with the larger intrinsic 2-adic height does; the two rules
   together resolve every height-2 curve at both sizes.  Above height 3
   the rules leave ties, which the walk never needs.
3. **The exclusion is not about 4-torsion (M2 falsified).**  Weak curves
   with no halvable 2-torsion point exist at height 1 in classes of
   depth ≥ 2 (213 at `p = 7`, 1,088 at `p = 11`), and depth-1 craters hold
   no weak curve in either halvability pattern.  A proof cannot go
   through rational 4-torsion; it has to use the class.
4. **M3 holds**: three halvable points if and only if height ≥ 2, the
   textbook `E[4] ⊂ E(F)` fact, which checks the instrument.
5. **The walk.**  The registered level-walk policy N failed its success
   prediction at `p = 7` in three disclosed versions (0.478, 0.537,
   0.561 against R's 0.825): a single volcano level is a small component,
   and most weak curves sit at height 1 in absolute numbers even though
   the density is highest at height 3.  Policy A, added after that
   failure and labelled so, keeps the refusal and the ascent and then
   walks freely.  Its numbers, and the random walk's, are in §4: the same success rate as R at every size, fewer steps at `p = 7`, about the same by `p = 13`.

## 2. Instrument checks

- Halving formula: for sampled points with a halvable abscissa, some
  half `Q` satisfied `2Q = P` on the curve in 69 of 69 trials
  (`p = 5, 7`).
- 3-isogenies (Vélu on a root of the 3-division polynomial found by
  Cantor–Zassenhaus over `F_{p⁶}`, codomain recovered from the images of
  the 2-torsion points): the codomain was in the same isogeny class, by
  the independently computed trace, and at the same height, in 70 of 70
  trials; zero trace mismatches in every walk below.
- The intrinsic height test (height 1 iff a 2-isogeny descends to the
  floor; height 2 iff a codomain has that property; else ≥ 3) is what
  the walk uses; the census heights are used only to score.

## 3. Mechanism tables

### depth 1 as a congruence on the trace

Classes by `(v₂(f) = 1 ?, |t| mod 16)`, read from `results/census_p{p}.jsonl`:

| p | p mod 8 | depth 1 | depth ≥ 2 |
|--:|--:|:--|:--|
| 7 | 7 | t ≡ 6: 37, t ≡ 10: 37 | t ≡ 2: 37, t ≡ 14: 37 |
| 11 | 3 | t ≡ 2: 152, t ≡ 14: 151 | t ≡ 6: 152, t ≡ 10: 151 |
| 13 | 5 | t ≡ 2: 254, t ≡ 14: 253 | t ≡ 6: 254, t ≡ 10: 254 |
| 17 | 1 | t ≡ 6: 578, t ≡ 10: 578 | t ≡ 2: 579, t ≡ 14: 578 |

Derivation: write `t = 2t′` (full 2-torsion forces `t ≡ 2 mod 4`), so
`t² − 4q³ = 4(t′² − q³)` with `t′` odd and `q³ ≡ 1 (mod 8)`; depth 1 is
`v₂(t′² − q³) = 3`, and `q³ ≡ 1 (mod 16)` iff `p ≡ ±1 (mod 8)`, which fixes
which odd residues of `t′` mod 8 give depth 1.  Every depth-1 class has
`D_K ≡ 0 (mod 8)`, as the tables of PR #1556 showed.

### patterns, p = 7

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

### patterns, p = 11

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

## 4. Walk tables

### walk, p = 7, q = 49, 20 starts per class, cap 3q = 147

| arm | starts | successes (rate) | refused | median steps (successes) | median steps (all) | ascent failures | max ascent steps | missing 3-codomains | trace mismatches |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| R (random) on classes with a weak curve | 1280 | 1056 (0.825) | 0 | 6.0 | 9.0 | 0 | 0 | 0 | 0 |
| N (refuse, ascend, walk level, lower on exhaustion) on classes with a weak curve | 1280 | 718 (0.561) | 0 | 2.0 | 5.0 | 0 | 2 | 0 | 0 |
| A (refuse, ascend, then walk freely) on classes with a weak curve | 1280 | 1056 (0.825) | 0 | 3.0 | 7.0 | 0 | 2 | 0 | 0 |
| R (random) on classes with no weak curve | 1680 | 0 (0.000) | 0 | NaN | 147.0 | 0 | 0 | 0 | 0 |
| N (refuse, ascend, walk level, lower on exhaustion) on classes with no weak curve | 1680 | 0 (0.000) | 1480 | NaN | 0.0 | 0 | 2 | 0 | 0 |
| A (refuse, ascend, then walk freely) on classes with no weak curve | 1680 | 0 (0.000) | 1480 | NaN | 0.0 | 0 | 2 | 0 | 0 |

W-1: starts on weak-holding classes refused: N 0, A 0
W-3 (N (refuse, ascend, walk level, lower on exhaustion) vs R): starts where both succeed 715; median steps R 4.0, this 2.0, ratio 0.500
W-3 (A (refuse, ascend, then walk freely) vs R): starts where both succeed 1052; median steps R 6.0, this 3.0, ratio 0.500

### walk, p = 11, q = 121, 8 starts per class, cap 3q = 363

| arm | starts | successes (rate) | refused | median steps (successes) | median steps (all) | ascent failures | max ascent steps | missing 3-codomains | trace mismatches |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| R (random) on classes with a weak curve | 2168 | 1268 (0.585) | 0 | 14.0 | 75.0 | 0 | 0 | 0 | 0 |
| N (refuse, ascend, walk level, lower on exhaustion) on classes with a weak curve | 2168 | 789 (0.364) | 0 | 3.0 | 154.0 | 0 | 2 | 0 | 0 |
| A (refuse, ascend, then walk freely) on classes with a weak curve | 2168 | 1265 (0.583) | 0 | 11.0 | 65.0 | 0 | 2 | 0 | 0 |
| R (random) on classes with no weak curve | 2680 | 0 (0.000) | 0 | NaN | 206.0 | 0 | 0 | 0 | 0 |
| N (refuse, ascend, walk level, lower on exhaustion) on classes with no weak curve | 2680 | 0 (0.000) | 2424 | NaN | 0.0 | 0 | 2 | 0 | 0 |
| A (refuse, ascend, then walk freely) on classes with no weak curve | 2680 | 0 (0.000) | 2424 | NaN | 0.0 | 0 | 2 | 0 | 0 |

W-1: starts on weak-holding classes refused: N 0, A 0
W-3 (N (refuse, ascend, walk level, lower on exhaustion) vs R): starts where both succeed 782; median steps R 9.0, this 3.0, ratio 0.333
W-3 (A (refuse, ascend, then walk freely) vs R): starts where both succeed 1244; median steps R 13.0, this 10.0, ratio 0.769

### walk, p = 13, q = 169, 8 starts per class, cap 3q = 507

| arm | starts | successes (rate) | refused | median steps (successes) | median steps (all) | ascent failures | max ascent steps | missing 3-codomains | trace mismatches |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| R (random) on classes with a weak curve | 3712 | 1938 (0.522) | 0 | 20.0 | 245.0 | 0 | 0 | 0 | 0 |
| N (refuse, ascend, walk level, lower on exhaustion) on classes with a weak curve | 3712 | 1167 (0.314) | 0 | 5.0 | 154.0 | 0 | 2 | 0 | 0 |
| A (refuse, ascend, then walk freely) on classes with a weak curve | 3712 | 1934 (0.521) | 0 | 18.0 | 248.0 | 0 | 2 | 0 | 0 |
| R (random) on classes with no weak curve | 4408 | 0 (0.000) | 0 | NaN | 206.0 | 0 | 0 | 0 | 0 |
| N (refuse, ascend, walk level, lower on exhaustion) on classes with no weak curve | 4408 | 0 (0.000) | 4056 | NaN | 0.0 | 0 | 2 | 0 | 0 |
| A (refuse, ascend, then walk freely) on classes with no weak curve | 4408 | 0 (0.000) | 4056 | NaN | 0.0 | 0 | 2 | 0 | 0 |

W-1: starts on weak-holding classes refused: N 0, A 0
W-3 (N (refuse, ascend, walk level, lower on exhaustion) vs R): starts where both succeed 1154; median steps R 15.0, this 5.0, ratio 0.333
W-3 (A (refuse, ascend, then walk freely) vs R): starts where both succeed 1906; median steps R 20.0, this 18.0, ratio 0.900

### Reading against the registration

- **W-1 passes**: neither N nor A refused a start on a class that holds a weak curve, at any size; both refused every start on a depth-1 class, which is 88% (`p = 7`) and 90% (`p = 11`) of the starts on classes with no weak curve, at zero cost (92% at `p = 13`), where R walked to its cap or exhaustion (147, 206 and 206 median steps).
- **W-2 fails for N** as registered, at every size (0.561, 0.364 and 0.314 against the registered 0.90), in three disclosed versions at `p = 7` (0.478 with a_est as the height proxy, 0.537 with down-up composites, 0.561 with the exact intrinsic height test and level lowering).  The registered clause on R also fails at `p = 7` (0.825, not below 0.75): this R has 3-isogenies, which the ledger's §17 walk also had, and its exhaustion rule and start distribution differ from the ledger's; the clause was mis-set.  At `p = 11` and `13` R reads 0.585 and 0.522.
- **W-3**: N reaches `0.33`–`0.50×` R's median steps where both succeed, but succeeds on half as many starts.  A, added post hoc, matches R's success rate within a few starts (1,056 vs 1,056; 1,265 vs 1,268; 1,934 vs 1,938) and takes `0.50×` (`p = 7`), `0.77×` (`p = 11`) and `0.90×` (`p = 13`) R's median steps where both succeed; over all starts `0.78×`, `0.87×` and `1.01×`.  The ascent prefix's gain fades with `p`.  The registered `0.40×` is met by neither at `p = 11, 13`.
- **W-4 passes**: the ascent never exceeded 2 steps, below every `v₂(f)`.

What the walk tables say, in order of weight: the refusal is exact and free and removes half of all classes before any step; the ascent prefix is a modest, cost-free gain that fades with `p` (`0.50×`, `0.77×`, `0.90×`); confining the walk to one level loses more in reach than it gains in density, because a level is a small component and most weak curves sit at height 1 in absolute numbers.  The success rate of every arm on weak-holding classes is bounded by the size of the 2,3-components, which only more isogeny degrees change.


## 5. What follows

- The statement to prove is the congruence in §1.1; the facts a proof
  must respect are §1.3 and the crater tables in §3.
- The refusal is exact and free, and is the first thing to wire into
  `jv_isogeny_walk.rs` on branch `research/jv-cover-end-to-end-20261007`.
  Whether the ascent prefix pays there depends on §4 and on the degrees
  5 and 7 that walk has and this instrument lacks.
