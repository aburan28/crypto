# Cycle 1 — negation-folded PairTable (`rootfind-negfold`)

Session: `sessions/cycle1` (`ECBS1hcd05f1dc18a3`), 48/48 verified, audit OK,
replay cert `ad54d3921a14e6ed7c253af3c9b24a634f83cb51cf030a80b7133c09d2330536`.

## Hypothesis

`PairTable::build_negation_folded` (~|F|²/4 set-up adds) cuts cold `S`
versus full `PairSumTable` (|F|(|F|+1)/2) without changing completeness.

## Mean cold S (n=13, 2 targets × 2 rounds)

| arm | mean S | S / rho | setup GAE |
|---|---:|---:|---:|
| rho-strong | 13.961 | 1.000 | — |
| ic-f4-m3-d5 | 16.120 | 1.155 | 0 (lower-bound S) |
| **ic-negfold-m2-d5** | **23.515** | 1.684 | 305 |
| ic-window-negfold-m3-d5 | 26.826 | 1.922 | 305 |
| ic-negfold-m2-d4 | 29.251 | 2.095 | 109 |
| ic-rootfind-m2-d5 (full) | 29.236 | 2.094 | 561 |
| ic-rootfind-m2-d4 (full) | 31.061 | 2.225 | 190 |
| ic-rootfind-window-m3-d5 | 32.546 | 2.331 | 561 |

Paired ops: `ic-negfold-m2-d5 / ic-f4-m3-d5` = 1.459 [1.404, 1.493] (**bounded**).

Relations GAE for negfold vs full is **identical** (789.8 at d=5); the entire
cold win is set-up (561 → 305).

## Decision

1. **Keep `rootfind-negfold` m=2 d=5** as cold candidate (23.5 vs prior best 28.5).
2. Window+negfold helps vs full window (26.8 vs 32.5) but loses to m=2 negfold.
3. Parity still **not** achieved (need ≤16.32). Remaining gap ≈ residual set-up
   (~305 GAE) plus m=2 online yield vs F4 m=3.
4. Secondary amS@64 for negfold-m2-d5 ≈ 18.1 (still above F4 ≈ 15.4).

## Next (cycle 2)

Parameterize `window_div` on `rootfind-window-negfold` and smoke-sweep
�

## Next (cycle 2)

Parameterize `window_div` on `rootfind-window-negfold` and smoke-sweep
⌈|F|/2⌉, ⌈|F|/3⌉, ⌈|F|/4⌉ to hunt a cold-S win over m=2 negfold.
