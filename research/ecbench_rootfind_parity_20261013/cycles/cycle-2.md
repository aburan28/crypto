# Cycle 2 — `window_div` sweep on `rootfind-window-negfold` (dead end)

Session: `sessions/cycle2` (`ECBS1hd0002abcd05a`), 36/36 verified, audit OK.

## Hypothesis

A longer or shorter cyclic window (`⌈|F|/div⌉`, div ∈ {2,3,4}) on the
negation-folded table could beat cold `S` of `rootfind-negfold` m=2 d=5.

## Mean cold S

| arm | mean S | S / rho |
|---|---:|---:|
| rho-strong | 14.341 | 1.000 |
| ic-f4-m3-d5 | 15.914 | 1.110 |
| **ic-negfold-m2-d5** | **22.481** | 1.568 |
| ic-win-neg-d4 (div=4) | 25.629 | 1.787 |
| ic-win-neg-d3 (div=3) | 26.468 | 1.846 |
| ic-win-neg-d2 (div=2) | 27.459 | 1.915 |

Narrower windows help slightly among m=3 arms, but **all lose to m=2 negfold**.

## Decision

1. **`window_div` sweep is a dead end for cold S** versus m=2 negfold.
2. Keep cycle-1 selection: `pdp3-koblitz:m=2,engine=rootfind-negfold` + d=5.
3. Secondary: amS@64 for negfold-m2-d5 ≈ **16.44** vs F4 ≈ 15.76 — near
   amortized parity, still labelled secondary.

## Next (cycle 3)

Claim panel (8 public × 3 rounds) with negfold m=2 d=5 vs F4 m=3 d=5 and
prior full-rootfind baselines; document cold parity Y/N and amortized secondary.
