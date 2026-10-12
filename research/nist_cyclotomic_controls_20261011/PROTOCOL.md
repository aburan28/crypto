# Protocol: controlled re-test of the comparative claims in RESEARCH_NIST_CYCLOTOMIC_POINT_HEIGHTS.md

**Status: frozen before any data was collected (2026-10-11).** Follows
`.claude/skills/controlled-experiment/SKILL.md`. Any change after the first
run is a new protocol version, reported alongside this one.

## 0. Claim classification of the original note

| original claim (note §) | class | action here |
|:--|:--|:--|
| `f ± 1` factor in `Z[t]` with `Φ_1, Φ_2, Φ_4` (§1) | identity (certified) | none |
| full factorisations of `p ± 1`, `n ± 1`, twist, `t`, `4p − t²` (§1, §2) | exact computation | none |
| embedding degrees `(n−1)/10, (n−1)/3, (n−1)/3` (§2.1) | exact computation | none |
| Cheon cost tables (§2.2) | exact computation given the factorisation | none |
| "P-256's `n − 1` is an outlier / 1-in-30 accident" (§0, §2.2) | **comparative** | **E1** |
| "`N_d`, twist: nothing small enough to matter" (§2.3) | **comparative** (twist part) | **E1** (twist) |
| "canonical lift has no small points" (§3.1) | **"nothing found" without positive control** | **E2** |
| "lifted generator ĥ ≈ log p; nowhere near small" (§3.2) | **comparative** | **E3** |
| "two-point lift independent / orthogonal" (§3.3) | **comparative** | **E4** |
| "the cyclotomic structure never reaches the points" (§0) | out of scope as stated | reword to the exact statement (MOV needs `k`; `k` is exact) |

## E1. Small-factor content of `n − 1`, `n + 1` and the twist order

- **Targets.** P-192, P-224, P-256: `n − 1`, `n + 1`, `T = 2p + 2 − n` (9 values).
- **Null.** For each curve's `p`, group orders of prime-order curves over `F_p` are modelled as **uniform random primes `n'` in the Hasse interval `[p + 1 − 2√p, p + 1 + 2√p]`**. Control values are `n' − 1`, `n' + 1`, `2p + 2 − n'`, computed identically. `N = 400` controls per curve. Seed `20261011 + curve index`.
- **Instrument (identical on every arm).** `curve_traits::arith::factor(x, 2^18)`: trial division to 2^16 then Brent rho, 2^18 iterations per composite.
- **Metric M1.** `log2` of the product, with multiplicity, of the returned prime factors that are `< 2^32`. M1 is defined on what the instrument returns. It is incomplete, but the incompleteness is identical across arms.
- **Metric M2.** The number of distinct returned prime factors `< 2^32`.
- **Positive control.** For each curve, 20 constructed values `x = s · q`. Here `s` is a random 2^20-smooth product of about 120 bits, and `q` is a random prime that brings `x` to the size of `n`. They must score M1 above the 97.5th percentile of the controls. Expected: yes.
- **Decision rule.** The empirical upper-tail fraction `(1 + #{control ≥ target}) / (N + 1)` and the matching lower tail. Two-sided `p = 2·min(tails)`. Bonferroni over 9 targets × 2 metrics = 18 tests, α = 0.05 / 18. A target is **unusual** only if `p < α`. Otherwise it is **within the control distribution (rank r/N)**.

## E2. Small-point census power on the canonical lift

- **Targets.** `y² = x³ − 3x + b` with the NIST `b`, for all three curves.
- **Instrument.** The census in `examples/nist_cyclotomic_point_heights.rs`, unchanged. It scans integer `|x| ≤ 10^6`, and rational `x = u/d²` with `2 ≤ d ≤ 40` and `|u| ≤ 10^5`, behind the residue filter.
- **Positive control.** For each curve's `p`, plant points on 30 curves with `b` of the same size as the NIST `b`. 10 curves get an integer point (`x0` uniform in `|x0| ≤ 10^6`, `y0` uniform below `√p`, `b = y0² − x0³ + 3x0`). 10 get a rational point with `d = 2` and 10 with `d = 3` (`u` uniform in `|u| ≤ 10^5` with `gcd(u, d) = 1`, `v` near `√(p d⁶)`, `b = (v² − u³ + 3u d⁴)/d⁶` when divisible). **Pass:** every planted point is found.
- **Negative control.** 30 curves with `b` uniform in `[0, p)` and no planted point. Expect zero finds.
- **Verdict.** If the positive control passes, "no point in this census" is reported for the NIST curves with the census bounds stated. It never becomes "no small point".

## E3. Canonical height of the lifted generator vs lifted random points

- **Target.** `ĥ(G̃)/ln p` on the one-point lift `y² = x³ − 3x + B` through `(G_x, G_y) ∈ Z²`, using the same height routine as the note (5 doublings).
- **Null.** `N = 200` uniform random points `R ∈ E(F_p)` per curve. Each is lifted identically through `(x_R, y_R)` with representatives in `[0, p)`. Seed fixed.
- **Positive control.** 20 points per curve with `x_R < 2^32` (`x` uniform in `[1, 2^32)` until it lands on the curve), lifted identically. Expected: `h_x` is much smaller. ĥ is reported as found.
- **Metrics.** `ĥ/ln p` and `h_x − ĥ`.
- **Decision rule.** Rank of the target in the null, two-sided, Bonferroni over 3 curves × 2 metrics.

## E4. Height pairing of two-point lifts

- **Target.** `|cos ∠(G̃, R̃)|` and the normalised regulator `det / (ĥ(G̃)ĥ(R̃))` on the two-point lift through `G̃` and `R̃` (`R = [0x1234567890abcdef]G`).
- **Null.** `N = 100` two-point lifts through pairs of uniform random points per curve, constructed identically.
- **Positive control.** On the one-point lift `E_1` through `G̃`, the genuinely dependent rational pairs `(G̃, [2]G̃)` and `(G̃, [3]G̃)`. Expected: `|cos| = 1` and normalised regulator 0 to numerical precision.
- **Decision rule.** Rank in the null, two-sided, Bonferroni over 3 × 2.

## What this protocol cannot show

Nothing here is a security statement. E1 measures one incomplete small-factor
metric. E2 is bounded by its census box. E3 and E4 compare against random
points under one lift construction, and say nothing about lifts chosen
adversarially.
