# Prime-order curve search in O(√p) (2026-10-17)

`find_prime_order_curve` (`src/cryptanalysis/ic_boundary.rs`) built every
`prime_search` workload by counting points with a Legendre symbol at all `p`
abscissae for each candidate curve, about `p` Jacobi symbols per candidate
and roughly `ln p` candidates per accepted curve. `ecbench plan`, `run`,
`profile-input`, `verify` and every CI replay of a `prime_search` spec paid
it: 11 s per 24-bit curve and 218 s per 26-bit curve on this host.

`PrimeCurve::prime_order` answers the only question the search asks, "is
`#E(F_p)` prime, and what is it", by stepping one point across the Hasse
interval `[p+1−t, p+1+t]`, `t = ⌊2√p⌋`: about `4√p` affine additions. If
`N = #E` is prime, every point `P ≠ ∞` has order `N`, `2N` lies above the
interval (`p > 37`), so `N` is the only `m` in the interval with
`[m]P = ∞`. One prime hit is therefore `N`; no point, several hits, or one
composite hit all mean `N` is not prime. It draws no randomness, so the
search consumes the same random numbers and returns the same curves.

## Results

`ecbench plan` on a one-curve `prime_search` spec (`bits*.json`, seed
59297), median of 3 runs (main at 26 bits: one run).

| bits | main (s) | this change (s) | speedup |
|--:|--:|--:|--:|
| 16 | 0.212 | 0.040 | 5× |
| 18 | 0.484 | 0.041 | 12× |
| 20 | 0.997 | 0.036 | 28× |
| 22 | 6.901 | 0.064 | 109× |
| 24 | 10.956 | 0.046 | 239× |
| 26 | 217.607 | 0.137 | 1,592× |

Below 22 bits most of the remaining time is process start-up and the rest
of planning. Raw timings: `timing-main.txt`, `timing-new.txt`. Figure and
report: `plan-time.svg`, `report.html` (`python3 make_figure.py`), copied to
[`benchmarks/20261017-prime-order-search/`](../../benchmarks/20261017-prime-order-search/).

## Identity of the output

- All 21 committed sessions on main whose spec uses `prime_search` re-plan
  to exactly the workload IDs in their `plan.json`: 288 workloads
  (`replan-check.txt`).
- Unit test `prime_order_matches_the_naive_count`: `prime_order` equals
  `point_count` + primality on 400 random curves (`p` < 20,000, prime,
  composite and even orders), and `find_prime_order_curve` equals the
  naive search at 8–16 bits for seeds 1, 59297 and 20261003.

## Scope

Engineering change to workload construction only: no solver, cost unit,
record or measured result changes, and no comparison or IC claim.
`point_count` stays for callers that need a composite order.

Host: Intel Xeon @ 2.10 GHz, 4 vCPU cloud container, shared, unpinned, L0
wall clock.
