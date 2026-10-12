#!/usr/bin/env python3
"""Erratum 1 to DECOMPOSITION-SURVEY.md: the double-large-prime exponent next to the
survey's no-large-prime model.  Pure arithmetic; derived, not measured.

  E_base(c, m) = 2(1 + c)/(m + 1)  for c <= 1/m, else 1/m + c   (survey section 2)
  E_dlp(c, m)  = 2(m - 1)/m^2 + c                              (Theorem 3 carried over, no rebalancing)

E_dlp at c = 0 is the q^{2-2/m} of the double-large-prime variant written in the unit
r = q^m.  Adding c without re-optimising the base is conservative: rebalancing can
only lower it.
"""


def e_base(c, m):
    return 2 * (1 + c) / (m + 1) if c <= 1 / m else 1 / m + c


def e_dlp(c, m):
    return 2 * (m - 1) / m ** 2 + c


def bar(f, m):
    """Largest c with f(c, m) < 1/2, by bisection on [0, 1]; 0.0 if none."""
    if f(0.0, m) >= 0.5:
        return 0.0
    lo, hi = 0.0, 1.0
    for _ in range(60):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if f(mid, m) < 0.5 else (lo, mid)
    return lo


print(f"{'m':>3} {'E_base(0,m)':>12} {'E_dlp(0,m)':>11} {'c* base':>8} {'c* dlp':>7} {'best bar':>9}  which")
for m in (2, 3, 4, 5, 6, 8):
    cb, cd = bar(e_base, m), bar(e_dlp, m)
    print(f"{m:>3} {e_base(0, m):>12.4f} {e_dlp(0, m):>11.4f} {cb:>8.4f} {cd:>7.4f} {max(cb, cd):>9.4f}  "
          f"{'dlp' if cd > cb else 'base'}")
print()
print("m = 3: E_dlp(0,3) = 4/9 < 1/2, so a free oracle does NOT tie rho at best once two large")
print("primes are allowed; the bar becomes c < 1/18 ~ 0.0556.  For m >= 4 the double-large-prime")
print("bar (1/2 - 2(m-1)/m^2) is below the survey's c*(m), so it loosens no m >= 4 gate.")
