"""Statistics shared by the programme's rounds (IC_TOOL_PROGRAM.md §5).

The primary interval is §22's: the geometric mean of paired ratios with
a Student `t` interval on their logarithms.  The `t` quantile is
computed here, from the regularised incomplete beta function, so that
no round depends on a table's range or on SciPy.
"""
from __future__ import annotations

import math
import statistics


def _betacf(a: float, b: float, x: float) -> float:
    """Continued fraction for the incomplete beta (Numerical Recipes, betacf)."""
    tiny, eps = 1e-300, 3e-16
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c, d = 1.0, 1.0 - qab * x / qap
    d = 1.0 / (d if abs(d) > tiny else tiny)
    h = d
    for m in range(1, 400):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        d = 1.0 / (d if abs(d) > tiny else tiny)
        c = 1.0 + aa / c
        c = c if abs(c) > tiny else tiny
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        d = 1.0 / (d if abs(d) > tiny else tiny)
        c = 1.0 + aa / c
        c = c if abs(c) > tiny else tiny
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < eps:
            break
    return h


def betai(a: float, b: float, x: float) -> float:
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    front = math.exp(math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b) + a * math.log(x) + b * math.log1p(-x))
    if x < (a + 1.0) / (a + b + 2.0):
        return front * _betacf(a, b, x) / a
    return 1.0 - front * _betacf(b, a, 1.0 - x) / b


def t_cdf(t: float, dof: int) -> float:
    x = dof / (dof + t * t)
    tail = 0.5 * betai(dof / 2.0, 0.5, x)
    return 1.0 - tail if t >= 0 else tail


def t_quantile(p: float, dof: int) -> float:
    """The `p` quantile of Student's t with `dof` degrees of freedom, by bisection."""
    lo, hi = -1e3, 1e3
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if t_cdf(mid, dof) < p:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def geo_ci(ratios: list[float], level: float = 0.95) -> dict:
    """Geometric mean of ratios, with a `t` interval on the logarithms."""
    logs = [math.log(x) for x in ratios]
    out: dict = {"n": len(logs)}
    if not logs:
        return out
    m = statistics.fmean(logs)
    out["geomean"] = math.exp(m)
    out["min"], out["max"] = min(ratios), max(ratios)
    if len(logs) > 1:
        h = t_quantile(0.5 + level / 2, len(logs) - 1) * statistics.stdev(logs) / math.sqrt(len(logs))
        out.update({"lo": math.exp(m - h), "hi": math.exp(m + h)})
    return out


def fit(xs: list[float], ys: list[float], level: float = 0.95) -> dict | None:
    """Least squares with a `t` interval on the slope; None below three points."""
    n = len(xs)
    if n < 3:
        return None
    xm, ym = statistics.fmean(xs), statistics.fmean(ys)
    sxx = sum((x - xm) ** 2 for x in xs)
    beta = sum((x - xm) * (y - ym) for x, y in zip(xs, ys)) / sxx
    alpha = ym - beta * xm
    dof = n - 2
    s2 = sum((y - alpha - beta * x) ** 2 for x, y in zip(xs, ys)) / dof
    h = t_quantile(0.5 + level / 2, dof) * math.sqrt(s2 / sxx)
    return {"beta": beta, "lo": beta - h, "hi": beta + h, "alpha": alpha, "points": n, "dof": dof}


def median_rep(rep: dict, key: str) -> float:
    return statistics.median(r[key] for r in rep["repetitions"])


def ic_cold_ns(rep: dict) -> float:
    """The index calculus's cold time: set-up plus the online interval, per
    repetition, at the median repetition (IC_TOOL_PROGRAM.md §5)."""
    return statistics.median(r["setup_ns"] + r["ic_online"]["wall_ns"] for r in rep["repetitions"])


def ic_online_ns(rep: dict) -> float:
    return statistics.median(r["ic_online"]["wall_ns"] for r in rep["repetitions"])


def rho_online_ns(rep: dict) -> float:
    return statistics.median(r["rho_online"]["wall_ns"] for r in rep["repetitions"])


def setup_phase_ns(rep: dict) -> dict:
    keys = rep["repetitions"][0]["setup_phases_ns"].keys()
    return {k: statistics.median(r["setup_phases_ns"][k] for r in rep["repetitions"]) for k in keys}


def online_phase_ns(rep: dict) -> dict:
    keys = rep["repetitions"][0]["ic_online"]["phases_ns"].keys()
    return {k: statistics.median(r["ic_online"]["phases_ns"][k] for r in rep["repetitions"]) for k in keys}
