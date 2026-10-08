#!/usr/bin/env python3
"""Ledger §23's declared prediction, computed before any §23 run.

Everything here comes from frozen files and formulas:

- the index calculus's per-target descent and reusable set-up, in units,
  from §22's isolated candidate records (`research/ic_descent_20260930/
  runs-isolated/main`, the median over the 20 files of each size);
- at the two sizes §22 did not run, a log-log interpolation in r between
  the two neighbouring sizes, marked as such;
- rho's floor `E = sqrt(pi r / 4n)` steps, priced at two step costs: the
  canonical step the thread has used since §19 (one batched addition and
  the table canonicalisation, 2.30 units at n = 41 in the declaration's
  probe) and the step the declaration's probe measured for the batched
  one-target walk (4.3 units at n = 41);
- Bernstein and Lange's precomputation trade-off (2012) as a model
  boundary: a generic walk that spends P steps building a table of
  distinguished points solves a new target in O steps with
  O * P ~ 1.93 * 1.21 * r / A, here with A = 2n and P the index
  calculus's own reusable set-up converted at the canonical step.

    python3 predict.py > prediction.json
"""
from __future__ import annotations

import glob
import json
import math
import statistics
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
S22 = ROOT / "research" / "ic_descent_20260930" / "runs-isolated" / "main"

# The six sizes and their recipes (PROTOCOL.md, "Recipe").
SIZES = [
    # (a, n, r, columns, descent summands, recipe source)
    (1, 47, 106781081677, 56, 3, "§20 sweep"),
    (0, 57, 275295876199, 64, 2, "§20 model optimum"),
    (0, 41, 549756390943, 80, 2, "§20 sweep"),
    (0, 53, 21044858204113, 264, 2, "§20 sweep"),
    (1, 59, 25179555920633, 256, 2, "§20 model optimum"),
    (0, 61, 162888033982417, 320, 2, "§20 sweep"),
]
CANONICAL_STEP = 2.30  # units; the declaration's probe, n = 41 (§22: 2.27)
PROBE_STEP = 4.3  # units; the declaration's probe of the batched walk, n = 41
BL_PRODUCT = 1.93 * 1.21  # Bernstein-Lange: online * precomputation / (r/A)


def floor_steps(r: int, n: int) -> float:
    return math.sqrt(math.pi * r / (4 * n))


def s22(a: int, n: int) -> dict | None:
    """Median per-target descent and reusable set-up over §22's files."""
    rows = []
    for path in sorted(glob.glob(str(S22 / f"k{a}n{n}" / "M*" / "r*-candidate.price.json"))):
        d = json.loads(Path(path).read_text())
        if d.get("status") != "complete":
            continue
        ph, k = d["median"]["phases_units"], d["targets"]
        online = (ph["descent"] + ph["verify_final"]) / k
        setup = sum(v for key, v in ph.items() if key not in ("descent", "verify_final"))
        rows.append((online, setup))
    if not rows:
        return None
    return {
        "files": len(rows),
        "ic_online_units": statistics.median(x[0] for x in rows),
        "setup_units": statistics.median(x[1] for x in rows),
    }


def interpolate(r: int, lo: tuple[int, dict], hi: tuple[int, dict], key: str) -> float:
    """Log-log interpolation of key/sqrt(r) between two sizes."""
    (r0, d0), (r1, d1) = lo, hi
    y0, y1 = math.log(d0[key] / math.sqrt(r0)), math.log(d1[key] / math.sqrt(r1))
    t = (math.log(r) - math.log(r0)) / (math.log(r1) - math.log(r0))
    return math.exp(y0 + t * (y1 - y0)) * math.sqrt(r)


def main() -> None:
    known = {(a, n): s22(a, n) for a, n, *_ in SIZES}
    rows = []
    for a, n, r, columns, m, source in SIZES:
        base = known[(a, n)]
        if base is None:
            # The neighbouring sizes in r, below and above.
            below = max((rr, known[(aa, nn)]) for aa, nn, rr, *_ in SIZES if rr < r and known[(aa, nn)])
            above = min((rr, known[(aa, nn)]) for aa, nn, rr, *_ in SIZES if rr > r and known[(aa, nn)])
            ic = interpolate(r, below, above, "ic_online_units")
            setup = interpolate(r, below, above, "setup_units")
            how = "interpolated in log-log between the neighbouring §22 sizes"
        else:
            ic, setup = base["ic_online_units"], base["setup_units"]
            how = f"§22 isolated candidate records, median of {base['files']} files"
        e = floor_steps(r, n)
        sqrt_r = math.sqrt(r)
        rho_probe = e * PROBE_STEP
        rho_model = e * CANONICAL_STEP
        a_order = 2 * n
        bl_online = BL_PRODUCT * r * CANONICAL_STEP**2 / (a_order * setup)
        rows.append({
            "curve": f"K_{a}/GF(2^{n})", "a": a, "n": n, "r": r, "log2_r": round(math.log2(r), 2),
            "recipe": {"columns": columns, "descent_summands": m, "source": source},
            "ic_source": how,
            "floor_steps": round(e),
            "s_floor_steps": math.sqrt(math.pi / (4 * n)),
            "ic_online_units": round(ic),
            "s_ic_online": ic / sqrt_r,
            "setup_units": round(setup),
            "s_setup": setup / sqrt_r,
            "rho_online_units_probe_step": round(rho_probe),
            "rho_online_units_canonical_step": round(rho_model),
            "online_speedup_probe_step": rho_probe / ic,
            "online_speedup_canonical_step": rho_model / ic,
            "cold_ratio_probe_step": (setup + ic) / rho_probe,
            "break_even_targets_probe_step": setup / (rho_probe - ic) if rho_probe > ic else None,
            "bl_online_units_model": round(bl_online),
            "ic_online_over_bl_model": ic / bl_online,
        })
    doc = {
        "what_this_is": "Ledger §23's declared prediction, computed from frozen files and formulas before any §23 run.",
        "canonical_step_units": CANONICAL_STEP,
        "probe_step_units": PROBE_STEP,
        "bernstein_lange_product": BL_PRODUCT,
        "rows": rows,
    }
    print(json.dumps(doc, indent=1))


if __name__ == "__main__":
    main()
