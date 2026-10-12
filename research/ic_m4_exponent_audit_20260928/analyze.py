#!/usr/bin/env python3
"""Pre-registered analysis of the m = 4 exponent audit (PREREGISTRATION.md §4-§7).

    python3 research/ic_m4_exponent_audit_20260928/analyze.py [RUNS_DIR]

Reads RUNS_DIR/{semaev,enumerate,null,degree}/K{a}n{n}l{ell}.jsonl (default: runs/ next to
this file), written by examples/m4_exponent_audit.rs, and prints the registered readout.
Pure Python (no numpy), deterministic: the bootstrap uses random.Random(BOOT_SEED).

Every rule here is fixed in PREREGISTRATION.md before any audit cell ran; change it only by
an additive, dated amendment there.
"""
import json
import math
import pathlib
import random
import statistics
import sys

HERE = pathlib.Path(__file__).resolve().parent
CURVES = (0, 1)
# Cells the survey names that do not exist in the tooling: KoblitzCurve::new admits a curve
# only when #E = h*r with r prime, r^2 not dividing #E, and r > h (PREREGISTRATION §1).
EXCLUDED = {(0, 11): "r = 23 <= h = 92 (23^2 | #E)", (0, 17): "r = 239 <= h = 548",
            (1, 13): "r = 79 <= h = 106"}
PRIME_N = (11, 13, 17, 19)
SIZES = (9, 11, 13, 15, 17, 19)
DEGREE_SIZES = (9, 11, 13)
TARGETS = {"semaev": 16, "enumerate": 16, "null": 8, "degree": 16}
C_STAR = 0.25          # survey §2 / §5.4: c*(m = 4)
S_STAR = 0.05          # survey §2: degree slope keeping c < c* in the H5 model
# Full-tree point additions of the enumeration null are a function of |F| alone; with the
# registered subspaces' realised |F| (PREREGISTRATION §1 table) their LS slope on n is 0.828.
ENUM_MODEL = 0.828
ENUM_RANGE = (0.60, 1.00)
B = 10_000
BOOT_SEED = 20_260_928
INF = math.inf


def ell_of(n):
    return (n + 2) // 4


def cell_name(a, n):
    return f"K{a}n{n}l{ell_of(n)}"


def load(runs, arm, a, n):
    path = runs / arm / f"{cell_name(a, n)}.jsonl"
    if not path.exists():
        return None
    rows = {}
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        d = json.loads(line)
        assert d["arm"] == arm and d["a"] == a and d["n"] == n, path
        assert d["label"] != "smoke", f"smoke line in {path}"
        rows[d["target"]] = d
    return rows


def cost(arm, d):
    """Per-target cost, INF when censored (a lower bound, never a measurement)."""
    if d is None:
        return INF  # never emitted: watchdog/CPU-limit kill -> censored
    if arm == "semaev":
        return INF if d["verdict"] == "censored" else d["word_ops"]
    if arm == "null":
        return INF if d["verdict"] == "censored" else d["word_ops"]
    if arm == "enumerate":
        return d["group_adds"]
    raise ValueError(arm)


def lower_median(values):
    """The ceil(N/2)-th smallest value: finite iff at most half are censored."""
    v = sorted(values)
    return v[(len(v) + 1) // 2 - 1] if v else INF


def ols(points):
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    mx, my = statistics.fmean(xs), statistics.fmean(ys)
    sxx = sum((x - mx) ** 2 for x in xs)
    return sum((x - mx) * (y - my) for x, y in points) / sxx


def fit(cells):
    """cells: {(a, n): [cost, ...]} -> (slope, points, dropped)."""
    points, dropped = [], []
    for (a, n), vals in sorted(cells.items()):
        med = lower_median(vals)
        if med == INF or med <= 0:
            dropped.append((a, n))
        else:
            points.append((n, math.log2(med)))
    if len({p[0] for p in points}) < 3:
        return None, points, dropped
    return ols(points), points, dropped


def bootstrap(cells, rng):
    slopes, invalid = [], 0
    keys = sorted(cells)
    for _ in range(B):
        pts = []
        for k in keys:
            vals = cells[k]
            if not vals:
                continue
            res = [vals[rng.randrange(len(vals))] for _ in vals]
            med = lower_median(res)
            if med != INF and med > 0:
                pts.append((k[1], math.log2(med)))
        if len({p[0] for p in pts}) < 3:
            invalid += 1
            continue
        slopes.append(ols(pts))
    slopes.sort()
    if not slopes:
        return None, None, invalid
    lo = slopes[int(0.025 * len(slopes))]
    hi = slopes[min(len(slopes) - 1, int(math.ceil(0.975 * len(slopes))) - 1)]
    return lo, hi, invalid


def arm_cells(runs, arm, subset=None):
    cells, meta = {}, {}
    for a in CURVES:
        for n in SIZES:
            rows = load(runs, arm, a, n)
            if rows is None:
                continue
            want = range(TARGETS[arm])
            ds = [rows.get(t) for t in want]
            if subset is not None:
                ds = [d for d in ds if d is not None and subset(d)]
            cells[(a, n)] = [cost(arm, d) for d in ds]
            meta[(a, n)] = ds
    return cells, meta


def report_arm(title, cells, rng, *, decisive=False):
    slope, points, dropped = fit(cells)
    print(f"\n## {title}")
    print(f"{'cell':<10} {'targets':>7} {'censored':>8} {'lower-median':>14} {'log2':>7}")
    for (a, n), vals in sorted(cells.items()):
        med = lower_median(vals)
        cen = sum(v == INF for v in vals)
        lg = f"{math.log2(med):7.2f}" if med not in (INF, 0) else "      -"
        print(f"{cell_name(a, n):<10} {len(vals):>7} {cen:>8} {med:>14} {lg}")
    for k in dropped:
        print(f"dropped from fit (more than half censored, or empty): {cell_name(*k)}")
    if slope is None:
        print("slope: not fitted (fewer than 3 distinct n retained)")
        return None, None, None
    lo, hi, invalid = bootstrap(cells, rng)
    print(f"slope c_hat = {slope:.4f} bits per unit n; bootstrap 95% band "
          f"[{lo:.4f}, {hi:.4f}] (B = {B}, invalid replicates {invalid})"
          if lo is not None else f"slope c_hat = {slope:.4f}; bootstrap: no valid replicate")
    if decisive and lo is not None:
        if hi < C_STAR:
            verdict = "ALIVE (upper end below c* = 0.25)"
        elif lo > C_STAR:
            verdict = "CLOSED for this engine at these sizes (lower end above c* = 0.25)"
        else:
            verdict = "INCONCLUSIVE (band contains c* = 0.25)"
        print(f"verdict reading (before control gating): {verdict}")
    return slope, lo, hi


def main():
    runs = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / "runs"
    rng = random.Random(BOOT_SEED)
    print(f"# m = 4 exponent audit readout  (runs: {runs})")

    sem, sem_meta = arm_cells(runs, "semaev")
    if not sem:
        print("no semaev cells found")
        return
    # Instrument checks on the primary arm.
    disagree = [(k, d["target"]) for k, ds in sem_meta.items() for d in ds
                if d is not None and d.get("agree") is False]
    unverified = [(k, d["target"]) for k, ds in sem_meta.items() for d in ds
                  if d is not None and d.get("verified") is False]
    print(f"\noracle/enumeration disagreements: {len(disagree)} {disagree[:10]}")
    print(f"returned decompositions failing the group check: {len(unverified)} {unverified[:10]}")

    for k, why in EXCLUDED.items():
        print(f"structurally excluded (no such curve in the tooling): {cell_name(*k)}: {why}")
    primary = report_arm("Primary: Semaev m = 4, all targets (word XORs)", sem, rng, decisive=True)
    report_arm("Sensitivity (not decisive): Semaev, prime n only",
               {k: v for k, v in sem.items() if k[1] in PRIME_N}, rng, decisive=True)
    ref, _ = arm_cells(runs, "semaev", subset=lambda d: d["verdict"] == "refuted")
    report_arm("Secondary: Semaev, refuted-only (oracle verdict)", ref, rng)
    sat, _ = arm_cells(runs, "semaev", subset=lambda d: d["verdict"] == "satisfiable")
    print("\n## Secondary: Semaev, satisfiable-only (oracle verdict)")
    for (a, n), vals in sorted(sat.items()):
        print(f"{cell_name(a, n):<10} satisfiable {len(vals):>2}  lower-median "
              f"{lower_median(vals) if vals else '-'}")
    if sat and all(len(v) >= 3 for v in sat.values()) and len(sat) >= 3:
        report_arm("Secondary: Semaev, satisfiable-only fit", sat, rng)
    else:
        print("satisfiable-only slope: not fitted (a cell has fewer than 3 satisfiable targets)")

    en, _ = arm_cells(runs, "enumerate")
    e = report_arm("Control: enumeration null (point additions)", en, rng)
    enum_ok = None
    if e[0] is not None:
        ok = enum_ok = ENUM_RANGE[0] <= e[0] <= ENUM_RANGE[1]
        print(f"enumeration control: c_hat_enum = {e[0]:.3f}, model {ENUM_MODEL}, "
              f"registered range {ENUM_RANGE}: {'PASS' if ok else 'FAIL (instrument suspect)'}")

    nu, nu_meta = arm_cells(runs, "null")
    report_arm("Control: random-system null (word XORs)", nu, rng)
    print("\nnull vs Semaev per retained Semaev cell (lower medians):")
    null_status = "pass"
    both_measured = 0
    for k in sorted(sem):
        s = lower_median(sem[k])
        if s == INF:
            continue
        v = nu.get(k)
        if v is None:
            print(f"{cell_name(*k)}: null not run")
            null_status = "incomplete"
            continue
        m = lower_median(v)
        rel = "censored (dearer than its budget)" if m == INF else ("dearer" if m > s else "NOT dearer")
        if m != INF:
            both_measured += 1
        if m != INF and m <= s:
            null_status = "fail"
        print(f"{cell_name(*k)}: semaev {s}  null {m if m != INF else 'censored'}  -> {rel}")
    print(f"null control (dearer at every retained size): {null_status.upper()}; "
          f"sizes where both are measured: {both_measured} "
          f"(growth comparison needs >= 3; otherwise 'growth not measurable within budget')")

    # The registered verdict (PREREGISTRATION.md §5.2-§6), gating applied.
    _, lo, hi = primary
    if unverified:
        verdict = "VOID (a returned decomposition failed the group check, §5.4)"
    elif lo is None:
        verdict = "INCONCLUSIVE (insufficient uncensored sizes)"
    elif hi < C_STAR:
        verdict = "ALIVE"
        if null_status == "fail":
            verdict = "INCONCLUSIVE (null control failed)"
    elif lo > C_STAR:
        verdict = "CLOSED for this engine at these sizes"
    else:
        verdict = "INCONCLUSIVE"
    if enum_ok is False and not verdict.startswith("VOID"):
        verdict = f"INSTRUMENT CHECK FAILED (enumeration control); reading: {verdict}"
    elif enum_ok is None and not verdict.startswith("VOID"):
        verdict = f"{verdict} [enumeration control not evaluable]"
    print(f"\n## REGISTERED VERDICT: {verdict}")

    # Degree secondary.
    print("\n## Secondary: refutation degree D(l) on Boolean-unsatisfiable targets")
    pts, faults = [], 0
    for a in CURVES:
        for n in DEGREE_SIZES:
            rows = load(runs, "degree", a, n)
            if rows is None:
                continue
            outs = []
            for t in range(TARGETS["degree"]):
                d = rows.get(t)
                if d is None:
                    continue
                if d["count_check"] is False:
                    faults += 1
                o = d["outcome"]
                if o["kind"] == "resolved" and o["refuted"]:
                    outs.append(str(o["degree"]))
                    pts.append((ell_of(n), o["degree"]))
                elif o["kind"] == "at_least":
                    outs.append(f">={o['degree']}")
                elif o["kind"] == "caps_hit":
                    outs.append("caps")
                elif o["kind"] == "resolved":
                    outs.append(f"{o['degree']}p")
            nsat = sum(1 for d in rows.values() if d["boolean_solutions"] > 0)
            print(f"{cell_name(a, n):<10} Boolean-sat {nsat:>2}/{len(rows)}  unsat outcomes: {' '.join(outs) or '-'}")
    print(f"exact-count vs solver-root-count faults: {faults} (any fault voids this secondary arm)")
    if len({p[0] for p in pts}) >= 2:
        s = ols(pts)
        print(f"refutation-degree slope s over resolved draws = {s:.3f} (s* ~ {S_STAR}); "
              "lower bounds (>=) are not in the fit")
    else:
        print("refutation-degree slope: not fitted (resolved degrees at fewer than 2 distinct l)")


if __name__ == "__main__":
    main()
