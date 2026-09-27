#!/usr/bin/env python3
"""Tables and fits for the PKM tower pilot.

Reads the JSONL files written by `examples/pkm_tower_pilot.rs` and prints,
per (kind, m, control, p):

- one row per N: runs, timeouts, and the median and range of the solving
  degree (`solving_degree_max`, the framework's definition), of the widest
  matrix and the F4 step count up to the last productive step, and of the F4
  wall time;
- the pre-registered fit of RESEARCH_PKM_TOWER_ORACLE.md section 5.2,
  D = alpha + beta*N by least squares, with a 95% bootstrap interval for beta
  that resamples targets within each N;
- because that interval collapses when D is constant within each N, the
  leave-one-N-out slopes, the slope over the upper half of the N range, and
  the final plateau of D (section 10.8's proposed criterion, post hoc);
- the growth rate of the matrix width, as the least-squares slope of
  log2(max_cols_to_solution) against N.

It exits non-zero if a repeated system disagrees with itself, a planted
solution was lost, or the two engines (`f4_fp` and `f4_fp_tower`) disagree on
whether a system has a solution: each is a bug, not data. Their solving
degrees are compared and printed, not enforced: a different degree is data.

Runs that hit the budget are reported, and excluded from the fits: their
degree is a lower bound, never a value.

Files may overlap. Every cell draws its tower, curve and targets from a seed
fixed by (seed, kind, m, t, g), so a cell that two runs share is the same
system measured twice. Such a system is counted once, as the copy read
first, unless that copy timed out and a later one finished: the finished copy
then stands for the system, and its degree must not fall below the lower
bound of the copy that stopped. Finished copies must agree on every
deterministic field, and any disagreement is printed. A copy halted by the
staircase stop is compared only on the solving degree and the width. Rows
written before `solving_degree_max` existed are skipped and counted.

    python3 analyze.py runs/*.jsonl
"""

import json
import math
import random
import statistics
import sys
from collections import defaultdict

# Fields that depend only on the system and the algorithm, not on the machine.
DETERMINISTIC = (
    "solving_degree_max",
    "last_productive_degree",
    "degree_reached",
    "max_cols_to_solution",
    "steps_to_solution",
    "inconsistent",
    "basis_len",
)
# The fields a run halted by the staircase stop shares with a full run of the
# same system.
STOP_COMPARABLE = ("solving_degree_max", "max_cols_to_solution")


def engine(r):
    """Rows written before the tower engine existed are all `f4_fp`'s."""
    return r.get("engine", "f4_fp")


def instance_key(r):
    """The system itself, whichever engine measured it."""
    return (
        r["p"], r["kind"], r["m"], r["control"], r["t"], r["g"], r["target"],
        r["target_index"], r["x_r"], r["curve"]["a"], r["curve"]["b"],
        json.dumps(r["tower"], sort_keys=True),
    )


def default_bound(r):
    """The example's degree bound without `--cap`: `n + d + 6`, or `2d + 8`
    for the naive control."""
    d = max(r["input_degrees"])
    return 2 * d + 8 if r["control"] == "naive" else r["n_vars"] + d + 6


def capped(r):
    """A run with a degree bound below the default: a confirmation run (note
    section 11.4), which is expected to stop at the bound and is not a
    measurement of D."""
    bound = r.get("max_degree_bound")
    return bound is not None and bound < default_bound(r)


def system_key(r):
    """One measurement: a system, the engine that measured it, and the
    degree bound it ran under."""
    return instance_key(r) + (engine(r), r.get("max_degree_bound"))


def load(paths, report=True):
    rows, seen = [], {}
    duplicates, mismatches, old_format = 0, [], 0
    for path in paths:
        with open(path) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                r = json.loads(line)
                if "summary" in r:
                    continue
                if "solving_degree_max" not in r:
                    old_format += 1
                    continue
                k = system_key(r)
                if k in seen:
                    duplicates += 1
                    first = seen[k]
                    # A timed-out copy stopped early; only finished copies compare.
                    if not (r["timed_out"] or first["timed_out"]):
                        # A copy halted by the staircase stop (`--stop-below`) ends
                        # before the basis is certified. It shares with a full run
                        # only the degree and width reached by then.
                        stopped = r.get("staircase_at_stop") is not None or first.get(
                            "staircase_at_stop"
                        ) is not None
                        fields = STOP_COMPARABLE if stopped else DETERMINISTIC
                        diff = [f for f in fields if r.get(f) != first.get(f)]
                        if diff:
                            mismatches.append((path, r["kind"], r["m"], r["control"], r["N"], diff))
                    elif r["timed_out"] != first["timed_out"]:
                        # One copy finished and one stopped early, whose degree
                        # is only a lower bound: the finished copy stands for the
                        # system (round 3 finished what round 2's D2 could not),
                        # and its degree must not fall below that bound.
                        done, cut = (first, r) if r["timed_out"] else (r, first)
                        if done["solving_degree_max"] < cut["solving_degree_max"]:
                            mismatches.append(
                                (path, r["kind"], r["m"], r["control"], r["N"],
                                 ["solving_degree_max below a stopped copy's lower bound"])
                            )
                        if done is r:
                            rows[rows.index(first)] = r
                            seen[k] = r
                    continue
                seen[k] = r
                rows.append(r)
    if report:
        print(
            f"{len(rows)} distinct systems; {duplicates} repeated measurements "
            f"({len(mismatches)} disagreeing on a deterministic field); "
            f"{old_format} rows in the format before `solving_degree_max`, skipped."
        )
        for m in mismatches:
            print("MISMATCH", m)
    return rows, len(mismatches)


def lsq(xs, ys):
    n = len(xs)
    mx = sum(xs) / n
    my = sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    if sxx == 0:
        return None
    b = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sxx
    return my - b * mx, b


def bootstrap_slope(by_n, key, reps=4000, seed=20260924):
    """Resample targets within each N; refit the slope each time."""
    rng = random.Random(seed)
    ns = sorted(by_n)
    if len(ns) < 2:
        return None
    slopes = []
    for _ in range(reps):
        xs, ys = [], []
        for n in ns:
            vals = by_n[n]
            for _ in range(len(vals)):
                xs.append(n)
                ys.append(key(rng.choice(vals)))
        fit = lsq(xs, ys)
        if fit is not None:
            slopes.append(fit[1])
    slopes.sort()
    lo = slopes[int(0.025 * len(slopes))]
    hi = slopes[int(0.975 * len(slopes)) - 1]
    return lo, hi


def jackknife_slopes(by_n, key):
    """Least-squares slopes with each N left out in turn, over per-N medians.

    The pre-registered bootstrap resamples targets within each N. When every
    target at an N gives the same degree, as the solving degree mostly does,
    that interval collapses to a point and says nothing about the shape of
    D(N). Leaving whole N values out does show it.
    """
    ns = sorted(by_n)
    if len(ns) < 3:
        return None
    med = {n: statistics.median(key(r) for r in by_n[n]) for n in ns}
    out = []
    for drop in ns:
        keep = [n for n in ns if n != drop]
        fit = lsq(keep, [med[n] for n in keep])
        if fit is not None:
            out.append(fit[1])
    return min(out), max(out)


def upper_half_slope(by_n, key):
    ns = sorted(by_n)
    if len(ns) < 4:
        return None
    upper = ns[len(ns) // 2 :]
    fit = lsq(upper, [statistics.median(key(r) for r in by_n[n]) for n in upper])
    return fit[1] if fit else None, upper[0], upper[-1]


def final_plateau(by_n, key):
    """The N range since D last rose, for the plateau criterion proposed in
    section 10.8 of the note (post hoc; for Stage A, not a re-grading)."""
    ns = sorted(by_n)
    med = [statistics.median(key(r) for r in by_n[n]) for n in ns]
    start = ns[0]
    for i in range(1, len(ns)):
        if med[i] > med[i - 1]:
            start = ns[i]
    return start, ns[-1], med[-1]


def a1_reading(by_n, key):
    """Amendment A1 of the note (section 10.8), as section 11.6 adopted it for
    round 2: H1a if the final plateau is longer than 10 in N; H0 if D rises at
    least every 4 in N over the upper half of the range, read as a rate (the
    least-squares slope over the upper half is at least 1/4, the rate of
    section 5.2's H0 threshold); inconclusive otherwise."""
    ns = sorted(by_n)
    lo_n, hi_n, _ = final_plateau(by_n, key)
    if hi_n - lo_n > 10:
        return f"H1a (final plateau L = {hi_n - lo_n} > 10)"
    up = upper_half_slope(by_n, key)
    if up is None or up[0] is None:
        return f"inconclusive (only {len(ns)} values of N)"
    reading = "H0" if up[0] >= 0.25 else "inconclusive"
    return f"{reading} (upper-half slope {up[0]:.3f} over N = {up[1]}…{up[2]}, final plateau L = {hi_n - lo_n})"


def compare_engines(rows):
    """Systems both engines finished: do their solving degrees agree?

    Returns the number of systems on which the verdicts (refuted or not)
    disagree. That is a bug in one engine; a different degree is data."""
    by_instance = defaultdict(dict)
    for r in rows:
        if not r["timed_out"]:
            by_instance[instance_key(r)][engine(r)] = r
    pairs = [(v["f4_fp"], v["f4_fp_tower"]) for v in by_instance.values()
             if "f4_fp" in v and "f4_fp_tower" in v]
    if not pairs:
        return 0
    agree = [a for a, b in pairs if a["solving_degree_max"] == b["solving_degree_max"]]
    verdicts = [a for a, b in pairs if a["inconsistent"] == b["inconsistent"]]
    print("\n### The two engines on the same systems\n")
    print(
        f"{len(pairs)} systems finished by both: the solving degree agrees on "
        f"{len(agree)}; the verdict (refuted or not) agrees on {len(verdicts)}."
    )
    diff = defaultdict(int)
    for a, b in pairs:
        if a["solving_degree_max"] != b["solving_degree_max"]:
            diff[(a["kind"], a["m"], a["control"], a["N"], a["target"],
                  a["solving_degree_max"], b["solving_degree_max"])] += 1
    if diff:
        print("\n| kind | m | control | N | target | D f4_fp | D f4_fp_tower | systems |")
        print("|:--|--:|:--|--:|:--|--:|--:|--:|")
        for (kind, m, control, n, target, d1, d2), c in sorted(diff.items()):
            print(f"| {kind} | {m} | {control} | {n} | {target} | {d1} | {d2} | {c} |")
    split = [(a, b) for a, b in pairs if a["inconsistent"] != b["inconsistent"]]
    for a, b in split:
        print("VERDICT MISMATCH", a["kind"], a["m"], a["control"], a["N"], a["target"],
              a["target_index"], "f4_fp:", a["inconsistent"], "f4_fp_tower:", b["inconsistent"])
    return len(split)


def verdict(lo, hi):
    """The decision rule of section 5.2, on beta alone (H1 also needs the null)."""
    if lo is None:
        return "no fit"
    if lo > 0.25:
        return "H0 (beta lower bound > 0.25)"
    if hi < 0.10:
        return "H1 candidate (beta upper bound < 0.10; also needs D below the null)"
    return "inconclusive"


def confirmations(capped_rows, rows):
    """The confirmation runs of note section 11.4, beside the full runs of the
    same systems. Returns how many contradict their full run: a capped run
    that refutes although the full run needed a step above the bound. The
    algorithm is deterministic and the steps up to the bound are the same, so
    that is a bug, not data."""
    if not capped_rows:
        return 0
    full = {(instance_key(r), engine(r)): r for r in rows}
    print("\n### Confirmation runs (degree bound below the default, section 11.4)\n")
    print("| kind | m | control | p | N | target | bound | D, full run | pairs above the bound | refuted | confirms |")
    print("|:--|--:|:--|--:|--:|:--|--:|--:|--:|:--|:--|")
    bad = 0
    for r in sorted(capped_rows, key=lambda r: (r["kind"], r["m"], r["control"], r["p"], r["N"], r["target_index"])):
        f = full.get((instance_key(r), engine(r)))
        d_full = f["solving_degree_max"] if f and not f["timed_out"] else None
        if r["timed_out"]:
            status = "no (timed out)"
        elif r["inconsistent"]:
            status = "no: refuted under the bound"
            if d_full is not None and d_full > r["max_degree_bound"]:
                bad += 1
                status += " (CONTRADICTS the full run)"
        elif r.get("staircase_at_stop") is not None:
            status = "no (staircase stop)"
        elif r["pairs_above_bound"] > 0:
            status = "yes"
        else:
            status = "no (finished below the bound)"
        print(f"| {r['kind']} | {r['m']} | {r['control']} | {r['p']} | {r['N']} | {r['target']} {r['target_index']} "
              f"| {r['max_degree_bound']} | {d_full if d_full is not None else '—'} | {r['pairs_above_bound']} "
              f"| {'yes' if r['inconsistent'] else 'no'} | {status} |")
    return bad


def main(paths):
    rows, mismatches = load(paths)
    capped_rows = [r for r in rows if capped(r)]
    rows = [r for r in rows if not capped(r)]
    groups = defaultdict(list)
    for r in rows:
        if r["control"] == "ladder":
            continue
        groups[(engine(r), r["kind"], r["m"], r["control"], r["p"])].append(r)
    for (eng, kind, m, control, prime), rs in sorted(groups.items()):
        print(f"\n### {kind}, m = {m}, {control}, p = {prime}, {eng}\n")
        print("| N | runs | timeouts | D (median, range) | width to solution (median) | F4 steps to solution (median) | F4 ms (median) |")
        print("|--:|--:|--:|:--|--:|--:|--:|")
        by_n = defaultdict(list)
        for r in rs:
            by_n[r["N"]].append(r)
        fit_by_n = {}
        for n in sorted(by_n):
            vals = by_n[n]
            done = [r for r in vals if not r["timed_out"]]
            to = len(vals) - len(done)
            if done:
                ds = [r["solving_degree_max"] for r in done]
                ws = [r["max_cols_to_solution"] for r in done]
                ss = [r["steps_to_solution"] for r in done]
                ms = [r["ms"] for r in done]
                d_cell = f"{statistics.median(ds):g} ({min(ds)}–{max(ds)})"
                w_cell = f"{statistics.median(ws):g}"
                s_cell = f"{statistics.median(ss):g}"
                t_cell = f"{statistics.median(ms):.1f}"
                fit_by_n[n] = done
            else:
                lower = max(r["solving_degree_max"] for r in vals)
                d_cell = f"≥ {lower} (all timed out)"
                w_cell = s_cell = t_cell = "—"
            print(f"| {n} | {len(vals)} | {to} | {d_cell} | {w_cell} | {s_cell} | {t_cell} |")
        pts = [(n, r["solving_degree_max"]) for n, vs in fit_by_n.items() for r in vs]
        if len({n for n, _ in pts}) >= 2:
            a, b = lsq([p[0] for p in pts], [p[1] for p in pts])
            ci = bootstrap_slope(fit_by_n, lambda r: r["solving_degree_max"])
            wpts = [(n, math.log2(max(1, r["max_cols_to_solution"]))) for n, vs in fit_by_n.items() for r in vs]
            _, wb = lsq([p[0] for p in wpts], [p[1] for p in wpts])
            wci = bootstrap_slope(fit_by_n, lambda r: math.log2(max(1, r["max_cols_to_solution"])))
            degenerate = ci[0] == ci[1]
            print(
                f"\nD = {a:.2f} + {b:.3f}·N over N = {min(fit_by_n)}…{max(fit_by_n)} "
                f"({len(fit_by_n)} values). Pre-registered bootstrap (targets within N): "
                f"beta 95% [{ci[0]:.3f}, {ci[1]:.3f}]"
                + (" (degenerate: no variance within any N)" if degenerate else "")
                + f"; rule on it: {verdict(*ci)}."
            )
            jk = jackknife_slopes(fit_by_n, lambda r: r["solving_degree_max"])
            if jk:
                print(f"Leave-one-N-out slopes of D: {jk[0]:.3f} to {jk[1]:.3f}.")
            up = upper_half_slope(fit_by_n, lambda r: r["solving_degree_max"])
            if up and up[0] is not None:
                print(f"Slope of D over the upper half, N = {up[1]}…{up[2]}: {up[0]:.3f}.")
            lo_n, hi_n, d_last = final_plateau(fit_by_n, lambda r: r["solving_degree_max"])
            print(f"Final plateau: D = {d_last:g} over N = {lo_n}…{hi_n}, length L = {hi_n - lo_n}.")
            print(f"A1 (section 11.6): {a1_reading(fit_by_n, lambda r: r['solving_degree_max'])}.")
            print(
                f"log2(width) slope {wb:.3f} per unit N over the whole range"
                + (f"; {upper_half_slope(fit_by_n, lambda r: math.log2(max(1, r['max_cols_to_solution'])))[0]:.3f} over the upper half." if up else ".")
            )
    split_verdicts = compare_engines(rows)
    contradicted = confirmations(capped_rows, rows)
    ladder = [r for r in rows if r["control"] == "ladder"]
    if ladder:
        print("\n### Generator-count ladder (planted targets)\n")
        print("| kind | N | g | runs | D (median, range) | width to solution (median) |")
        print("|:--|--:|--:|--:|:--|--:|")
        cells = defaultdict(list)
        for r in ladder:
            cells[(r["kind"], r["N"], r["g"])].append(r)
        for (kind, n, g), rs in sorted(cells.items()):
            done = [r for r in rs if not r["timed_out"]]
            if done:
                ds = [r["solving_degree_max"] for r in done]
                ws = [r["max_cols_to_solution"] for r in done]
                print(f"| {kind} | {n} | {g} | {len(rs)} | {statistics.median(ds):g} ({min(ds)}–{max(ds)}) | {statistics.median(ws):g} |")
            else:
                print(f"| {kind} | {n} | {g} | {len(rs)} | all timed out | — |")
    bad = [r for r in rows if r.get("planted_ok") is False]
    print(f"\nPlanted-solution violations: {len(bad)} of {len(rows)} rows.")
    # A repeat that disagrees, a lost planted solution, two engines that
    # disagree on whether a system has a solution, or a confirmation run that
    # contradicts its full run: each is a bug, not data.
    return 1 if mismatches or bad or split_verdicts or contradicted else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
