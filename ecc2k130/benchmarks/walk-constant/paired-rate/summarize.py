#!/usr/bin/env python3
"""The paired sigma / table (rule v2) rate, and what it makes the cost to solve.

    python3 benchmarks/walk-constant/paired-rate/summarize.py RUN_DIR [--freeze DIR]
    python3 benchmarks/walk-constant/paired-rate/summarize.py --check-deadline

RUN_DIR is what a launcher brought back from gpujob.sh (its results/ files
directly, or a results/ subdirectory).  Per geometry this prints the median
rate of each binary, the sigma/table ratio of every repetition (the two runs
of one round are the pair), and the cost to solve, table over sigma:

    cost = (sigma rate / table rate) x K,
    K    = c_table (1 + delta/2) / c_sigma x table trap loss / sigma guard loss,

with every factor of K taken from the round-2 files by
benchmarks/walk-constant/summarize.py (WALK-CONSTANT.md section 11, item 6),
not copied here.  The verdict follows the rule WALK-CONSTANT.md section 11.5
declared before any run: the table walk pays if the cost bracket's upper end
is below 1 in both geometries, does not if its lower end is above 1 in
either, and is undecided otherwise.  A failed build or verification, or a
bench sample that printed no rate, or a nonzero job exit, voids it.  It also prints the switch
deadline: switching walks abandons the sigma corpus, so a table walk that is
cheaper per solve pays only while that corpus holds less than a fraction f of
sigma's expected work (section 11.5).  --freeze copies the raw files and
summary.json into DIR.
"""
import contextlib
import io
import json
import math
import pathlib
import re
import shutil
import statistics
import sys

HERE = pathlib.Path(__file__).resolve().parent
BINARIES = ("sigma", "table")
GEOMETRIES = ("385024", "auto")
# Section 6's paired sessions for the rule-v1 table kernel, sigma / table.
OLD_RATIO = {"auto": 14.41 / 16.56, "385024": 14.98 / 16.35}


def cost_factor():
    """(K at the sigma bracket's upper end, K at its lower end), from the frozen files."""
    sys.path.insert(0, str(HERE.parent))
    cwd = pathlib.Path.cwd()
    try:
        import os
        os.chdir(HERE.parent)
        with contextlib.redirect_stdout(io.StringIO()):
            import summarize as S
    finally:
        os.chdir(cwd)
    base = S.c_table_v2 * (1 + S.delta_v2 / 2) * S.loss_table_v2 / S.loss_sigma
    return base / S.sigma_hi, base / S.sigma_lo


def deadline(cost):
    """The fraction of sigma's expected work after which a switch at this cost no longer pays.

    With collisions arriving as in a birthday walk, the work W to the first
    one is Rayleigh distributed.  After w0 without one, the expected work
    left is E[W] exp(x^2) erfc(x), x = w0 / (sqrt 2 s), and w0 / E[W] = 2x /
    sqrt(pi).  A fresh table walk costs cost x E[W], so it pays while
    exp(x^2) erfc(x) > cost."""
    if cost >= 1:
        return 0.0
    lo, hi = 0.0, 20.0
    for _ in range(100):
        mid = (lo + hi) / 2
        if math.exp(mid * mid) * math.erfc(mid) > cost:
            lo = mid
        else:
            hi = mid
    return 2 * lo / math.sqrt(math.pi)


def parse_bench(text):
    """({geometry: {binary: [(rep, rate B/s, MHz, W, C)]}}, [samples that printed no rate])."""
    out = {}
    missing = []
    current = None
    for line in text.splitlines():
        m = re.match(r"=== bench (\S+) (\S+) rep (\d+)", line)
        if m:
            if current:
                missing.append(current)
            current = (m.group(1), m.group(2), int(m.group(3)))
            continue
        m = re.search(r"finished: ([\d.]+) M it/s.*?(\d+) MHz, ([\d.]+) W, (\d+)", line)
        if m and current:
            geometry, binary, rep = current
            out.setdefault(geometry, {}).setdefault(binary, []).append(
                (rep, float(m.group(1)) / 1000, int(m.group(2)), float(m.group(3)), int(m.group(4))))
            current = None
    if current:
        missing.append(current)
    return out, missing


def check_deadline(draws=2_000_000):
    """Monte Carlo against deadline(): Rayleigh work, the mean left past f E[W] must equal the cost."""
    import random
    rng = random.Random(1)
    mean = math.sqrt(math.pi / 2)
    work = [math.sqrt(-2 * math.log(1 - rng.random())) for _ in range(draws)]
    worst = 0.0
    for cost in (0.811, 0.863, 0.95):
        w0 = deadline(cost) * mean
        left = [w - w0 for w in work if w > w0]
        got = sum(left) / len(left) / mean
        worst = max(worst, abs(got - cost))
        print("cost %.3f: deadline %.4f of E[W]; simulated work left there %.4f of E[W]" % (cost, deadline(cost), got))
    return 0 if worst < 0.002 else 1


def main(argv):
    if argv == ["--check-deadline"]:
        return check_deadline()
    if not argv or argv[0].startswith("-"):
        print(__doc__)
        return 2
    run = pathlib.Path(argv[0])
    res = run / "results" if (run / "results").is_dir() else run
    freeze = pathlib.Path(argv[argv.index("--freeze") + 1]) if "--freeze" in argv else None
    failures = (res / "failures.txt").read_text().split("\n") if (res / "failures.txt").exists() else []
    failures = [f for f in failures if f.strip()]
    # Every launcher records the job's exit status beside the results.
    if (res / "exit-code").exists() and (res / "exit-code").read_text().strip() != "0":
        failures.append("job exit %s" % (res / "exit-code").read_text().strip())
    verified = {b: "(300 verified against the reference, 0 dropped)" in
                ((res / ("verify-%s.log" % b)).read_text() if (res / ("verify-%s.log" % b)).exists() else "")
                and "MISMATCH" not in (res / ("verify-%s.log" % b)).read_text()
                for b in BINARIES}
    bench, missing = parse_bench((res / "bench.txt").read_text()) if (res / "bench.txt").exists() else ({}, [])
    k_lo, k_hi = cost_factor()
    summary = {"host": (res / "host.txt").read_text().strip() if (res / "host.txt").exists() else None,
               "verified": verified, "failures": failures, "missing_samples": ["%s %s rep %d" % m for m in missing],
               "K": [k_lo, k_hi], "geometries": {}}
    print("K (cost per unit sigma/table ratio) = %.4f - %.4f, from the round-2 files" % (k_lo, k_hi))
    verdicts = []
    for geometry in GEOMETRIES:
        rows = bench.get(geometry, {})
        reps = sorted(set(r[0] for r in rows.get("sigma", [])) & set(r[0] for r in rows.get("table", [])))
        by = {b: {r[0]: r for r in rows.get(b, [])} for b in BINARIES}
        ratios = [by["sigma"][i][1] / by["table"][i][1] for i in reps]
        if not ratios:
            continue
        med = {b: statistics.median(r[1] for r in rows[b]) for b in BINARIES}
        lo, hi = min(ratios) * k_lo, max(ratios) * k_hi
        verdict = "pays" if hi < 1 else "does not pay" if lo > 1 else "undecided"
        verdicts.append(verdict)
        old = OLD_RATIO.get(geometry)
        rule_cost = statistics.median(ratios) / old - 1 if old else None
        summary["geometries"][geometry] = {
            "reps": len(reps), "median_B_per_s": med, "ratios_sigma_over_table": ratios,
            "cost_table_over_sigma": [lo, hi], "verdict": verdict,
            "switch_deadline_fraction_of_sigma_work": [deadline(hi), deadline(lo)],
            "rule_v2_rate_cost_vs_v1_session": rule_cost}
        print("%-7s %d paired reps: sigma %.3f, table %.3f B/s (medians); sigma/table %.3f - %.3f; "
              "cost table/sigma %.3f - %.3f -> %s, switch deadline %.1f%% - %.1f%% of sigma's work%s" % (
                  geometry, len(reps), med["sigma"], med["table"], min(ratios), max(ratios), lo, hi, verdict,
                  100 * deadline(hi), 100 * deadline(lo),
                  "" if old is None else "; v2 rate cost against the v1 session %+.1f%%" % (100 * rule_cost)))
    valid = all(verified.values()) and not failures and not missing and len(verdicts) == len(GEOMETRIES)
    overall = ("pays" if all(v == "pays" for v in verdicts) else
               "does not pay" if "does not pay" in verdicts else "undecided") if valid else "INVALID"
    if not valid:
        print("INVALID: verified %s, failures %s, samples with no rate %s, geometries %d of %d" % (
            verified, failures, summary["missing_samples"], len(verdicts), len(GEOMETRIES)))
    summary["verdict"] = overall
    print("verdict:", overall)
    if freeze:
        freeze.mkdir(parents=True, exist_ok=True)
        for f in res.iterdir():
            if f.suffix in (".txt", ".log", ".json") or f.name == "exit-code":
                shutil.copy(f, freeze / f.name)
        if (run / "launch.json").exists():
            shutil.copy(run / "launch.json", freeze / "launch.json")
        (freeze / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
