#!/usr/bin/env python3
"""WALK-CONSTANT.md section 11.6: the device harness's seed-to-seed scatter.

    python3 benchmarks/walk-constant/scatter.py

Reads scatter-v2.jsonl and scatter-v1.jsonl (scatter.sh writes them) and the
emulation's matching rows from matrix-v3.jsonl (rule v2) and matrix-v2.jsonl
(rule v1), and applies the tests section 11.6 declared before the runs:

1. per rule, the chi-square of the 16 constants about their inverse-variance
   mean on 15 degrees of freedom: scatter established at p < 0.01,
   suggestive at p < 0.05, absent otherwise;
2. the reading of the two verdicts together;
3. per rule, the pooled constant against the emulation's, with both errors
   (the pooled error scaled by sqrt(chi2 / dof) where test 1 established
   scatter): the harnesses disagree if |z| > 2.58.
"""
import json
import math
import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
DECLARED = {"v2": list(range(240, 256)), "v1": list(range(260, 276))}
SHAPE = {"n": 23, "walk": "device-table", "dist": "ecc2k130", "branches": 8, "walks": 8}
EMULATION = {"v2": "matrix-v3.jsonl", "v1": "matrix-v2.jsonl"}


def chi2_sf(x, dof):
    """P(chi2 with dof degrees of freedom > x), by the regularized incomplete gamma."""
    a, x = dof / 2, x / 2
    if x <= 0:
        return 1.0
    front = math.exp(-x + a * math.log(x) - math.lgamma(a))
    if x < a + 1:
        term = total = 1 / a
        n = a
        while term > 1e-15 * total:
            n += 1
            term *= x / n
            total += term
        return 1 - front * total
    b, c, d = x + 1 - a, 1e300, 1 / (x + 1 - a)
    h = d
    for i in range(1, 1000):
        an = -i * (i - a)
        b += 2
        d = 1 / (an * d + b)
        c = b + an / c
        h *= d * c
    return front * h


def rows(name):
    return [json.loads(line) for line in (HERE / name).read_text().splitlines() if line.strip()]


def emulation(rule):
    for r in rows(EMULATION[rule]):
        if (r.get("walk") == "table" and r.get("n") == 23 and r.get("dist") == "ecc2k130"
                and r.get("branches") == 8 and r.get("walks") == 8):
            return r["c"], r["c_se"]
    raise SystemExit("no emulation row for %s in %s" % (rule, EMULATION[rule]))


def main():
    result, problems = {}, []
    for rule in ("v2", "v1"):
        runs = rows("scatter-%s.jsonl" % rule)
        seeds = [r["seed"] for r in runs]
        if sorted(seeds) != DECLARED[rule]:
            problems.append("%s: seeds %s, declared %s" % (rule, seeds, DECLARED[rule]))
        for r in runs:
            bad = {k: r.get(k) for k, v in SHAPE.items() if r.get(k) != v}
            if bad:
                problems.append("%s seed %s: %s" % (rule, r["seed"], bad))
        w = [1 / r["c_se"] ** 2 for r in runs]
        mean = sum(wi * r["c"] for wi, r in zip(w, runs)) / sum(w)
        se = 1 / math.sqrt(sum(w))
        chi2 = sum(wi * (r["c"] - mean) ** 2 for wi, r in zip(w, runs))
        dof = len(runs) - 1
        p = chi2_sf(chi2, dof)
        verdict = "established" if p < 0.01 else "suggestive" if p < 0.05 else "absent"
        scaled = se * math.sqrt(chi2 / dof) if verdict == "established" else se
        emu, emu_se = emulation(rule)
        z = (mean - emu) / math.hypot(scaled, emu_se)
        result[rule] = verdict
        print("rule %s, seeds %d-%d:" % (rule, min(seeds), max(seeds)))
        print("  seed  trials       c      se   z about the mean")
        for r in runs:
            print("  %4d  %6d  %.5f  %.5f  %+5.2f" % (r["seed"], r["trials"], r["c"], r["c_se"],
                                                  (r["c"] - mean) / r["c_se"]))
        print("  pooled c = %.5f +- %.5f over %d trials" % (mean, se, sum(r["trials"] for r in runs)))
        print("  test 1: chi2 = %.2f on %d dof, p = %.4f -> scatter %s" % (chi2, dof, p, verdict))
        print("  test 3: against the emulation's %.5f +- %.5f (%s): z = %+.2f -> %s" % (
            emu, emu_se, EMULATION[rule], z, "the harnesses disagree" if abs(z) > 2.58 else "agree"))
    v2, v1 = result["v2"] == "established", result["v1"] == "established"
    reading = ("the scatter is new with rule v2" if v2 and not v1 else
               "the harness's standard error understates its scatter under either rule" if v2 and v1 else
               "rule v1 scatters and rule v2 does not; reported as found" if v1 else
               "no scatter under either rule: seeds 230 and 231's gap and the 13.3 were chance "
               "and the n = 41 reference outlier")
    print("test 2:", reading)
    if problems:
        print("NOT AS DECLARED:\n  " + "\n  ".join(problems))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
