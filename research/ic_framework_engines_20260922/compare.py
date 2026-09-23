#!/usr/bin/env python3
"""Compare the engines of one suite run, against the targets declared
in the ledger's §17.2 before the run.

    python3 research/ic_framework_engines_20260922/compare.py --run baseline_v1
    python3 research/ic_framework_engines_20260922/compare.py --run baseline_v1 --freeze-manifest
    python3 research/ic_framework_engines_20260922/compare.py --run candidate_x --manifest manifest.json

Writes ``results/<run>/compare.json`` and ``compare.md``.  It reads
reports and never computes a number the reports do not support: every
ratio is between two measured wall times of the same target in the same
process (stage) or two verified whole runs on the same planted target
(bench), and every interval is a percentile bootstrap over those pairs.

Checks it fails on, rather than reports around:

- an engine whose decided answer differs from the reference's anywhere;
- with ``--manifest``: a system fingerprint or a cell's verdict digest
  that differs from the frozen run's — a candidate that solved other
  inputs, or decided the same inputs differently, is not comparable;
- a whole-pipeline row whose logarithm did not verify.
"""

import argparse
import json
import math
import random
import statistics
import sys
import zlib
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent


def bootstrap_median(values, resamples, seed, interval):
    """Percentile bootstrap of the median."""
    if not values:
        return None
    rng = random.Random(seed)
    k = len(values)
    meds = sorted(statistics.median(rng.choices(values, k=k)) for _ in range(resamples))
    lo = meds[int((1 - interval) / 2 * resamples)]
    hi = meds[min(resamples - 1, int((1 + interval) / 2 * resamples))]
    return {"median": statistics.median(values), "lo": lo, "hi": hi, "pairs": k}


def decided(run):
    return run["verdict"] in ("solved", "unsatisfiable")


def target_wall(runs):
    """Median wall (ns) over repetitions if every repetition decided,
    else None; and whether the engine hit the budget on this target."""
    budget = any(r["verdict"] == "budget" for r in runs)
    if not runs or budget or not all(decided(r) for r in runs):
        return None, budget
    return statistics.median(r["wall_ns"] for r in runs), budget


def stage_rows(report, suite, boot):
    """Per cell: per-engine paired ratios to the baseline and to the reference."""
    out = []
    for cell in report["cells"]:
        ref = cell["reference_engine"]
        base = suite["baseline_engine"]
        engines = [e["engine"] for e in cell["engines"]]
        per = {e: [target_wall(t["runs"].get(e, [])) for t in cell["per_target"]] for e in engines}
        summaries = {e["engine"]: e for e in cell["engines"]}
        row = {
            "family": cell["family"], "n": cell["n"], "n_prime": cell["n_prime"], "m": cell["summands"],
            "vars": cell["n_vars"], "reference": ref, "verdict_digest": cell["verdict_digest"],
            "d_sr": [cell["d_semireg_min"], cell["d_semireg_max"]], "engines": {},
        }
        for e in engines:
            s = summaries[e]
            walls = per[e]
            vs_ref = [w / r for (w, _), (r, _) in zip(walls, per.get(ref, [])) if w and r]
            vs_base = [w / b for (w, _), (b, _) in zip(walls, per.get(base, [])) if w and b] if base in per else []
            # Targets the baseline could not decide but this engine did:
            # the ratio there is at most this engine's wall over the budget.
            beyond = sum(1 for (w, _), (b, bb) in zip(walls, per.get(base, [])) if w and bb) if base in per else 0
            # A stable per-(cell, engine) stream: Python's own `hash` of a
            # string is salted per process and would not reproduce.
            seed = zlib.crc32(f"{cell['family']}/{cell['n']}/{cell['summands']}/{e}".encode())
            row["engines"][e] = {
                "workload": s["workload"],
                "declined": s["declined"],
                "decided": s["decided"],
                "targets": cell["targets"],
                "over_budget": s["over_budget"],
                "agrees": s["agrees_with_reference"],
                "d_learn_mean": s["d_learn_mean"], "d_learn_max": s["d_learn_max"],
                "d_reach_mean": s["d_reach_mean"], "d_reach_max": s["d_reach_max"],
                "ms_mean": s["ms_mean"], "op_unit": s["op_unit"], "ops_mean": s["ops_mean"],
                "vs_reference": bootstrap_median(vs_ref, boot["resamples"], boot["seed"] ^ seed, boot["interval"]) if e != ref else None,
                "vs_baseline": bootstrap_median(vs_base, boot["resamples"], boot["seed"] ^ seed ^ 1, boot["interval"]) if e != base else None,
                "baseline_budget_but_decided_here": beyond,
                "_vs_ref_pairs": vs_ref,
            }
        out.append(row)
    return out


def bench_pairs(reports, engines):
    """(n, seed, recovered) -> engine -> row, verified rows only."""
    table = defaultdict(dict)
    unverified = []
    for part, rep in reports:
        for row in rep["rows"]:
            spec = row["spec"]
            name = spec.get("solver") if spec["oracle"] == "descent-algebraic" else "pair-table"
            if not row.get("verified"):
                unverified.append((part, row["label"]))
                continue
            key = (row["instance"], spec["seed"], row["recovered"])
            table[key][name] = row
    return table, unverified


def per_relation(row):
    rel = row["decomposition"]["relations_found"]
    return row["total_gae"] / rel if rel else None


def fmt_ci(ci, digits=2):
    if not ci:
        return "—"
    f = lambda v: f"{v:,.0f}" if v >= 100 else f"{v:.{digits + 1}g}"
    return f"{f(ci['median'])} [{f(ci['lo'])}, {f(ci['hi'])}] (k={ci['pairs']})"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", required=True)
    ap.add_argument("--suite", default="suite.json",
                    help="the frozen suite file in this directory (default suite.json)")
    ap.add_argument("--manifest", help="check input fingerprints and verdict digests against this frozen manifest")
    ap.add_argument("--freeze-manifest", action="store_true", help="write manifest.json from this run")
    args = ap.parse_args()

    suite = json.loads((HERE / args.suite).read_text())
    boot = suite["analysis"]["bootstrap"]
    run_dir = HERE / "results" / args.run
    problems = []

    stage, bench = [], []
    for kind, key in (("stage", "stage_parts"), ("bench", "bench_parts")):
        for part in suite.get(key, []):
            path = run_dir / f"{part['id']}.json"
            if not path.exists():
                problems.append(f"missing part {part['id']}")
                continue
            rep = json.loads(path.read_text())
            if rep.get("status") != "complete":
                problems.append(f"part {part['id']} status {rep.get('status')}")
            (stage if kind == "stage" else bench).append((part["id"], rep))

    # ── Input identity ──────────────────────────────────────────────
    manifest = {"suite_version": suite["version"], "run": args.run, "cells": {}}
    for pid, rep in stage:
        for cell in rep["cells"]:
            key = f"{pid}/{cell['family']}/{cell['n']}:{cell['n_prime']}:{cell['summands']}"
            manifest["cells"][key] = {
                "verdict_digest": cell["verdict_digest"],
                "targets": [{"x_r": t["x_r"], "system_blake3": t["system_blake3"]} for t in cell["per_target"]],
            }
    if args.manifest:
        frozen = json.loads(Path(args.manifest).read_text())
        for key, cell in frozen["cells"].items():
            mine = manifest["cells"].get(key)
            if mine is None:
                problems.append(f"manifest cell {key} not in this run")
            elif mine["targets"] != cell["targets"]:
                problems.append(f"cell {key}: input fingerprints differ from the manifest")
            elif mine["verdict_digest"] != cell["verdict_digest"]:
                problems.append(f"cell {key}: verdict digest differs from the manifest")
    if args.freeze_manifest:
        path = HERE / "manifest.json"
        if path.exists():
            sys.exit("manifest.json exists; it is frozen")
        path.write_text(json.dumps(manifest, indent=1) + "\n")

    # ── Stage ───────────────────────────────────────────────────────
    cells = []
    for pid, rep in stage:
        for row in stage_rows(rep, suite, boot):
            row["part"] = pid
            cells.append(row)
            for e, r in row["engines"].items():
                if r["decided"] and not r["agrees"]:
                    problems.append(f"{pid} {row['family']} {row['n']}:{row['n_prime']}:{row['m']}: {e} disagrees with {row['reference']}")

    two = [c for c in cells if c["m"] == 2]
    f4 = suite["f4_family"]
    alg = suite["algebraic_engines"]

    # Target 1: reach.
    reach = {}
    for e in f4:
        span = [c for c in two if 16 <= c["vars"] <= 22 and e in c["engines"]]
        reach[e] = {
            "cells": len(span),
            "all_decided": bool(span) and all(c["engines"][e]["decided"] == c["engines"][e]["targets"] for c in span),
            "all_agree": all(c["engines"][e]["agrees"] for c in span),
        }
    t1 = any(v["all_decided"] and v["all_agree"] and v["cells"] for v in reach.values()) \
        if any(v["cells"] for v in reach.values()) else None

    # Target 3: a stage crossover against fes-f2.
    crossings = []
    for c in two:
        if c["reference"] not in ("fes-f2", "fes-f2-wide"):
            continue
        for e in alg:
            ci = c["engines"].get(e, {}).get("vs_reference")
            if ci and ci["median"] < 1 and ci["hi"] < 1:
                crossings.append({"part": c["part"], "family": c["family"], "n": c["n"], "vars": c["vars"], "engine": e, "ci": ci})
    t3 = bool(crossings)

    # The trajectory: pooled over families and parts at each size.
    traj = {}
    for e in alg:
        by_vars = defaultdict(list)
        for c in two:
            if c["reference"] in ("fes-f2", "fes-f2-wide") and e in c["engines"]:
                by_vars[c["vars"]] += c["engines"][e]["_vs_ref_pairs"]
        traj[e] = {v: bootstrap_median(p, boot["resamples"], boot["seed"] ^ v, boot["interval"]) for v, p in sorted(by_vars.items()) if p}
    gaining = {e: (t[26]["median"] < t[16]["median"]) if (t.get(26) and t.get(16)) else None for e, t in traj.items()}
    abandon = (not any(gaining.values())) if any(g is not None for g in gaining.values()) else None

    # Extrapolation, marked as such: log-linear fit of the ratio in the
    # number of unknowns over the sizes each engine decided, 16 and up.
    extrapolation = {}
    for e, t in traj.items():
        pts = [(v, math.log(ci["median"])) for v, ci in t.items() if v >= 16 and ci and ci["median"] > 0]
        if len(pts) >= 4:
            xs, ys = zip(*pts)
            mx, my = statistics.fmean(xs), statistics.fmean(ys)
            sxx = sum((x - mx) ** 2 for x in xs)
            slope = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sxx
            icpt = my - slope * mx
            cross = (-icpt / slope) if slope < 0 else None
            extrapolation[e] = {
                "sizes": list(xs), "log_ratio_per_unknown": slope,
                "ratio_factor_per_two_unknowns": math.exp(2 * slope),
                "crossover_unknowns_extrapolated": cross,
                "marked": "extrapolation from a log-linear fit over the sizes listed, not a measurement",
            }

    # ── Whole method ────────────────────────────────────────────────
    table, unverified = bench_pairs(bench, None)
    for part, label in unverified:
        problems.append(f"{part}: row `{label}` did not verify")
    gate = {}
    by_n = defaultdict(lambda: defaultdict(list))
    s_rows = defaultdict(lambda: defaultdict(list))
    for (inst, seed, _), rows in table.items():
        n = int(inst.split("-n")[1].split("-")[0]) if "-n" in inst else None
        for name, row in rows.items():
            s_rows[n][name].append(row)
        base = rows.get(suite["baseline_engine"])
        if not base:
            continue
        for e in f4 + ["crossbred-f2", "fes-f2", "exhaustive"]:
            if e in rows:
                a, b = per_relation(rows[e]), per_relation(base)
                if a and b:
                    by_n[n][e].append(a / b)
    for n, engines in sorted(by_n.items()):
        gate[n] = {e: bootstrap_median(v, boot["resamples"], boot["seed"] ^ n, boot["interval"]) for e, v in engines.items()}
    t2_by_n = {n: any(ci and ci["median"] <= 0.8 and ci["hi"] < 1 for e, ci in g.items() if e in f4) for n, g in gate.items()}
    t2 = all(t2_by_n.get(n, False) for n in (13, 15)) if t2_by_n else None

    whole = {}
    for n, engines in sorted(s_rows.items()):
        whole[n] = {}
        for name, rows in engines.items():
            ss = [r["s"] for r in rows]
            rho = [r["s_over_rho"] for r in rows if r.get("s_over_rho") is not None]
            whole[n][name] = {
                "rows": len(rows),
                "s_median": statistics.median(ss),
                "s_min": min(ss), "s_max": max(ss),
                "s_over_rho_median": statistics.median(rho) if rho else None,
                "per_relation_gae_median": statistics.median(pr) if (pr := [per_relation(r) for r in rows if per_relation(r)]) else None,
            }
        pt = whole[n].get("pair-table")
        if pt:
            for name, w in whole[n].items():
                w["s_over_pair_table"] = w["s_median"] / pt["s_median"]

    for c in cells:
        for r in c["engines"].values():
            r.pop("_vs_ref_pairs", None)

    result = {
        "run": args.run,
        "suite_version": suite["version"],
        "problems": problems,
        "targets": {
            "1_reach": {"met": t1, "by_engine": reach},
            "2_gate": {"met": t2, "by_n": {str(k): v for k, v in t2_by_n.items()}, "ratios": {str(k): v for k, v in gate.items()}},
            "3_question": {"met": t3, "crossings": crossings},
            "abandon": {"triggered": abandon, "gaining_16_to_26": gaining},
        },
        "trajectory_vs_fes": {e: {str(v): ci for v, ci in t.items()} for e, t in traj.items()},
        "extrapolation": extrapolation,
        "whole_method": {str(k): v for k, v in whole.items()},
        "cells": cells,
    }
    (run_dir / "compare.json").write_text(json.dumps(result, indent=1) + "\n")

    md = [f"# Engine suite: run `{args.run}`\n"]
    md.append("Problems: " + ("none" if not problems else "; ".join(problems)) + "\n")
    md.append("## Declared targets (ledger §17.2)\n")
    verdict = lambda v, yes="met", no="not met": "no data in this run" if v is None else (yes if v else no)
    md.append(f"- **1 reach**: {verdict(t1)} — " + ", ".join(f"{e}: {v['cells']} cells, all decided {v['all_decided']}, agree {v['all_agree']}" for e, v in reach.items()))
    md.append(f"- **2 gate**: {verdict(t2)} — " + "; ".join(f"n={n}: {v}" for n, v in t2_by_n.items()))
    md.append(f"- **3 question**: {'met' if t3 else 'not met'} — {len(crossings)} crossing cell(s)" + "".join(
        f"; {c['engine']} at {c['vars']} unknowns ({c['family']}, {c['part']}): {fmt_ci(c['ci'])}" for c in crossings))
    md.append(f"- **abandon**: {verdict(abandon, 'triggered', 'not triggered')} — gaining 16→26: {gaining}\n")
    md.append("## Stage: wall / fes-f2, pooled over families, median [95% CI] (pairs)\n")
    sizes = sorted({v for t in traj.values() for v in t})
    md.append("| engine | " + " | ".join(f"{v}" for v in sizes) + " |")
    md.append("|:--|" + "--:|" * len(sizes))
    for e, t in traj.items():
        md.append(f"| {e} | " + " | ".join(fmt_ci(t.get(v)) for v in sizes) + " |")
    md.append("\n## Extrapolation (marked: not a measurement)\n")
    for e, x in extrapolation.items():
        md.append(f"- {e}: ×{x['ratio_factor_per_two_unknowns']:.2f} per two unknowns over {x['sizes']}; crossover at {x['crossover_unknowns_extrapolated'] and round(x['crossover_unknowns_extrapolated'], 1)} unknowns")
    md.append("\n## Stage cells\n")
    md.append("| part | E | n:n':m | vars | engine | finds | decided | budget | D_learn | D_reach | D_sr | ms | vs reference | vs baseline | agrees |")
    md.append("|:--|:--|:--|--:|:--|:--|--:|--:|--:|--:|:--|--:|--:|--:|:--|")
    for c in cells:
        for e, r in c["engines"].items():
            if r["declined"]:
                continue
            dsr = f"{c['d_sr'][0]}" if c["d_sr"][0] == c["d_sr"][1] else f"{c['d_sr'][0]}–{c['d_sr'][1]}"
            dl = "—" if r["d_learn_mean"] is None else f"{r['d_learn_mean']:.2f}"
            dr = "—" if r["d_reach_mean"] is None else f"{r['d_reach_mean']:.2f}"
            vr = "ref" if e == c["reference"] else fmt_ci(r["vs_reference"])
            vb = "base" if e == suite["baseline_engine"] else fmt_ci(r["vs_baseline"])
            if r["baseline_budget_but_decided_here"]:
                vb += f" (+{r['baseline_budget_but_decided_here']} where the baseline hit the budget)"
            ag = "—" if not r["decided"] else ("yes" if r["agrees"] else "NO")
            md.append(f"| {c['part']} | {c['family']} | {c['n']}:{c['n_prime']}:{c['m']} | {c['vars']} | {e} | {'first' if r['workload'] == 'first solution' else 'all'} | {r['decided']}/{r['targets']} | {r['over_budget']} | {dl} | {dr} | {dsr} | {r['ms_mean']:.2f} | {vr} | {vb} | {ag} |")
    md.append("\n## Whole method (S; per-relation GAE; paired over curve × planted target)\n")
    for n, w in whole.items():
        md.append(f"### n = {n}\n")
        md.append("| row | runs | S median | S range | S / pair table | S / rho | GAE per relation |")
        md.append("|:--|--:|--:|:--|--:|--:|--:|")
        for name, x in sorted(w.items(), key=lambda kv: kv[1]["s_median"]):
            rho = "—" if x["s_over_rho_median"] is None else f"{x['s_over_rho_median']:,.2f}"
            md.append(f"| {name} | {x['rows']} | {x['s_median']:,.1f} | {x['s_min']:,.1f}–{x['s_max']:,.1f} | {x.get('s_over_pair_table', float('nan')):,.2f} | {rho} | {'—' if x['per_relation_gae_median'] is None else format(x['per_relation_gae_median'], ',.0f')} |")
        if n in gate:
            md.append("\nPer verified relation over buchberger-f2, median [95% CI]: " + "; ".join(f"{e} {fmt_ci(ci, 4)}" for e, ci in gate[n].items()) + "\n")
    (run_dir / "compare.md").write_text("\n".join(md) + "\n")
    print("\n".join(md[:8]))
    if problems:
        sys.exit(1)


if __name__ == "__main__":
    main()
