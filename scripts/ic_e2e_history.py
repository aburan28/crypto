#!/usr/bin/env python3
"""Run history for the end-to-end index-calculus benchmark: paired reruns on
matched hardware, and the interval AGENTS.md section 8 asks a runtime claim to
carry.

``add``     extracts one compact record from a CI run artifact (the directory
            ``ic_e2e_benchmark.py run`` wrote, as uploaded by the workflow) into
            ``docs/ic/ci/runs/<run-id>.json``.  Never overwrites.
``report``  reads every record, checks that the pinned counters are identical
            across runs and equal to the frozen reference (otherwise the runs
            are not the same algorithm and must not be pooled), and reports per
            rung: the per-run whole-process ratios rho / IC, their mean, sample
            standard deviation, a t-based 95 % interval, and the same for the
            paired difference rho - IC in seconds.  A rung's runtime claim
            "faster than rho end to end" holds only when the ratio interval's
            lower bound exceeds 1 (equivalently the difference interval excludes
            0) with at least three runs.

Wall time stays a practicality note: this satisfies the *runtime* half of
section 8 on matched hardware and says nothing about S, which remains null
for the index-calculus side (no measured unit conversion).
"""
from __future__ import annotations

import argparse
import json
import math
import statistics
import sys
from pathlib import Path

SCHEMA_VERSION = 1

# Two-sided 97.5 % Student t quantiles, df = 1..30 (then normal).
T975 = [None, 12.706, 4.303, 3.182, 2.776, 2.571, 2.447, 2.365, 2.306, 2.262, 2.228,
        2.201, 2.179, 2.160, 2.145, 2.131, 2.120, 2.110, 2.101, 2.093, 2.086,
        2.080, 2.074, 2.069, 2.064, 2.060, 2.056, 2.052, 2.048, 2.045, 2.042]


def t975(df: int) -> float:
    if df < 1:
        raise ValueError("need at least two samples")
    return T975[df] if df <= 30 else 1.960


def interval(values: list[float]) -> dict:
    n = len(values)
    mean = statistics.fmean(values)
    if n < 2:
        return {"n": n, "mean": mean, "sd": None, "ci95": None}
    sd = statistics.stdev(values)
    half = t975(n - 1) * sd / math.sqrt(n)
    return {"n": n, "mean": mean, "sd": sd, "ci95": [mean - half, mean + half]}


def extract(run_dir: Path, run_id: str, url: str | None) -> dict:
    manifest = json.loads((run_dir / "manifest.json").read_text())
    summary = json.loads((run_dir / "summary.json").read_text())
    if not summary.get("ok"):
        raise SystemExit(f"{run_dir}: the gate did not pass in this run; a failed run is not a rerun")
    rungs = {}
    for r in summary["rungs"]:
        report = json.loads((run_dir / r["params"].split("/")[-1].replace(".json", "") / "report.json").read_text()) \
            if (run_dir / r["params"].split("/")[-1].replace(".json", "") / "report.json").exists() else None
        counters = None
        if report is not None:
            sys.path.insert(0, str(Path(__file__).resolve().parent))
            import ic_e2e_benchmark as bench  # noqa: E402
            counters = bench.counters_of(report)
        w = r["wall"]["now"]
        rungs[r["params"]] = {
            "degree": r["degree"], "subgroup_bits": r["subgroup_bits"], "targets": r["targets"],
            "ic_verified": r["ic_verified"], "rho_verified": r["rho_verified"],
            "counters": counters,
            "ic_whole_process_seconds": w["ic_whole_process_seconds"], "rho_seconds_total": w["rho_seconds_total"],
            "whole_process_ratio": w["whole_process_ratio"], "charged_ratio": w["charged_ratio"],
            "whole_process_crossover": w["whole_process_crossover"],
        }
    return {
        "schema_version": SCHEMA_VERSION, "run_id": run_id, "url": url,
        "started_at": manifest.get("started_at"), "host": manifest.get("host"),
        "ic_binary_sha256": manifest.get("ic_binary", {}).get("sha256"), "git": manifest.get("git"),
        "reference": summary.get("reference", {}).get("path"),
        "rungs": rungs,
    }


def cmd_add(args: argparse.Namespace) -> int:
    out = Path(args.out_dir) / f"{args.run_id}.json"
    if out.exists():
        print(f"refusing to overwrite {out}", file=sys.stderr)
        return 2
    rec = extract(Path(args.run_dir), args.run_id, args.url)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(rec, indent=1, sort_keys=True) + "\n")
    print(f"recorded run {args.run_id}: {len(rec['rungs'])} rung(s)")
    return 0


def build_report(records: list[dict], reference: dict | None) -> dict:
    hosts = sorted({json.dumps(r["host"], sort_keys=True) for r in records})
    report = {"schema_version": SCHEMA_VERSION, "runs": [{k: r[k] for k in ("run_id", "url", "started_at", "ic_binary_sha256")} for r in records],
              "distinct_hosts": len(hosts), "host": json.loads(hosts[0]) if len(hosts) == 1 else None, "rungs": {}, "problems": []}
    if len(hosts) > 1:
        report["problems"].append("runs come from more than one host; matched-hardware intervals are not defined across them")
    params = sorted({p for r in records for p in r["rungs"]})
    for p in params:
        rows = [(r["run_id"], r["rungs"][p]) for r in records if p in r["rungs"]]
        counters = [json.dumps(x["counters"], sort_keys=True) for _, x in rows if x["counters"] is not None]
        identical = len(set(counters)) <= 1
        ref_ok = True
        if reference is not None and p in reference["rungs"] and counters:
            ref_ok = json.loads(counters[0]) == reference["rungs"][p]["counters"]
        if not identical:
            report["problems"].append(f"{p}: counters differ between runs; not one algorithm")
        if not ref_ok:
            report["problems"].append(f"{p}: counters differ from the frozen reference")
        ratios = [x["whole_process_ratio"] for _, x in rows]
        diffs = [x["rho_seconds_total"] - x["ic_whole_process_seconds"] for _, x in rows]
        ri, di = interval(ratios), interval(diffs)
        claim = bool(ri["ci95"] and len(ratios) >= 3 and ri["ci95"][0] > 1.0 and di["ci95"] and di["ci95"][0] > 0.0)
        report["rungs"][p] = {
            "degree": rows[0][1]["degree"], "subgroup_bits": rows[0][1]["subgroup_bits"], "targets": rows[0][1]["targets"],
            "runs": [{"run_id": rid, "ratio": x["whole_process_ratio"], "ic_seconds": x["ic_whole_process_seconds"],
                      "rho_seconds": x["rho_seconds_total"], "charged_ratio": x["charged_ratio"]} for rid, x in rows],
            "all_verified": all(x["ic_verified"] == x["targets"] == x["rho_verified"] for _, x in rows),
            "counters_identical_across_runs": identical, "counters_equal_reference": ref_ok,
            "ratio": ri, "rho_minus_ic_seconds": di,
            "runtime_claim_faster_than_rho_end_to_end": claim if (identical and ref_ok) else False,
            "claim_rule": "n >= 3 matched runs, 95% t-interval of the per-run rho/IC ratio above 1 and of rho - IC seconds above 0, counters identical",
        }
    return report


def render_markdown(report: dict) -> str:
    lines = ["## e2e benchmark run history (matched hardware)", ""]
    h = report["host"]
    lines.append(f"{len(report['runs'])} run(s) on {report['distinct_hosts']} host" + ("" if report["distinct_hosts"] == 1 else "s") +
                 (f" ({h['platform']}, {h['cpu_count']} vCPU)" if h else "") + ". Wall time is a practicality note; `S` stays null.")
    lines.append("")
    lines.append("| rung | log₂ r | runs | rho/IC per run | mean ± sd | 95% CI | rho − IC seconds, 95% CI | counters | verified | runtime claim |")
    lines.append("|:--|--:|--:|:--|:--|:--|:--|:--|:--|:--|")
    for p, r in report["rungs"].items():
        ri, di = r["ratio"], r["rho_minus_ic_seconds"]
        lines.append("| `%s` | %.1f | %d | %s | %.3f ± %s | %s | %s | %s | %s | %s |" % (
            p.split("/")[-1], r["subgroup_bits"], ri["n"], ", ".join("%.3f" % x["ratio"] for x in r["runs"]),
            ri["mean"], ("%.3f" % ri["sd"]) if ri["sd"] is not None else "—",
            ("[%.3f, %.3f]" % tuple(ri["ci95"])) if ri["ci95"] else "n < 2",
            ("[%.1f, %.1f]" % tuple(di["ci95"])) if di["ci95"] else "n < 2",
            "identical" if r["counters_identical_across_runs"] and r["counters_equal_reference"] else "DIFFER",
            "all" if r["all_verified"] else "NOT ALL",
            "**faster than rho, end to end**" if r["runtime_claim_faster_than_rho_end_to_end"] else "not established"))
    lines.append("")
    lines.append("Runs: " + ", ".join(f"[{x['run_id']}]({x['url']})" if x.get("url") else x["run_id"] for x in report["runs"]) + ".")
    for pr in report["problems"]:
        lines.append(f"- **problem**: {pr}")
    lines.append("")
    return "\n".join(lines)


def cmd_report(args: argparse.Namespace) -> int:
    records = sorted((json.loads(Path(p).read_text()) for p in args.runs), key=lambda r: r["started_at"] or "")
    reference = json.loads(Path(args.reference).read_text()) if args.reference else None
    report = build_report(records, reference)
    text = render_markdown(report)
    if args.summary_json:
        Path(args.summary_json).write_text(json.dumps(report, indent=1, sort_keys=True) + "\n")
    if args.summary_markdown:
        Path(args.summary_markdown).write_text(text)
    print(text)
    return 1 if report["problems"] else 0


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="command", required=True)
    a = sub.add_parser("add", help="record one CI run artifact")
    a.add_argument("--run-dir", required=True, help="the artifact directory (manifest.json, summary.json, <rung>/report.json)")
    a.add_argument("--run-id", required=True)
    a.add_argument("--url")
    a.add_argument("--out-dir", default="docs/ic/ci/runs")
    a.set_defaults(func=cmd_add)
    r = sub.add_parser("report", help="intervals over recorded runs")
    r.add_argument("--runs", nargs="+", required=True)
    r.add_argument("--reference", help="frozen reference whose counters every run must match")
    r.add_argument("--summary-json")
    r.add_argument("--summary-markdown")
    r.set_defaults(func=cmd_report)
    return ap


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    sys.exit(main())
