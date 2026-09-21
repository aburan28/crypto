#!/usr/bin/env python3
"""End-to-end index-calculus benchmark: run the ledger rungs, gate, freeze.

The unit of this check is one ``ic workflow`` run per parameter file: select
→ collect → logs → solve, followed by the signed-Frobenius ρ baseline on the
same known-answer targets in the same process.  That is the whole method, cold,
from setup to recovered logarithm, so the quantities gated here are end-to-end
quantities.  A stage number never passes this gate on its own.

Three subcommands:

``run``     builds nothing; it invokes an existing ``ic`` binary once per
            parameter file, keeps every report and stderr under the output
            directory and writes a manifest (git commit, binary hash, host).
``check``   compares the reports against a frozen reference and fails closed.
``freeze``  writes a new reference from a run directory.  It refuses to
            overwrite: a superseded reference stays in the tree as the
            "before" mark, per AGENTS.md §7.

What ``check`` gates, per rung:

1. **Correctness.**  Status ``complete``, every target recovered and verified
   as ``[d]G = Q`` (and equal to the planted scalar where one was planted),
   and every ρ walk verified on the same targets.  A row without a verified
   answer is not a result.
2. **Pinned counters.**  Collection trials, summands scanned, relations,
   pair-table pairs, factor-base points and columns, descent trials, ρ
   iterations and ρ group additions must equal the reference exactly.  The
   probe sequence, the descent and the ρ walk are all seeded, so a drift means
   the algorithm changed; that is allowed, but only through a deliberate
   ``freeze`` that records the new numbers next to the old ones.
3. **End-to-end wall ratio.**  ``rho_seconds_total / (select + collect + logs
   + pair table + descent)`` — ρ against the whole method, measured on the same
   host in the same process — must not fall below the reference by more than
   the tolerance, and a rung whose reference already crosses ρ end to end must
   still cross.  Wall time is a practicality note under AGENTS.md §6, which is
   why it is gated only as a paired same-host ratio with a wide tolerance and
   never as an absolute.
4. **Baseline sanity.**  ρ's measured group additions per √r per target must
   sit within a factor of the expected ``√(π/(4n))`` for the signed-Frobenius
   walk; a broken baseline makes every ratio meaningless.

What it does not do: it does not compute an ``S`` for the index-calculus side.
The counters it pins live in different units (group additions for the pair
table and ρ, table lookups for the scan, scalar multiplications for the probes)
and this repository records no measured conversion between them, so the IC
``S`` is reported as null with that reason (AGENTS.md §8).  Passing this check
is at most engineering in the §3 sense; it is never an advance.
"""
from __future__ import annotations

import argparse
import datetime as _dt
import hashlib
import json
import math
import os
import platform
import subprocess
import sys
from pathlib import Path
from typing import Any

SCHEMA_VERSION = 1
MANIFEST_FILE = "manifest.json"
REPORT_FILE = "report.json"
STDERR_FILE = "stderr.txt"
RUN_DIR = "run"

# Report fields pinned exactly.  Each entry is (label, path into the report).
# ``stages.<name>`` addresses the stage report with that ``stage`` value.
COUNTER_PATHS: tuple[tuple[str, tuple[str, ...]], ...] = (
    ("factor_base_points", ("factor_base", "points")),
    ("factor_base_columns", ("factor_base", "columns")),
    ("pair_table_stored_pairs", ("factor_base", "pair_table_stored_pairs")),
    ("collection_trials", ("stages.collect", "trials_total")),
    ("collection_summands_scanned", ("stages.collect", "summands_scanned_total")),
    ("collection_relations", ("stages.collect", "relations_total")),
    ("descent_trials", ("stages.baseline", "vs_rho", "ic", "descent_trials_total")),
    ("rho_iterations", ("stages.baseline", "vs_rho", "rho", "iterations_total")),
)


class CheckFailure(Exception):
    """A gate failed; the message names the rung and the reason."""


# ── report access ──────────────────────────────────────────────────────────


def _stage(report: dict[str, Any], name: str) -> dict[str, Any] | None:
    for stage in report.get("stages", []):
        if stage.get("stage") == name:
            return stage
    return None


def _lookup(report: dict[str, Any], path: tuple[str, ...]) -> Any:
    node: Any = report
    for key in path:
        if key.startswith("stages."):
            node = _stage(report, key[len("stages."):])
        elif isinstance(node, dict):
            node = node.get(key)
        else:
            node = None
        if node is None:
            return None
    return node


def tier_of(report: dict[str, Any]) -> str:
    """Which representation of the pair table the run actually built.

    Pinned beside the counters because none of them can see it.
    ``pair_table_stored_pairs`` is ``|F|(|F|+1)/2`` for the full tier and
    for the compact one alike — they hold the same pairs and differ only
    in how a pair is stored — so a run that switched between them reads
    as "counters identical" while the algorithm has changed.  That is the
    one thing this gate exists to refuse, and it went through it once.
    """
    tier = _lookup(report, ("factor_base", "pair_table_tier"))
    if not isinstance(tier, str) or not tier:
        raise CheckFailure("factor_base.pair_table_tier missing: the report predates tier reporting")
    return tier


def counters_of(report: dict[str, Any]) -> dict[str, int]:
    out: dict[str, int] = {}
    for label, path in COUNTER_PATHS:
        value = _lookup(report, path)
        if not isinstance(value, int) or isinstance(value, bool):
            raise CheckFailure(f"counter {label!r} missing or not an integer in report")
        out[label] = value
    rows = _lookup(report, ("stages.baseline", "vs_rho", "targets_detail")) or []
    try:
        out["rho_group_additions"] = sum(int(r["walk_group_additions"]) for r in rows)
    except (KeyError, TypeError, ValueError) as exc:
        raise CheckFailure(f"rho targets_detail lacks walk_group_additions: {exc}") from exc
    return out


def wall_of(report: dict[str, Any]) -> dict[str, float | bool]:
    vs = _lookup(report, ("stages.baseline", "vs_rho"))
    if not isinstance(vs, dict):
        raise CheckFailure("baseline stage (vs_rho) missing: run the rung with baseline.rho = true")
    ic = vs["ic"]
    ic_whole = float(ic["precompute_seconds"]) + float(ic["pair_table_seconds"]) + float(ic["descent_seconds_total"])
    rho_total = float(vs["rho"]["seconds_total"])
    # Derived from the timings gated here, not read off the report's verdict:
    # correctness is enforced separately, so the crossover is the timing alone.
    crossover = ic_whole < rho_total
    return {
        "ic_whole_process_seconds": ic_whole,
        "ic_precompute_seconds": float(ic["precompute_seconds"]) + float(ic["pair_table_seconds"]),
        "ic_descent_seconds_total": float(ic["descent_seconds_total"]),
        "rho_seconds_total": rho_total,
        "whole_process_ratio": (rho_total / ic_whole) if ic_whole > 0 else math.inf,
        "charged_ratio": float(vs["ratio"]["charged"]),
        "amortised_ratio": float(vs["ratio"]["amortised"]),
        "whole_process_crossover": crossover,
        "reported_whole_process_crossover": bool(vs["verdict"]["whole_process_crossover"]),
    }


def rho_s_of(report: dict[str, Any], counters: dict[str, int]) -> dict[str, float]:
    """ρ's cost in the repository unit: group additions per target per √r.

    The signed-Frobenius walk on ``A = 2n`` classes is expected at
    ``√(πr/2)/√(2n)`` steps, i.e. ``S = √(π/(4n))``.
    """
    vs = _lookup(report, ("stages.baseline", "vs_rho"))
    r = int(vs["subgroup_order"])
    n = int(vs["n"])
    targets = int(vs["targets"])
    measured = counters["rho_group_additions"] / (targets * math.sqrt(r))
    expected = math.sqrt(math.pi / (4.0 * n))
    return {"measured": measured, "expected": expected, "measured_over_expected": measured / expected}


def verify_correctness(report: dict[str, Any]) -> tuple[int, int, int]:
    """Return (targets, ic_verified, rho_verified) or raise."""
    if report.get("operation") != "workflow":
        raise CheckFailure("report is not an ic workflow report")
    if report.get("status") != "complete":
        raise CheckFailure(f"workflow status {report.get('status')!r}, failure={report.get('failure')!r}")
    if report.get("failure"):
        raise CheckFailure(f"workflow recorded a failure: {report['failure']!r}")
    sol = report.get("solutions") or {}
    items = sol.get("items") or []
    vs = _lookup(report, ("stages.baseline", "vs_rho"))
    if not isinstance(vs, dict):
        raise CheckFailure("baseline stage (vs_rho) missing: run the rung with baseline.rho = true")
    targets = int(vs["targets"])
    if sol.get("count") != targets or len(items) != targets:
        raise CheckFailure(f"{sol.get('count')} solutions for {targets} targets")
    for item in items:
        if not item.get("verified"):
            raise CheckFailure(f"target {item.get('index')} not verified by the descent")
        expected = item.get("expected")
        if expected not in (None, "not_constructed") and item.get("recovered") != expected:
            raise CheckFailure(f"target {item.get('index')}: recovered {item.get('recovered')!r} != planted {expected!r}")
    ic_verified = sum(1 for item in items if item.get("verified"))
    if int(vs["ic"]["verified"]) != targets or ic_verified != targets:
        raise CheckFailure(f"ic verified {vs['ic']['verified']} of {targets}")
    rows = vs.get("targets_detail") or []
    if len(rows) != targets:
        raise CheckFailure(f"rho has {len(rows)} rows for {targets} targets")
    for row in rows:
        if not row.get("verified"):
            raise CheckFailure(f"target {row.get('index')} not verified by rho")
    rho_verified = int(vs["rho"]["verified"])
    if rho_verified != targets:
        raise CheckFailure(f"rho verified {rho_verified} of {targets}")
    return targets, ic_verified, rho_verified


# ── run ────────────────────────────────────────────────────────────────────


def _git_commit(cwd: Path) -> dict[str, Any]:
    try:
        head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=cwd, capture_output=True, text=True, check=True).stdout.strip()
        dirty = subprocess.run(["git", "status", "--porcelain"], cwd=cwd, capture_output=True, text=True, check=True).stdout.strip() != ""
        return {"commit": head, "dirty": dirty}
    except (OSError, subprocess.CalledProcessError):
        return {"commit": None, "dirty": None}


def _blake3_or_sha256(path: Path) -> dict[str, str]:
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 16), b""):
            digest.update(chunk)
    return {"sha256": digest.hexdigest()}


def _host() -> dict[str, Any]:
    return {
        "platform": platform.platform(),
        "machine": platform.machine(),
        "cpu_count": os.cpu_count(),
        "python": platform.python_version(),
    }


def cmd_run(args: argparse.Namespace) -> int:
    out = Path(args.output)
    if out.exists() and any(out.iterdir()):
        print(f"refusing to write into non-empty {out}", file=sys.stderr)
        return 2
    out.mkdir(parents=True, exist_ok=True)
    ic = Path(args.ic)
    if not ic.is_file():
        print(f"ic binary not found: {ic}", file=sys.stderr)
        return 2
    manifest: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "started_at": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        "ic_binary": {"path": str(ic), **_blake3_or_sha256(ic)},
        "git": _git_commit(Path.cwd()),
        "host": _host(),
        "rungs": [],
    }
    failures = 0
    for params in args.params:
        params_path = Path(params)
        name = params_path.stem
        rung_dir = out / name
        rung_dir.mkdir()
        report = rung_dir / REPORT_FILE
        cmd = [str(ic), "workflow", "--params", str(params_path), "--dir", str(rung_dir / RUN_DIR), "--out", str(report), "--json"]
        print(f"[ic-e2e] {name}: {' '.join(cmd)}", flush=True)
        started = _dt.datetime.now(_dt.timezone.utc)
        proc = subprocess.run(cmd, capture_output=True, text=True)
        (rung_dir / STDERR_FILE).write_text(proc.stderr)
        elapsed = (_dt.datetime.now(_dt.timezone.utc) - started).total_seconds()
        entry = {
            "name": name,
            "params": str(params_path),
            "command": cmd,
            "exit_code": proc.returncode,
            "wall_seconds": elapsed,
            "report": str(report.relative_to(out)),
            "stderr": str((rung_dir / STDERR_FILE).relative_to(out)),
        }
        if proc.returncode != 0 or not report.is_file():
            failures += 1
            print(f"[ic-e2e] {name}: exit {proc.returncode}; stderr tail:\n{proc.stderr[-2000:]}", file=sys.stderr)
        else:
            print(f"[ic-e2e] {name}: complete in {elapsed:.1f}s", flush=True)
        manifest["rungs"].append(entry)
    manifest["finished_at"] = _dt.datetime.now(_dt.timezone.utc).isoformat()
    (out / MANIFEST_FILE).write_text(json.dumps(manifest, indent=1, sort_keys=True) + "\n")
    return 1 if failures else 0


# ── measurement extraction ────────────────────────────────────────────────


def measure(report: dict[str, Any]) -> dict[str, Any]:
    targets, ic_verified, rho_verified = verify_correctness(report)
    counters = counters_of(report)
    vs = _lookup(report, ("stages.baseline", "vs_rho"))
    return {
        "name": report.get("name"),
        "params_digest": report.get("params_digest"),
        "degree": report.get("degree"),
        "subgroup_order": vs["subgroup_order"],
        "subgroup_bits": math.log2(int(vs["subgroup_order"])),
        "targets": targets,
        "claim_boundary": vs.get("claim_boundary"),
        "ic_verified": ic_verified,
        "rho_verified": rho_verified,
        "counters": counters,
        "pair_table_tier": tier_of(report),
        "wall": wall_of(report),
        "rho_S": rho_s_of(report, counters),
        "ic_S": None,
        "ic_S_null_reason": (
            "pinned IC counters are in mixed units (group additions, table lookups, scalar "
            "multiplications) with no measured conversion recorded; AGENTS.md §8 says leave it null"
        ),
    }


def load_run(output: Path) -> tuple[dict[str, Any], dict[str, dict[str, Any]]]:
    manifest_path = output / MANIFEST_FILE
    if not manifest_path.is_file():
        raise CheckFailure(f"no {MANIFEST_FILE} under {output}")
    manifest = json.loads(manifest_path.read_text())
    reports: dict[str, dict[str, Any]] = {}
    for entry in manifest.get("rungs", []):
        report_path = output / entry["report"]
        if entry.get("exit_code") != 0 or not report_path.is_file():
            raise CheckFailure(f"{entry['params']}: ic exited {entry.get('exit_code')} without a report")
        reports[entry["params"]] = json.loads(report_path.read_text())
    if not reports:
        raise CheckFailure("run directory holds no rungs")
    return manifest, reports


# ── freeze ─────────────────────────────────────────────────────────────────


def cmd_freeze(args: argparse.Namespace) -> int:
    ref_out = Path(args.reference_out)
    if ref_out.exists():
        print(f"refusing to overwrite {ref_out}; a superseded reference stays as the 'before' mark", file=sys.stderr)
        return 2
    try:
        manifest, reports = load_run(Path(args.output))
        rungs = {params: measure(report) for params, report in reports.items()}
        for params, m in rungs.items():
            _rho_sanity(params, m, args.rho_sanity_factor)
    except CheckFailure as exc:
        print(f"freeze refused: {exc}", file=sys.stderr)
        return 1
    reference = {
        "schema_version": SCHEMA_VERSION,
        "frozen_at": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        "frozen_from": {
            "git": manifest.get("git"),
            "ic_binary": manifest.get("ic_binary"),
            "host": manifest.get("host"),
            "note": args.note,
        },
        "unit": (
            "counters are exact seeded counts pinned bit-for-bit; wall figures are same-host, "
            "same-process seconds and are gated only as the paired ratio rho/IC; rho_S is group "
            "additions per target per sqrt(r); ic_S is null (no measured unit conversion)"
        ),
        "rungs": rungs,
    }
    ref_out.parent.mkdir(parents=True, exist_ok=True)
    ref_out.write_text(json.dumps(reference, indent=1, sort_keys=True) + "\n")
    print(f"froze {len(rungs)} rung(s) into {ref_out}")
    return 0


# ── check ──────────────────────────────────────────────────────────────────


def _rho_sanity(params: str, m: dict[str, Any], factor: float) -> None:
    ratio = m["rho_S"]["measured_over_expected"]
    if not (1.0 / factor <= ratio <= factor):
        raise CheckFailure(
            f"{params}: rho baseline S = {m['rho_S']['measured']:.4f} is {ratio:.2f}x its expected "
            f"{m['rho_S']['expected']:.4f} (outside 1/{factor:g}..{factor:g}); a broken baseline makes every ratio meaningless"
        )


def check_rung(params: str, ref: dict[str, Any], m: dict[str, Any], tolerance: float, rho_factor: float) -> dict[str, Any]:
    problems: list[str] = []
    if m["params_digest"] != ref["params_digest"]:
        problems.append(f"parameter file changed (digest {m['params_digest'][:12]} != frozen {ref['params_digest'][:12]}); re-freeze deliberately")
    frozen_tier = ref.get("pair_table_tier")
    if frozen_tier is not None and m.get("pair_table_tier") != frozen_tier:
        problems.append(
            f"pair table tier changed ({frozen_tier} -> {m.get('pair_table_tier')}); "
            "the algorithm changed — freeze a new reference beside the old one"
        )
    drift = {}
    for label, frozen in ref["counters"].items():
        now = m["counters"].get(label)
        if now != frozen:
            drift[label] = {"frozen": frozen, "now": now, "ratio": (now / frozen) if (frozen and now is not None) else None}
    if drift:
        detail = ", ".join(f"{k}: {v['frozen']} -> {v['now']}" for k, v in drift.items())
        problems.append(f"counters drifted from the frozen reference ({detail}); the algorithm changed — freeze a new reference beside the old one")
    try:
        _rho_sanity(params, m, rho_factor)
    except CheckFailure as exc:
        problems.append(str(exc))
    ref_ratio = float(ref["wall"]["whole_process_ratio"])
    now_ratio = float(m["wall"]["whole_process_ratio"])
    floor = ref_ratio * (1.0 - tolerance)
    if now_ratio < floor:
        problems.append(
            f"end-to-end wall ratio rho/IC fell to {now_ratio:.3f} from frozen {ref_ratio:.3f} (floor {floor:.3f} at tolerance {tolerance:g})"
        )
    if ref["wall"]["whole_process_crossover"] and not m["wall"]["whole_process_crossover"]:
        problems.append("this rung crossed rho end to end when frozen and no longer does")
    if m["wall"]["whole_process_crossover"] != m["wall"]["reported_whole_process_crossover"]:
        problems.append(
            f"the report's whole_process_crossover verdict ({m['wall']['reported_whole_process_crossover']}) disagrees "
            f"with its own timings ({m['wall']['whole_process_crossover']}); the evidence is inconsistent"
        )
    return {
        "params": params,
        "name": m["name"],
        "degree": m["degree"],
        "subgroup_bits": m["subgroup_bits"],
        "targets": m["targets"],
        "ic_verified": m["ic_verified"],
        "rho_verified": m["rho_verified"],
        "counters_identical": not drift,
        "counter_drift": drift,
        "wall": {"frozen": ref["wall"], "now": m["wall"], "floor_ratio": floor, "tolerance": tolerance},
        "rho_S": m["rho_S"],
        "ic_S": None,
        "problems": problems,
        "ok": not problems,
    }


def _fmt_ratio(x: float) -> str:
    return "inf" if math.isinf(x) else f"{x:.2f}"


def render_markdown(summary: dict[str, Any]) -> str:
    lines = [
        "## ic end-to-end benchmark",
        "",
        f"Result: **{'PASS' if summary['ok'] else 'FAIL'}** — {summary['rungs_ok']} of {summary['rungs_total']} rung(s) passed.",
        "",
        "Rows are rungs; every figure is the whole method (select + collect + logs + pair table + descent) "
        "against the signed-Frobenius ρ on the same targets in the same process. Counters are pinned exactly; "
        "the wall ratio is a same-host paired ratio gated with a tolerance, never an absolute. "
        "Passing is engineering at most (AGENTS.md §3), never an advance.",
        "",
        "| rung | n | log₂ r | targets | IC verified | ρ verified | counters | ρ/IC whole-process (frozen → now, floor) | ρ/IC charged (advisory) | e2e crossover | result |",
        "|:--|--:|--:|--:|--:|--:|:--|:--|--:|:--|:--|",
    ]
    for r in summary["rungs"]:
        w = r["wall"]
        lines.append(
            f"| `{r['params']}` | {r['degree']} | {r['subgroup_bits']:.1f} | {r['targets']} | {r['ic_verified']}/{r['targets']} | "
            f"{r['rho_verified']}/{r['targets']} | {'identical' if r['counters_identical'] else 'DRIFT'} | "
            f"{_fmt_ratio(w['frozen']['whole_process_ratio'])} → {_fmt_ratio(w['now']['whole_process_ratio'])} (≥ {_fmt_ratio(w['floor_ratio'])}) | "
            f"{_fmt_ratio(w['now']['charged_ratio'])} | {'yes' if w['now']['whole_process_crossover'] else 'no'} | "
            f"{'pass' if r['ok'] else 'FAIL'} |"
        )
    lines.append("")
    lines.append(
        "ρ baseline sanity (group additions per target per √r, measured / expected √(π/4n)): "
        + "; ".join(f"`{r['params']}` {r['rho_S']['measured']:.3f} / {r['rho_S']['expected']:.3f}" for r in summary["rungs"])
        + "."
    )
    lines.append("")
    lines.append("IC `S` is null on every row: the pinned IC counters are in mixed units with no measured conversion recorded (AGENTS.md §8).")
    for r in summary["rungs"]:
        for p in r["problems"]:
            lines.append(f"- **{r['params']}**: {p}")
    if summary.get("errors"):
        lines.append("")
        for e in summary["errors"]:
            lines.append(f"- **error**: {e}")
    lines.append("")
    return "\n".join(lines)


def cmd_check(args: argparse.Namespace) -> int:
    summary: dict[str, Any] = {"schema_version": SCHEMA_VERSION, "rungs": [], "errors": [], "ok": False}
    try:
        reference = json.loads(Path(args.reference).read_text())
        if reference.get("schema_version") != SCHEMA_VERSION:
            raise CheckFailure(f"reference schema {reference.get('schema_version')} != {SCHEMA_VERSION}")
        manifest, reports = load_run(Path(args.output))
        summary["manifest"] = {k: manifest.get(k) for k in ("git", "ic_binary", "host", "started_at", "finished_at")}
        summary["reference"] = {"path": str(args.reference), "frozen_at": reference.get("frozen_at"), "frozen_from": reference.get("frozen_from")}
        missing = sorted(set(reference["rungs"]) - set(reports))
        if missing:
            raise CheckFailure(f"reference rung(s) not run: {', '.join(missing)}")
        extra = sorted(set(reports) - set(reference["rungs"]))
        if extra:
            raise CheckFailure(f"rung(s) run without a frozen reference: {', '.join(extra)}")
        for params in reference["rungs"]:
            try:
                m = measure(reports[params])
            except CheckFailure as exc:
                summary["errors"].append(f"{params}: {exc}")
                continue
            summary["rungs"].append(check_rung(params, reference["rungs"][params], m, args.wall_ratio_tolerance, args.rho_sanity_factor))
    except (CheckFailure, OSError, json.JSONDecodeError, KeyError, TypeError, ValueError) as exc:
        summary["errors"].append(str(exc))
    summary["rungs_total"] = len(summary["rungs"]) + len(summary["errors"])
    summary["rungs_ok"] = sum(1 for r in summary["rungs"] if r["ok"])
    summary["ok"] = not summary["errors"] and bool(summary["rungs"]) and all(r["ok"] for r in summary["rungs"])
    text = render_markdown(summary)
    if args.summary_markdown:
        Path(args.summary_markdown).write_text(text)
    if args.summary_json:
        Path(args.summary_json).write_text(json.dumps(summary, indent=1, sort_keys=True) + "\n")
    print(text)
    return 0 if summary["ok"] else 1


# ── cli ────────────────────────────────────────────────────────────────────


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)

    run = sub.add_parser("run", help="run ic workflow on each parameter file into a fresh output directory")
    run.add_argument("--ic", required=True, help="path to the ic binary")
    run.add_argument("--output", required=True, help="fresh directory for reports, run state and the manifest")
    run.add_argument("--params", nargs="+", required=True, help="workflow parameter files (baseline.rho must be true)")
    run.set_defaults(func=cmd_run)

    check = sub.add_parser("check", help="gate a run directory against a frozen reference")
    check.add_argument("--output", required=True, help="run directory written by `run`")
    check.add_argument("--reference", required=True, help="frozen reference JSON")
    check.add_argument("--wall-ratio-tolerance", type=float, default=0.5, help="allowed relative fall of rho/IC whole-process ratio (default 0.5)")
    check.add_argument("--rho-sanity-factor", type=float, default=3.0, help="rho S must lie within this factor of sqrt(pi/4n) (default 3)")
    check.add_argument("--summary-markdown", help="write the summary table here (for $GITHUB_STEP_SUMMARY)")
    check.add_argument("--summary-json", help="write the machine-readable summary here")
    check.set_defaults(func=cmd_check)

    freeze = sub.add_parser("freeze", help="write a new frozen reference from a run directory; never overwrites")
    freeze.add_argument("--output", required=True, help="run directory written by `run`")
    freeze.add_argument("--reference-out", required=True, help="path of the new reference (must not exist)")
    freeze.add_argument("--rho-sanity-factor", type=float, default=3.0)
    freeze.add_argument("--note", default="", help="why this reference supersedes the previous one")
    freeze.set_defaults(func=cmd_freeze)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    sys.exit(main())
