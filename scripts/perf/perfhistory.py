#!/usr/bin/env python3
"""Performance history: one snapshot per commit, a trend chart, a PR check.

    perfhistory.py record --binary perfbench --ecbench ecbench --out docs/perf/history
    perfhistory.py record --ref REV --out DIR            # build REV (perfbench + ecbench) and record it
    perfhistory.py render --history DIR --md HISTORY.md --svg history.svg --json history.json
    perfhistory.py check  --binary perfbench --ecbench ecbench --history DIR [--threshold 0.10]

`perfindex.py` answers "is B faster than A?" for one pair of revisions.
This script answers "how has each kernel moved over the commits we
measured?": every snapshot records, for each `perfbench` kernel and for
each frozen `ecbench exec` child in `docs/perf/ecbench-children.json`,

  * the callgrind instruction count of the measured region (`Ir`):
    deterministic, host-independent up to the CPU features the code
    detects at run time, the number the chart is drawn from;
  * the wall time of one native run (median of a few samples): noisy on a
    shared CI host, recorded beside `Ir`, never charted as progress;
  * the kernel's fingerprint: when it changes, the kernel computed
    something else, the step is not a speedup, and the chained index
    takes no credit for it ("re-based").

The chained index of a kernel at snapshot t is the product over
consecutive snapshot pairs of Ir(previous) / Ir(current), with factor 1
where the fingerprint changed or the kernel was absent.  Areas are the
geometric mean of their kernels' indices, the overall index the weighted
geometric mean of areas (`docs/perf/weights.json`), exactly the formula of
`docs/perf/PERFORMANCE_INDEX.md` applied to Ir along the history instead
of across one pair.  Only the Python standard library is used.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import math
import os
import shutil
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import perfindex  # noqa: E402  (same directory)

SCHEMA = "perf.history/v1"
ECBENCH_AREA = "ecbench"
MOVER = 0.05  # a single-step change below this is within build noise
PALETTE = {
    "index": "#111111",
    "gf2_la": "#4e79a7",
    "bool_gb": "#f28e2b",
    "fp_gb": "#e15759",
    "sat": "#76b7b2",
    "pdp": "#59a14f",
    "field_ec": "#edc948",
    "relation": "#b07aa1",
    "dlp": "#ff9da7",
    ECBENCH_AREA: "#9c755f",
}


# ----------------------------------------------------------------------------
# git


def git(root: Path, *args: str) -> str:
    return subprocess.run(
        ["git", "-C", str(root), *args], capture_output=True, text=True, check=True
    ).stdout.strip()


def rev_facts(root: Path, rev: str) -> dict:
    sha = git(root, "rev-parse", rev)
    parents = git(root, "rev-list", "--parents", "-n", "1", sha).split()[1:]
    committed = git(root, "show", "-s", "--format=%cI", sha)
    subject = git(root, "show", "-s", "--format=%s", sha)
    return {
        "rev": sha,
        "rev_short": sha[:12],
        "parents": parents,
        "committed_at": committed,
        "subject": subject,
    }


# ----------------------------------------------------------------------------
# measuring


def wall_one(binary: str, kid: str, samples: int, max_seconds: float) -> dict:
    env = dict(os.environ, RAYON_NUM_THREADS="1")
    cmd = [
        binary, "run", "--filter", kid, "--exact",
        "--samples", str(samples), "--max-seconds", str(max_seconds),
    ]
    (row,) = perfindex.run_json(cmd, env)
    return row


def measure_perfbench(binary: str, filters: list[str], full: bool, wall_samples: int,
                      max_seconds: float, with_wall: bool, log) -> list[dict]:
    rows = []
    kernels = perfindex.list_kernels(binary, full, filters)
    for i, k in enumerate(kernels, 1):
        t0 = time.time()
        ir_row = perfindex.instr_one(binary, k["id"])
        row = {
            "id": k["id"],
            "area": k["area"],
            "tier": k.get("tier", "Quick"),
            "fingerprint": ir_row["fingerprint"],
            "ir": ir_row["ir"],
            "wall_ns": None,
        }
        if with_wall:
            w = wall_one(binary, k["id"], wall_samples, max_seconds)
            if w["fingerprint"] != row["fingerprint"]:
                raise RuntimeError(f"{k['id']}: fingerprint differs between callgrind and native runs")
            row["wall_ns"] = w["median_ns"]
        rows.append(row)
        log(f"[{i}/{len(kernels)}] {k['id']}: Ir={row['ir']} wall={_ms(row['wall_ns'])} ({time.time() - t0:.1f}s)")
    return rows


def ecbench_child_once(ecbench: str, child: dict, workdir: Path, callgrind: bool) -> tuple[dict, dict | None]:
    """Run one measured child; return (child_output, callgrind_ir or None)."""
    inp = json.dumps(child["input"])
    env = dict(os.environ)
    if callgrind:
        env["ECBENCH_CALLGRIND_SOLVE"] = "1"
        prefix = workdir / "cg.out"
        for f in workdir.glob("cg.out*"):
            f.unlink()
        cmd = ["valgrind", "--tool=callgrind", f"--callgrind-out-file={prefix}", ecbench, "exec"]
    else:
        cmd = [ecbench, "exec"]
    proc = subprocess.run(cmd, input=inp, capture_output=True, text=True, env=env)
    if proc.returncode != 0:
        raise RuntimeError(f"{child['id']}: ecbench exec exited {proc.returncode}:\n{proc.stderr[-3000:]}")
    out = json.loads(proc.stdout.strip().splitlines()[-1])
    if out.get("error"):
        raise RuntimeError(f"{child['id']}: child error: {out['error']}")
    ir = None
    if callgrind:
        p = subprocess.run(
            [ecbench, "callgrind-ir", "--prefix", str(prefix)], capture_output=True, text=True, check=True
        )
        ir = json.loads(p.stdout)
    return out, ir


def child_fingerprint(out: dict) -> str:
    """What the child computed: the recovered scalar and every counted phase."""
    rep = out["report"]
    phases = [
        {k: v for k, v in sorted(ph.items()) if not k.endswith("_ns") and k != "wall_ns"}
        for ph in rep.get("phases", [])
    ]
    blob = json.dumps({"recovered": rep.get("recovered"), "phases": phases}, sort_keys=True)
    import hashlib

    return hashlib.sha256(blob.encode()).hexdigest()[:16]


def measure_ecbench(ecbench: str, children_path: Path, wall_samples: int, with_wall: bool, log) -> list[dict]:
    doc = json.loads(children_path.read_text())
    rows = []
    with tempfile.TemporaryDirectory() as td:
        work = Path(td)
        for child in doc["children"]:
            t0 = time.time()
            out, ir = ecbench_child_once(ecbench, child, work, callgrind=True)
            rec = out["report"].get("recovered")
            if rec != child["expected_recovered"]:
                raise RuntimeError(f"{child['id']}: recovered {rec}, expected {child['expected_recovered']}")
            row = {
                "id": child["id"],
                "area": ECBENCH_AREA,
                "tier": "Quick",
                "fingerprint": child_fingerprint(out),
                "ir": ir["solve_ir"],
                "prework_ir": ir["prework_ir"],
                "wall_ns": None,
            }
            if with_wall:
                times = []
                for _ in range(max(1, wall_samples)):
                    t = time.perf_counter()
                    o2, _ = ecbench_child_once(ecbench, child, work, callgrind=False)
                    times.append(int((time.perf_counter() - t) * 1e9))
                    if child_fingerprint(o2) != row["fingerprint"]:
                        raise RuntimeError(f"{child['id']}: counted output differs between runs")
                row["wall_ns"] = int(statistics.median(times))
            rows.append(row)
            log(f"{child['id']}: solve Ir={row['ir']} wall={_ms(row['wall_ns'])} ({time.time() - t0:.1f}s)")
    return rows


def _ms(ns):
    return "n/a" if ns is None else f"{ns / 1e6:.3f}ms"


def valgrind_version() -> str:
    try:
        return subprocess.run(["valgrind", "--version"], capture_output=True, text=True).stdout.strip()
    except OSError:
        return ""


# ----------------------------------------------------------------------------
# record


def build_ref(root: Path, ref: str, keep_dir: Path, log) -> tuple[str, str, Path]:
    """Build perfbench (this tree's harness) and ecbench at REF in a sparse worktree."""
    crate = perfindex.crate_dir(root)
    src = keep_dir / "src"
    target = keep_dir / "target"
    subprocess.run(["git", "-C", str(root), "worktree", "add", "--no-checkout", "--detach", str(src), ref], check=True)
    dirs = perfindex.sparse_dirs(crate)
    if dirs:
        subprocess.run(["git", "-C", str(src), "sparse-checkout", "set", *dirs, "research/ecbench_calibration_20261002"], check=True)
    subprocess.run(["git", "-C", str(src), "checkout", "--detach", ref], check=True)
    harness_rel = Path(crate) / "examples" / perfindex.EXAMPLE
    dst = src / harness_rel
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(root / harness_rel, dst)
    env = dict(os.environ, CARGO_TARGET_DIR=str(target))
    t0 = time.time()
    subprocess.run(
        ["cargo", "build", "--release", "--example", perfindex.EXAMPLE, "--bin", "ecbench"],
        cwd=src / crate, env=env, check=True,
    )
    log(f"built {ref} in {time.time() - t0:.0f}s")
    return str(target / "release" / "examples" / perfindex.EXAMPLE), str(target / "release" / "ecbench"), src


def snapshot_name(facts: dict) -> str:
    t = dt.datetime.fromisoformat(facts["committed_at"]).astimezone(dt.timezone.utc)
    return f"{t.strftime('%Y%m%dT%H%M%SZ')}-{facts['rev_short']}.json"


def cmd_record(a: argparse.Namespace) -> None:
    root = perfindex.repo_root()
    log = lambda s: print(s, file=sys.stderr, flush=True)  # noqa: E731
    out_dir = Path(a.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    tmp_build = None
    try:
        if a.ref:
            tmp_build = Path(tempfile.mkdtemp(prefix="perfhistory-build-"))
            perfbench, ecbench, src = build_ref(root, a.ref, tmp_build, log)
            facts = rev_facts(root, a.ref)
            harness = perfindex.sha256_tree(root / perfindex.crate_dir(root) / "examples" / perfindex.EXAMPLE)
        else:
            if not a.binary:
                raise SystemExit("--binary is required without --ref")
            perfbench, ecbench = a.binary, a.ecbench
            facts = rev_facts(root, a.rev or "HEAD")
            harness = perfindex.sha256_tree(root / perfindex.crate_dir(root) / "examples" / perfindex.EXAMPLE)
        host = perfindex.host_manifest()
        kernels = measure_perfbench(perfbench, a.filter, a.full, a.wall_samples, a.max_seconds, not a.no_wall, log)
        if ecbench and a.children and Path(a.children).exists():
            kernels += measure_ecbench(ecbench, Path(a.children), a.wall_samples, not a.no_wall, log)
        snap = {
            "schema": SCHEMA,
            **facts,
            "recorded_at": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
            "note": a.note or "",
            "harness_sha256": harness,
            "valgrind": valgrind_version(),
            "host": {k: v for k, v in host.items() if not k.startswith("loadavg")},
            "wall_samples": 0 if a.no_wall else a.wall_samples,
            "kernels": kernels,
        }
        path = out_dir / snapshot_name(facts)
        path.write_text(json.dumps(snap, indent=1) + "\n")
        log(f"wrote {path} ({len(kernels)} kernels)")
        print(path)
    finally:
        if tmp_build is not None and not a.keep:
            subprocess.run(["git", "-C", str(root), "worktree", "remove", "--force", str(tmp_build / "src")],
                           capture_output=True)
            shutil.rmtree(tmp_build, ignore_errors=True)


# ----------------------------------------------------------------------------
# the chained index


def load_history(history: Path) -> list[dict]:
    snaps = [json.loads(p.read_text()) for p in sorted(history.glob("*.json"))]
    snaps = [s for s in snaps if s.get("schema") == SCHEMA]
    snaps.sort(key=lambda s: (s["committed_at"], s["recorded_at"]))
    return snaps


def chained(snaps: list[dict], weights: dict) -> dict:
    """Per-kernel chained indices, area indices and the overall index at every snapshot."""
    ids: dict[str, str] = {}  # id -> area
    for s in snaps:
        for k in s["kernels"]:
            ids[k["id"]] = k["area"]
    per = {kid: [] for kid in ids}  # list of {ir, wall_ns, fingerprint, index, step, rebased}
    for t, s in enumerate(snaps):
        rows = {k["id"]: k for k in s["kernels"]}
        for kid in ids:
            prev = per[kid][-1] if per[kid] else None
            k = rows.get(kid)
            if k is None:
                per[kid].append({"present": False, "index": prev["index"] if prev else 1.0, "step": None,
                                 "rebased": False, "ir": None, "wall_ns": None, "fingerprint": None})
                continue
            step, rebased = None, False
            index = prev["index"] if prev else 1.0
            if prev and prev["present"]:
                if prev["fingerprint"] == k["fingerprint"] and prev["ir"] and k["ir"]:
                    step = prev["ir"] / k["ir"]
                    index *= step
                elif prev["fingerprint"] != k["fingerprint"]:
                    rebased = True
            per[kid].append({"present": True, "index": index, "step": step, "rebased": rebased,
                             "ir": k["ir"], "wall_ns": k.get("wall_ns"), "fingerprint": k["fingerprint"]})
    areas = sorted({a for a in ids.values()})
    area_index = {a: [] for a in areas}
    overall = []
    for t in range(len(snaps)):
        logs_by_area = {a: [] for a in areas}
        for kid, area in ids.items():
            logs_by_area[area].append(math.log(per[kid][t]["index"]))
        for a in areas:
            area_index[a].append(math.exp(statistics.fmean(logs_by_area[a])) if logs_by_area[a] else 1.0)
        present = [a for a in areas if logs_by_area[a]]
        w = {a: float(weights.get(a, 1.0)) for a in present}
        tot = sum(w.values()) or 1.0
        overall.append(math.exp(sum(w[a] / tot * math.log(area_index[a][t]) for a in present)))
    return {"ids": ids, "per": per, "areas": areas, "area_index": area_index, "overall": overall}


# ----------------------------------------------------------------------------
# render


def fmt_x(x: float) -> str:
    return f"{x:.3f}×"


def svg_chart(snaps: list[dict], ch: dict, title: str) -> str:
    n = len(snaps)
    W, H = 960, 420
    ml, mr, mt, mb = 70, 230, 40, 70
    pw, ph = W - ml - mr, H - mt - mb
    series = [("index", ch["overall"])] + [(a, ch["area_index"][a]) for a in ch["areas"]]
    vals = [math.log2(v) for _, ys in series for v in ys if v > 0]
    lo, hi = (min(vals), max(vals)) if vals else (0.0, 0.0)
    lo, hi = min(lo, 0.0), max(hi, 0.0)
    pad = max(0.25, (hi - lo) * 0.08)
    lo, hi = lo - pad, hi + pad

    def X(i):
        return ml + (pw * i / (n - 1) if n > 1 else pw / 2)

    def Y(v):
        return mt + ph * (1 - (math.log2(v) - lo) / (hi - lo))

    out = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" font-family="ui-sans-serif, system-ui, sans-serif" font-size="12">',
           f'<rect width="{W}" height="{H}" fill="#ffffff"/>',
           f'<text x="{ml}" y="22" font-size="15" font-weight="600" fill="#111">{title}</text>',
           f'<text x="{ml}" y="{mt - 6}" fill="#555">chained Ir index relative to the first snapshot (log scale; above 1 = fewer instructions). wall time is not charted.</text>']
    # gridlines at powers of two
    k0, k1 = math.floor(lo), math.ceil(hi)
    for k in range(k0, k1 + 1):
        y = mt + ph * (1 - (k - lo) / (hi - lo))
        if mt <= y <= mt + ph:
            out.append(f'<line x1="{ml}" x2="{ml + pw}" y1="{y:.1f}" y2="{y:.1f}" stroke="#e5e5e5"/>')
            label = f"{2 ** k:g}×" if k >= 0 else f"1/{2 ** (-k):g}×"
            out.append(f'<text x="{ml - 8}" y="{y + 4:.1f}" text-anchor="end" fill="#555">{label}</text>')
    out.append(f'<line x1="{ml}" x2="{ml + pw}" y1="{Y(1.0):.1f}" y2="{Y(1.0):.1f}" stroke="#999" stroke-dasharray="4 3"/>')
    # x labels
    step = max(1, math.ceil(n / 8))
    for i, s in enumerate(snaps):
        if i % step == 0 or i == n - 1:
            x = X(i)
            out.append(f'<line x1="{x:.1f}" x2="{x:.1f}" y1="{mt + ph}" y2="{mt + ph + 5}" stroke="#555"/>')
            date = s["committed_at"][:10]
            out.append(f'<text x="{x:.1f}" y="{mt + ph + 18}" text-anchor="middle" fill="#555">{date}</text>')
            out.append(f'<text x="{x:.1f}" y="{mt + ph + 32}" text-anchor="middle" fill="#888" font-size="10">{s["rev_short"][:8]}</text>')
    # series
    for name, ys in series:
        color = PALETTE.get(name, "#333")
        width = 3 if name == "index" else 1.6
        pts = " ".join(f"{X(i):.1f},{Y(v):.1f}" for i, v in enumerate(ys))
        if n == 1:
            out.append(f'<circle cx="{X(0):.1f}" cy="{Y(ys[0]):.1f}" r="4" fill="{color}"/>')
        else:
            out.append(f'<polyline points="{pts}" fill="none" stroke="{color}" stroke-width="{width}"/>')
            for i, v in enumerate(ys):
                out.append(f'<circle cx="{X(i):.1f}" cy="{Y(v):.1f}" r="2.5" fill="{color}"/>')
    # legend with latest values
    lx, ly = ml + pw + 16, mt
    for j, (name, ys) in enumerate(series):
        color = PALETTE.get(name, "#333")
        y = ly + j * 18
        out.append(f'<line x1="{lx}" x2="{lx + 18}" y1="{y + 5}" y2="{y + 5}" stroke="{color}" stroke-width="{3 if name == "index" else 2}"/>')
        label = "overall index" if name == "index" else name
        out.append(f'<text x="{lx + 24}" y="{y + 9}" fill="#111" font-weight="{600 if name == "index" else 400}">{label} {fmt_x(ys[-1])}</text>')
    out.append(f'<text x="{ml}" y="{H - 8}" fill="#888" font-size="10">{n} snapshot{"s" if n != 1 else ""}; generated by scripts/perf/perfhistory.py render</text>')
    out.append("</svg>")
    return "\n".join(out) + "\n"


def pct(a, b):
    """b relative to a, as a signed percentage of change in Ir (negative = fewer instructions)."""
    if not a or not b:
        return None
    return (b / a - 1.0) * 100.0


def render_md(snaps: list[dict], ch: dict, svg_rel: str | None) -> str:
    n = len(snaps)
    last = snaps[-1]
    L = []
    L.append("# Performance history\n")
    L.append("Generated by `scripts/perf/perfhistory.py render` from the snapshots in "
             "`docs/perf/history/`; do not edit by hand.  One snapshot per measured commit: "
             "callgrind instruction counts (`Ir`) of every quick `perfbench` kernel and of the frozen "
             "`ecbench` children in `docs/perf/ecbench-children.json`, with a wall-time reading beside "
             "each.  The index chains consecutive `Ir` ratios and takes no credit across a fingerprint "
             "change (the kernel computed something else: **re-based**).  Formula: "
             "`docs/perf/PERFORMANCE_INDEX.md`; what this does and does not mean: engineering class, "
             "`AGENTS.md` §3.  **Noise floor:** two builds of different sources move the `Ir` of "
             "untouched kernels by up to a few per cent (the release profile's 16 codegen units are "
             "partitioned afresh, and cross-module inlining changes with them), so a single step under "
             "about 5 % on one kernel is build noise, not a result; a trend over several snapshots, or "
             "an area index, is.\n")
    if svg_rel:
        L.append(f"![performance history]({svg_rel})\n")
    L.append(f"**Latest snapshot:** `{last['rev_short']}` ({last['committed_at'][:19]}Z) — {last['subject']}  ")
    L.append(f"**Snapshots:** {n}, from `{snaps[0]['rev_short']}` ({snaps[0]['committed_at'][:10]}).  ")
    L.append(f"**Overall chained Ir index since first snapshot:** **{fmt_x(ch['overall'][-1])}**  ")
    if last.get("note"):
        L.append(f"**Note on the latest snapshot:** {last['note']}  ")
    L.append("")
    # areas
    L.append("## Areas\n")
    L.append("| area | kernels | index since first | step from previous |")
    L.append("|---|---:|---:|---:|")
    for a in ch["areas"]:
        ks = [kid for kid, ar in ch["ids"].items() if ar == a]
        cur = ch["area_index"][a][-1]
        prev = ch["area_index"][a][-2] if n > 1 else cur
        L.append(f"| {a} | {len(ks)} | {fmt_x(cur)} | {fmt_x(cur / prev) if prev else 'n/a'} |")
    L.append(f"| **overall** | {len(ch['ids'])} | **{fmt_x(ch['overall'][-1])}** | "
             f"{fmt_x(ch['overall'][-1] / ch['overall'][-2]) if n > 1 else 'n/a'} |\n")
    # snapshots
    L.append("## Snapshots\n")
    L.append("| # | committed | commit | overall | " + " | ".join(ch["areas"]) + " | kernels | host |")
    L.append("|---:|---|---|---:|" + "---:|" * len(ch["areas"]) + "---:|---|")
    for t, s in enumerate(snaps):
        cpu = s.get("host", {}).get("cpu", "")
        cpu = (cpu[:28] + "…") if len(cpu) > 29 else cpu
        L.append(f"| {t + 1} | {s['committed_at'][:10]} | `{s['rev_short'][:8]}` | {fmt_x(ch['overall'][t])} | "
                 + " | ".join(fmt_x(ch["area_index"][a][t]) for a in ch["areas"])
                 + f" | {len(s['kernels'])} | {cpu} |")
    L.append("")
    # movers
    if n > 1:
        moves = []
        for kid, hist in ch["per"].items():
            h = hist[-1]
            if h["step"]:
                moves.append((h["step"], kid, hist[-2]["ir"], h["ir"]))
        moves.sort(reverse=True)
        up = [m for m in moves if m[0] > 1 + MOVER][:12]
        down = sorted([m for m in moves if m[0] < 1 / (1 + MOVER)])[:12]
        L.append("## Movers in the latest step\n")
        L.append(f"Kernels whose `Ir` changed by more than {MOVER * 100:.0f} % against the previous snapshot with the "
                 "same fingerprint (the build-noise floor is a few per cent, see above).\n")
        if up:
            L.append("Fewer instructions:\n")
            L.append("| kernel | previous Ir | latest Ir | step |")
            L.append("|---|---:|---:|---:|")
            for s_, kid, a, b in up:
                L.append(f"| `{kid}` | {a:,} | {b:,} | {fmt_x(s_)} |")
            L.append("")
        if down:
            L.append("More instructions:\n")
            L.append("| kernel | previous Ir | latest Ir | step |")
            L.append("|---|---:|---:|---:|")
            for s_, kid, a, b in down:
                L.append(f"| `{kid}` | {a:,} | {b:,} | {fmt_x(s_)} |")
            L.append("")
        if not up and not down:
            L.append("None.\n")
    # per kernel
    L.append("## Kernels at the latest snapshot\n")
    L.append("`Δ prev` is the latest step's change in `Ir` (negative = fewer instructions); `index` is the chained "
             "index since the first snapshot that has the kernel; `wall` is one native median on the recording "
             "host and is informational only.\n")
    for a in ch["areas"]:
        L.append(f"### {a}\n")
        L.append("| kernel | Ir | wall | Δ prev | index | note |")
        L.append("|---|---:|---:|---:|---:|---|")
        for kid in sorted(k for k, ar in ch["ids"].items() if ar == a):
            hist = ch["per"][kid]
            h = hist[-1]
            if not h["present"]:
                L.append(f"| `{kid}` | — | — | — | {fmt_x(h['index'])} | absent in latest snapshot |")
                continue
            prev = hist[-2] if n > 1 else None
            d = pct(prev["ir"], h["ir"]) if prev and prev["present"] and not h["rebased"] else None
            note = "re-based (fingerprint changed)" if h["rebased"] else ("new" if n > 1 and prev and not prev["present"] else "")
            L.append(f"| `{kid}` | {h['ir']:,} | {_ms(h['wall_ns'])} | {'' if d is None else f'{d:+.1f} %'} | {fmt_x(h['index'])} | {note} |")
        L.append("")
    L.append("## How to add a point\n")
    L.append("```bash\ncargo build --release --example perfbench --bin ecbench\n"
             "python3 scripts/perf/perfhistory.py record --binary target/release/examples/perfbench \\\n"
             "    --ecbench target/release/ecbench --out docs/perf/history\n"
             "python3 scripts/perf/perfhistory.py render --history docs/perf/history \\\n"
             "    --md docs/perf/HISTORY.md --svg docs/perf/history.svg --json docs/perf/history.json\n```\n")
    L.append("`.github/workflows/perf-history.yml` does this for every push to `main` and checks each pull "
             "request against the latest snapshot.\n")
    return "\n".join(L)


def render_json(snaps: list[dict], ch: dict) -> dict:
    """A compact machine-readable view for agents: latest per-kernel numbers and the area series."""
    return {
        "schema": "perf.history_summary/v1",
        "snapshots": [{"rev": s["rev"], "committed_at": s["committed_at"], "subject": s["subject"],
                       "overall": ch["overall"][t],
                       "areas": {a: ch["area_index"][a][t] for a in ch["areas"]}}
                      for t, s in enumerate(snaps)],
        "kernels": {kid: {"area": ch["ids"][kid],
                          "latest_ir": ch["per"][kid][-1]["ir"],
                          "latest_wall_ns": ch["per"][kid][-1]["wall_ns"],
                          "index": ch["per"][kid][-1]["index"],
                          "series_ir": [h["ir"] for h in ch["per"][kid]],
                          "rebased_at": [snaps[t]["rev"] for t, h in enumerate(ch["per"][kid]) if h["rebased"]]}
                    for kid in sorted(ch["ids"])},
    }


def cmd_render(a: argparse.Namespace) -> None:
    snaps = load_history(Path(a.history))
    if not snaps:
        raise SystemExit(f"no snapshots in {a.history}")
    ch = chained(snaps, perfindex.load_weights(a.weights))
    if a.svg:
        Path(a.svg).write_text(svg_chart(snaps, ch, "crypto: perfbench + ecbench instruction-count history"))
    svg_rel = None
    if a.svg and a.md:
        svg_rel = os.path.relpath(a.svg, Path(a.md).parent)
    if a.md:
        Path(a.md).write_text(render_md(snaps, ch, svg_rel))
    if a.json:
        Path(a.json).write_text(json.dumps(render_json(snaps, ch), indent=1) + "\n")
    print(f"{len(snaps)} snapshots; overall index {fmt_x(ch['overall'][-1])}; "
          + ", ".join(f"{a_}={fmt_x(ch['area_index'][a_][-1])}" for a_ in ch["areas"]))


# ----------------------------------------------------------------------------
# check


def cmd_check(a: argparse.Namespace) -> None:
    log = lambda s: print(s, file=sys.stderr, flush=True)  # noqa: E731
    snaps = load_history(Path(a.history))
    if not snaps:
        raise SystemExit(f"no snapshots in {a.history}")
    base = snaps[-1]
    rows = measure_perfbench(a.binary, a.filter, a.full, 1, 0.5, False, log)
    if a.ecbench and a.children and Path(a.children).exists():
        rows += measure_ecbench(a.ecbench, Path(a.children), 1, False, log)
    prev = {k["id"]: k for k in base["kernels"]}
    cur = {k["id"]: k for k in rows}
    regress, improve, rebased, new, gone, same = [], [], [], [], [], 0
    logs = []
    for kid in sorted(set(prev) | set(cur)):
        p, c = prev.get(kid), cur.get(kid)
        if p is None:
            new.append(kid)
        elif c is None:
            gone.append(kid)
        elif p["fingerprint"] != c["fingerprint"]:
            rebased.append(kid)
        else:
            r = c["ir"] / p["ir"]
            logs.append(math.log(r))
            if r > 1 + a.threshold:
                regress.append((r, kid, p["ir"], c["ir"]))
            elif r < 1 - a.threshold:
                improve.append((r, kid, p["ir"], c["ir"]))
            else:
                same += 1
    index = math.exp(-statistics.fmean(logs)) if logs else 1.0
    L = [f"## perf history check against `{base['rev_short']}` ({base['committed_at'][:10]})\n",
         f"Geometric-mean Ir index of this tree over the latest snapshot (same-fingerprint kernels): **{fmt_x(index)}** "
         f"(above 1 = fewer instructions). Threshold ±{a.threshold * 100:.0f} %.\n",
         f"- {len(improve)} kernels with fewer instructions, {len(regress)} with more, {same} within threshold, "
         f"{len(rebased)} re-based (fingerprint changed), {len(new)} new, {len(gone)} missing.\n"]
    for title, items in (("More instructions (regressions)", sorted(regress, reverse=True)),
                         ("Fewer instructions", sorted(improve))):
        if items:
            L.append(f"### {title}\n\n| kernel | snapshot Ir | this tree Ir | ratio |\n|---|---:|---:|---:|")
            for r, kid, pi, ci in items:
                L.append(f"| `{kid}` | {pi:,} | {ci:,} | {r:.3f} |")
            L.append("")
    if rebased:
        L.append("### Re-based (fingerprint changed; not a speed comparison)\n\n" + "\n".join(f"- `{k}`" for k in rebased) + "\n")
    if new:
        L.append("### New kernels\n\n" + "\n".join(f"- `{k}`" for k in new) + "\n")
    if gone:
        L.append("### Missing kernels\n\n" + "\n".join(f"- `{k}`" for k in gone) + "\n")
    md = "\n".join(L)
    if a.md:
        Path(a.md).write_text(md)
    print(md)
    if a.fail_on_regression and regress:
        sys.exit(1)


# ----------------------------------------------------------------------------


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("record", help="measure one revision and write a snapshot")
    r.add_argument("--binary", help="perfbench binary for the working tree")
    r.add_argument("--ecbench", help="ecbench binary for the working tree")
    r.add_argument("--ref", help="build and measure this revision instead (sparse worktree, this tree's harness)")
    r.add_argument("--rev", help="record the snapshot under this revision (default HEAD)")
    r.add_argument("--children", default="docs/perf/ecbench-children.json")
    r.add_argument("--out", required=True, help="snapshot directory")
    r.add_argument("--note", default="")
    r.add_argument("--filter", action="append", default=[])
    r.add_argument("--full", action="store_true", help="include Tier::Full kernels")
    r.add_argument("--wall-samples", type=int, default=5)
    r.add_argument("--max-seconds", type=float, default=1.0)
    r.add_argument("--no-wall", action="store_true", help="skip the native wall-time reading")
    r.add_argument("--keep", action="store_true", help="keep the --ref build directory")
    r.set_defaults(func=cmd_record)

    d = sub.add_parser("render", help="write HISTORY.md, history.svg and history.json from the snapshots")
    d.add_argument("--history", required=True)
    d.add_argument("--md")
    d.add_argument("--svg")
    d.add_argument("--json")
    d.add_argument("--weights")
    d.set_defaults(func=cmd_render)

    c = sub.add_parser("check", help="compare the working tree's Ir against the latest snapshot")
    c.add_argument("--binary", required=True)
    c.add_argument("--ecbench")
    c.add_argument("--children", default="docs/perf/ecbench-children.json")
    c.add_argument("--history", required=True)
    c.add_argument("--threshold", type=float, default=0.10)
    c.add_argument("--filter", action="append", default=[])
    c.add_argument("--full", action="store_true")
    c.add_argument("--md")
    c.add_argument("--fail-on-regression", action="store_true")
    c.set_defaults(func=cmd_check)

    a = p.parse_args()
    a.func(a)


if __name__ == "__main__":
    main()
