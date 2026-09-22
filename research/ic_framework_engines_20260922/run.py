#!/usr/bin/env python3
"""Run the matched engine suite frozen in ``suite.json``.

    python3 research/ic_framework_engines_20260922/run.py --run-id baseline_v1
    python3 research/ic_framework_engines_20260922/run.py --run-id baseline_v1 --parts A2-m2-20-22,H1-holdout

Every part writes one ``ic`` report into ``results/<run-id>/<part>.json``
beside a provenance record, and never overwrites anything: a part whose
report already exists is refused, so a run directory is append-only.
Several invocations with disjoint ``--parts`` may share one run
directory, which is how the frozen run spread its parts over the host's
cores.  The ``ic`` binary is built beforehand (``cargo build --release
--bin ic``); this script does not build, so the binary that ran is the
one whose hash the provenance records.

A candidate engine is compared against the frozen baseline by running
the same parts with the candidate added to each part's engine list
(``--add-engine name[:k=v,...]``) into a new run directory, then
``compare.py --run <new> --manifest manifest.json``.
"""

import argparse
import hashlib
import json
import os
import platform
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def sh(*cmd):
    try:
        return subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True, check=True).stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return None


def host():
    cpu = None
    try:
        for line in open("/proc/cpuinfo"):
            if line.startswith("model name"):
                cpu = line.split(":", 1)[1].strip()
                break
    except OSError:
        pass
    mem = None
    try:
        for line in open("/proc/meminfo"):
            if line.startswith("MemTotal"):
                mem = int(line.split()[1]) * 1024
                break
    except OSError:
        pass
    return {
        "cpu_model": cpu,
        "logical_cpus": os.cpu_count(),
        "memory_bytes": mem,
        "kernel": platform.release(),
        "python": platform.python_version(),
        "load_average_at_start": os.getloadavg() if hasattr(os, "getloadavg") else None,
    }


def stage_command(binary, part, defaults, out, extra_engines):
    cmd = [str(binary), "--json", "--out", str(out), "descent"]
    for engine in part["engines"] + extra_engines:
        cmd += ["--solver", engine]
    cmd += [
        "--cells", ",".join(part["cells"]),
        "--targets", str(part.get("targets", defaults["targets"])),
        "--repeats", str(part.get("repeats", defaults["repeats"])),
        "--families", ",".join(part.get("families", defaults["families"])),
        "--seed", str(part["seed"]),
        "--budget-seconds", str(part.get("budget_seconds", defaults["budget_seconds"])),
    ]
    return cmd


def bench_command(binary, part, out, extra_engines, scratch):
    sweep = HERE / part["sweep"]
    if extra_engines:
        # A candidate joins the frozen matrix; the frozen sweep file is
        # never edited, so the widened copy lives beside the report.
        doc = json.loads(sweep.read_text())
        doc["matrix"]["solver"] = doc["matrix"]["solver"] + extra_engines
        sweep = scratch / (part["id"] + ".sweep.json")
        sweep.write_text(json.dumps(doc, indent=1) + "\n")
    return [str(binary), "--json", "--out", str(out), "bench", "--sweep", str(sweep),
            "--rho-runs", str(part.get("rho_runs", 16))]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--parts", help="comma-separated part ids; default all, in suite order")
    ap.add_argument("--binary", default=str(ROOT / "target" / "release" / "ic"))
    ap.add_argument("--add-engine", action="append", default=[],
                    help="a candidate engine spec added to every part's engine list")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    suite = json.loads((HERE / "suite.json").read_text())
    set_vars = [v for v in suite["environment_must_be_unset"] if v in os.environ]
    if set_vars:
        sys.exit(f"refusing to run: {', '.join(set_vars)} set; these override the engines under test")
    binary = Path(args.binary)
    if not binary.exists():
        sys.exit(f"no binary at {binary}; build it with `cargo build --release --bin ic` first")

    parts = [("stage", p) for p in suite["stage_parts"]] + [("bench", p) for p in suite["bench_parts"]]
    if args.parts:
        wanted = args.parts.split(",")
        known = {p["id"] for _, p in parts}
        unknown = [w for w in wanted if w not in known]
        if unknown:
            sys.exit(f"unknown parts: {unknown}; known: {sorted(known)}")
        parts = [(k, p) for k, p in parts if p["id"] in wanted]

    out_dir = HERE / "results" / args.run_id
    out_dir.mkdir(parents=True, exist_ok=True)
    binary_sha256 = sha256_file(binary)
    commit = sh("git", "rev-parse", "HEAD")
    dirty = sh("git", "status", "--porcelain", "--untracked-files=no")
    for kind, part in parts:
        out = out_dir / f"{part['id']}.json"
        prov_path = out_dir / f"{part['id']}.provenance.json"
        if out.exists() or prov_path.exists():
            sys.exit(f"{out} exists; a run directory is append-only, use a new --run-id")
        if kind == "stage":
            cmd = stage_command(binary, part, suite["stage_defaults"], out, args.add_engine)
        else:
            cmd = bench_command(binary, part, out, args.add_engine, out_dir)
        print(" ".join(cmd), flush=True)
        if args.dry_run:
            continue
        started = time.time()
        log = out_dir / f"{part['id']}.stderr.log"
        with open(log, "x") as err:
            code = subprocess.run(cmd, cwd=ROOT, stdout=subprocess.DEVNULL, stderr=err).returncode
        prov = {
            "part": part["id"],
            "kind": kind,
            "suite_version": suite["version"],
            "suite_sha256": sha256_file(HERE / "suite.json"),
            "sweep_sha256": sha256_file(HERE / part["sweep"]) if kind == "bench" else None,
            "command": cmd,
            "added_engines": args.add_engine,
            "exit_code": code,
            "started_unix": started,
            "elapsed_seconds": time.time() - started,
            "source_commit": commit,
            "working_tree_dirty": bool(dirty),
            "binary_sha256": binary_sha256,
            "rustc": sh("rustc", "-V"),
            "cargo_profile": "release",
            "worker_threads": 1,
            "timeout_s": part.get("budget_seconds", suite["stage_defaults"]["budget_seconds"])
            if kind == "stage" else 120,
            "evidence_type": "measured_local",
            "host": host(),
            "concurrent_parts_note": "parts may run concurrently on separate cores; every stage cell is paired within one process, interleaved per target",
        }
        prov_path.write_text(json.dumps(prov, indent=1) + "\n")
        print(f"  -> exit {code} in {prov['elapsed_seconds']:.0f} s", flush=True)
        if code != 0:
            sys.exit(f"part {part['id']} failed; see {log}")


if __name__ == "__main__":
    main()
