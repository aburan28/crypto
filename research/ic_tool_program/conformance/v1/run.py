#!/usr/bin/env python3
"""The conformance suite's runner (IC_TOOL_PROGRAM.md §4, the C suite).

    python3 run.py --ic <binary> [--build-commit <sha>] [--out report.json]

Each case runs in a fresh temporary directory.  It passes only if:
- the process exits within its timeout;
- it never panics (exit status 101 is Rust's panic status);
- it meets every expectation in `cases.json`.

The report lists every case and the verdict.  The exit status is 0 only
if every case passed.
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
SUITE = HERE.parents[1] / "suite" / "v1"
PANIC_STATUS = 101


def expand(value: str, tmp: Path) -> str:
    return value.replace("{tmp}", str(tmp)).replace("{suite}", str(SUITE))


def set_path(doc: dict, dotted: str, value) -> None:
    *parents, last = dotted.split(".")
    node = doc
    for key in parents:
        node = node[key]
    node[last] = value


def materialise(files: dict, tmp: Path) -> None:
    for rel, spec in files.items():
        path = tmp / rel
        path.parent.mkdir(parents=True, exist_ok=True)
        if "copy" in spec:
            text = Path(expand(spec["copy"], tmp)).read_text()
            if "set" in spec:
                doc = json.loads(text)
                for dotted, value in spec["set"].items():
                    set_path(doc, dotted, value)
                text = json.dumps(doc, indent=1) + "\n"
            path.write_text(text)
        else:
            path.write_text(spec["text"])


def run_case(case: dict, ic: Path, build_commit: str | None) -> dict:
    tmp = Path(tempfile.mkdtemp(prefix=f"ic-conf-{case['id']}-"))
    try:
        materialise(case.get("files", {}), tmp)
        argv = [str(ic), *(expand(a, tmp) for a in case["argv"])]
        env = {**os.environ, **case.get("env", {})}
        try:
            proc = subprocess.run(argv, capture_output=True, text=True, env=env, timeout=case["timeout_s"])
        except subprocess.TimeoutExpired:
            return {"id": case["id"], "pass": False, "why": f"no exit within {case['timeout_s']} s"}
        expect, failures = case["expect"], []
        if proc.returncode == PANIC_STATUS:
            failures.append("panicked (exit status 101)")
        if expect.get("exit") == "zero" and proc.returncode != 0:
            failures.append(f"exit status {proc.returncode}, expected 0")
        if expect.get("exit") == "nonzero" and proc.returncode == 0:
            failures.append("exit status 0, expected a refusal")
        for needle in expect.get("stderr_contains", []):
            if needle not in proc.stderr:
                failures.append(f"stderr lacks {needle!r}")
        if "json_file" in expect:
            path = Path(expand(expect["json_file"], tmp))
            try:
                doc = json.loads(path.read_text())
            except (OSError, json.JSONDecodeError) as e:
                failures.append(f"no JSON report at {path.name}: {e}")
                doc = None
            if doc is not None:
                for key, want in expect.get("json_equals", {}).items():
                    if doc.get(key) != want:
                        failures.append(f"{key} = {doc.get(key)!r}, expected {want!r}")
                if expect.get("json_equals_build_commit") and build_commit is not None:
                    if doc.get("build_commit") != build_commit:
                        failures.append(f"build_commit = {doc.get('build_commit')!r}, expected {build_commit!r}")
        return {"id": case["id"], "pass": not failures, "why": failures, "exit": proc.returncode,
                "stderr_tail": proc.stderr.strip().splitlines()[-3:]}
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ic", required=True, type=Path)
    ap.add_argument("--build-commit")
    ap.add_argument("--out", type=Path)
    args = ap.parse_args()
    cases = json.loads((HERE / "cases.json").read_text())["cases"]
    results = [run_case(c, args.ic.resolve(), args.build_commit) for c in cases]
    report = {"binary": str(args.ic), "cases": len(results), "passed": sum(r["pass"] for r in results),
              "results": results}
    text = json.dumps(report, indent=1)
    if args.out:
        args.out.write_text(text + "\n")
    print(text)
    sys.exit(0 if report["passed"] == report["cases"] else 1)


if __name__ == "__main__":
    main()
