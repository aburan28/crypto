#!/usr/bin/env python3
"""Conformance suite v2's runner: v1's cases, then v2's
(`research/ic_tool_program/design/schema-v2.md` §9).

    python3 run.py --ic <binary> --through B1 [--build-commit <sha>] [--out report.json]

`--through` names the newest step whose cases must pass. v1's cases are
step B0's. A step is accepted only when its cases and every earlier
step's pass. The expectations are v1's, plus:
- `exit` as an exact number;
- `json_paths`, values at dotted paths;
- `json_contains`, an object some element of a list must contain;
- `same_outputs_as`, a second command whose report must agree on the
  listed paths.

The exit status is 0 only if every selected case passed.
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
V1 = HERE.parent / "v1"
SUITE = HERE.parents[1] / "suite" / "v1"
PANIC_STATUS = 101
STEPS = ("B0", "B1", "B2", "B3", "B4", "B5", "B6", "B7")
MISSING = object()


def expand(value: str, tmp: Path) -> str:
    return (value.replace("{tmp}", str(tmp)).replace("{suite}", str(SUITE))
            .replace("{cases}", str(HERE / "params")))


def set_path(doc: dict, dotted: str, value) -> None:
    *parents, last = dotted.split(".")
    node = doc
    for key in parents:
        node = node[key]
    node[last] = value


def get_path(doc, dotted: str):
    node = doc
    for key in dotted.split("."):
        if isinstance(node, list) and key.isdigit() and int(key) < len(node):
            node = node[int(key)]
        elif isinstance(node, dict) and key in node:
            node = node[key]
        else:
            return MISSING
    return node


def contains(want, have) -> bool:
    """Whether `have` holds everything in `want`: objects key by key, anything else by equality."""
    if isinstance(want, dict):
        return isinstance(have, dict) and all(k in have and contains(v, have[k]) for k, v in want.items())
    return want == have


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


def load_report(path: Path, failures: list[str]):
    try:
        return json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as e:
        failures.append(f"no JSON report at {path.name}: {e}")
        return None


def run(argv: list[str], ic: Path, tmp: Path, env: dict, timeout: int):
    return subprocess.run([str(ic), *(expand(a, tmp) for a in argv)], capture_output=True, text=True,
                          env={**os.environ, **env}, timeout=timeout)


def check_exit(want, status: int, failures: list[str]) -> None:
    if status == PANIC_STATUS:
        failures.append("panicked (exit status 101)")
    if want == "zero" and status != 0:
        failures.append(f"exit status {status}, expected 0")
    elif want == "nonzero" and status == 0:
        failures.append("exit status 0, expected a refusal")
    elif isinstance(want, int) and status != want:
        failures.append(f"exit status {status}, expected {want}")


def check_report(doc: dict, expect: dict, build_commit: str | None, failures: list[str]) -> None:
    for key, want in expect.get("json_equals", {}).items():
        if doc.get(key) != want:
            failures.append(f"{key} = {doc.get(key)!r}, expected {want!r}")
    if expect.get("json_equals_build_commit") and build_commit is not None:
        if doc.get("build_commit") != build_commit:
            failures.append(f"build_commit = {doc.get('build_commit')!r}, expected {build_commit!r}")
    for dotted, want in expect.get("json_paths", {}).items():
        have = get_path(doc, dotted)
        if have is MISSING:
            failures.append(f"{dotted} is absent, expected {want!r}")
        elif have != want:
            failures.append(f"{dotted} = {have!r}, expected {want!r}")
    for dotted, wants in expect.get("json_contains", {}).items():
        have = get_path(doc, dotted)
        if not isinstance(have, list):
            failures.append(f"{dotted} is not a list")
            continue
        for want in wants:
            if not any(contains(want, item) for item in have):
                failures.append(f"{dotted} has no element containing {want!r}")


def run_case(case: dict, ic: Path, build_commit: str | None) -> dict:
    tmp = Path(tempfile.mkdtemp(prefix=f"ic-conf-{case['id']}-"))
    failures: list[str] = []
    try:
        materialise(case.get("files", {}), tmp)
        env = case.get("env", {})
        try:
            proc = run(case["argv"], ic, tmp, env, case["timeout_s"])
        except subprocess.TimeoutExpired:
            return {"id": case["id"], "pass": False, "why": [f"no exit within {case['timeout_s']} s"]}
        expect = case["expect"]
        check_exit(expect.get("exit"), proc.returncode, failures)
        for needle in expect.get("stderr_contains", []):
            if needle not in proc.stderr:
                failures.append(f"stderr lacks {needle!r}")
        doc = None
        if "json_file" in expect:
            doc = load_report(Path(expand(expect["json_file"], tmp)), failures)
            if doc is not None:
                check_report(doc, expect, build_commit, failures)
        same = expect.get("same_outputs_as")
        if same and doc is not None:
            try:
                other = run(same["argv"], ic, tmp, env, case["timeout_s"])
            except subprocess.TimeoutExpired:
                other = None
                failures.append(f"the comparison run did not exit within {case['timeout_s']} s")
            if other is not None:
                check_exit("zero", other.returncode, failures)
                ref = load_report(Path(expand(same["json_file"], tmp)), failures)
                if ref is not None:
                    for dotted in same["paths"]:
                        a, b = get_path(doc, dotted), get_path(ref, dotted)
                        if a is MISSING or a != b:
                            failures.append(f"{dotted} differs from the comparison run's")
        return {"id": case["id"], "step": case.get("step", "B0"), "pass": not failures, "why": failures,
                "exit": proc.returncode, "stderr_tail": proc.stderr.strip().splitlines()[-3:]}
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def cases_through(step: str) -> list[dict]:
    v1 = [{**c, "step": "B0"} for c in json.loads((V1 / "cases.json").read_text())["cases"]]
    v2 = json.loads((HERE / "cases.json").read_text())["cases"]
    last = STEPS.index(step)
    return [c for c in v1 + v2 if STEPS.index(c["step"]) <= last]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ic", required=True, type=Path)
    ap.add_argument("--through", required=True, choices=STEPS)
    ap.add_argument("--build-commit")
    ap.add_argument("--out", type=Path)
    args = ap.parse_args()
    results = [run_case(c, args.ic.resolve(), args.build_commit) for c in cases_through(args.through)]
    report = {"binary": str(args.ic), "through": args.through, "cases": len(results),
              "passed": sum(r["pass"] for r in results), "results": results}
    text = json.dumps(report, indent=1)
    if args.out:
        args.out.write_text(text + "\n")
    print(text)
    sys.exit(0 if report["passed"] == report["cases"] else 1)


if __name__ == "__main__":
    main()
