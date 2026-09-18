#!/usr/bin/env python3
"""100 then 2000 GF(2^23) trials for the poly12 selector against selected."""
from pathlib import Path
import json, re, subprocess, time
from experiment import ROOT, OUT, sha

SRC = OUT / "collision_study_poly12.cpp"
BIN = ROOT / "build/collision-poly12"


def parse(text):
    rows = {}
    for line in text.splitlines():
        match = re.search(r"mode=(\S+) .*solved=(\d+) bad=(\d+) invariant_bad=(\d+) overdue=(\d+).*iterations_mean=([0-9.]+)", line)
        if match:
            rows[match.group(1)] = {
                "solved": int(match.group(2)),
                "bad": int(match.group(3)),
                "invariant_bad": int(match.group(4)),
                "overdue": int(match.group(5)),
                "iterations_mean": float(match.group(6)),
                "line": line,
            }
    return rows


def main():
    compile_log = OUT / "poly15-collision-compile.log"
    assert not compile_log.exists()
    started = time.monotonic()
    with compile_log.open("x") as stream:
        result = subprocess.run(
            ["g++", "-O3", "-std=c++17", "-I", str(ROOT.parent), str(SRC), "-o", str(BIN)],
            cwd=ROOT.parent, stdout=stream, stderr=subprocess.STDOUT, timeout=120,
        )
    assert result.returncode == 0, compile_log.read_text()[-2000:]
    rows = {"compile_seconds": time.monotonic() - started, "compile_log_sha256": sha(compile_log)}
    for trials in (100, 2000):
        log = OUT / f"poly15-collision-{trials}.txt"
        assert not log.exists()
        started = time.monotonic()
        with log.open("x") as stream:
            run = subprocess.run(
                [str(BIN), str(trials)],
                cwd=ROOT.parent, stdout=stream, stderr=subprocess.STDOUT, timeout=3600,
            )
        text = log.read_text().strip()
        parsed = parse(text)
        selected = parsed["selected"]
        candidate = parsed["poly12-sparse-bridge3"]
        ratio = candidate["iterations_mean"] / selected["iterations_mean"]
        clean = (
            run.returncode == 0
            and selected["solved"] == trials
            and candidate["solved"] == trials
            and selected["bad"] == 0 and candidate["bad"] == 0
            and selected["invariant_bad"] == 0 and candidate["invariant_bad"] == 0
            and selected["overdue"] == 0 and candidate["overdue"] == 0
        )
        row = {
            "returncode": run.returncode,
            "elapsed_seconds": time.monotonic() - started,
            "log_sha256": sha(log),
            "output": text,
            "parsed": parsed,
            "work_ratio": ratio,
            "passed": clean and ratio <= 1.10,
        }
        rows[f"trials_{trials}"] = row
        print(text, flush=True)
        print(json.dumps({"trials": trials, "work_ratio": ratio, "passed": row["passed"]}), flush=True)
        if not row["passed"]:
            (OUT / "poly15-collision.json").write_text(json.dumps(rows, indent=2) + "\n")
            raise SystemExit(1)
        if trials == 100:
            continue
    (OUT / "poly15-collision.json").write_text(json.dumps(rows, indent=2) + "\n")
    print(json.dumps(rows, indent=2), flush=True)


if __name__ == "__main__":
    main()
