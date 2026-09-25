#!/usr/bin/env python3
"""One cold, bounded producer/verifier run with exact child transcripts."""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import resource
import subprocess
import sys
import time
from pathlib import Path

from audit import HERE, canonical, preflight


def utc() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="microseconds")


def sha_file(path: Path) -> dict:
    raw = path.read_bytes()
    return {"sha256": hashlib.sha256(raw).hexdigest(), "bytes": len(raw)}


def rss_bytes(raw: int) -> int:
    return raw if sys.platform == "darwin" else raw * 1024


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--output-dir", required=True, type=Path)
    args = p.parse_args()
    start_utc, start = utc(), time.monotonic()
    frozen, inp = preflight()
    cap = inp["caps"]
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    assert not any(out.iterdir()), "output directory must be empty; never overwrite a first run"
    records = []
    status, reason = "FAIL", "not started"

    def child(label: str, argv: list[str]) -> None:
        nonlocal reason
        before = resource.getrusage(resource.RUSAGE_CHILDREN)
        child_start_utc, child_start = utc(), time.monotonic()
        remaining = cap["wall_seconds"] - (child_start - start)
        assert remaining > 0, "total wall cap exhausted before child"
        timed_out = False
        try:
            cp = subprocess.run(argv, cwd=HERE, capture_output=True, timeout=remaining)
            stdout, stderr, exit_code = cp.stdout, cp.stderr, cp.returncode
        except subprocess.TimeoutExpired as e:
            timed_out = True
            stdout, stderr, exit_code = e.stdout or b"", e.stderr or b"", -999
        child_end_utc, wall = utc(), time.monotonic() - child_start
        after = resource.getrusage(resource.RUSAGE_CHILDREN)
        stdout_path, stderr_path = out / f"{label}.stdout.txt", out / f"{label}.stderr.txt"
        stdout_path.write_bytes(stdout)
        stderr_path.write_bytes(stderr)
        record = {
            "label": label, "argv": argv, "utc_start": child_start_utc,
            "utc_end": child_end_utc, "exit_code": exit_code,
            "timed_out": timed_out, "wall_seconds": wall,
            "cpu_user_seconds": after.ru_utime - before.ru_utime,
            "cpu_system_seconds": after.ru_stime - before.ru_stime,
            "children_ru_maxrss_bytes_upper": rss_bytes(after.ru_maxrss),
            "stdout": {"path": stdout_path.name, **sha_file(stdout_path)},
            "stderr": {"path": stderr_path.name, **sha_file(stderr_path)},
        }
        records.append(record)
        if timed_out or exit_code != 0 or record["children_ru_maxrss_bytes_upper"] > cap["rss_bytes"]:
            reason = f"{label}: timeout={timed_out}, exit={exit_code}, rss_upper={record['children_ru_maxrss_bytes_upper']}"
            raise RuntimeError(reason)

    try:
        producer = out / "producer.json"
        verifier = out / "verify.json"
        child("producer", [sys.executable, str(HERE / "audit.py"), "--run", "--output", str(producer)])
        assert producer.is_file(), "producer receipt missing"
        child("independent-verifier", [sys.executable, str(HERE / "verify.py"),
                                       "--producer", str(producer), "--output", str(verifier)])
        assert verifier.is_file(), "independent verifier receipt missing"
        produced = json.loads(producer.read_text())
        checked = json.loads(verifier.read_text())
        assert produced["toy_mu4_map"] == checked["toy_mu4_map"]
        assert produced["width"]["status"] == "NOT_ADMITTED"
        status, reason = "PASS", "both child transcripts and receipts matched"
    except Exception as e:
        reason = f"{type(e).__name__}: {e}"
    end_utc, wall = utc(), time.monotonic() - start
    files = {path.name: sha_file(path) for path in sorted(out.iterdir()) if path.is_file()}
    receipt = {
        "schema": "symbolic-oaware-width-gate-run-v1", "status": status,
        "reason": reason, "frozen_sha256": sha_file(HERE / "FROZEN.json")["sha256"],
        "source_commit": frozen["parent_commit"], "utc_start": start_utc,
        "utc_end": end_utc, "wall_seconds": wall,
        "wall_cap_seconds": cap["wall_seconds"], "rss_cap_bytes": cap["rss_bytes"],
        "children": records, "files": files,
    }
    blob = canonical(receipt)
    assert len(blob) <= cap["receipt_bytes"]
    (out / "receipt.json").write_bytes(blob)
    print(json.dumps({"status": status, "receipt_sha256": hashlib.sha256(blob).hexdigest(),
                      "receipt_bytes": len(blob), "reason": reason}, sort_keys=True))
    if status != "PASS" or wall > cap["wall_seconds"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
