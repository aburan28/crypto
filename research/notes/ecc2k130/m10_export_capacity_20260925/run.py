#!/usr/bin/env python3
"""Held, externally bounded two-arm representation-size experiment runner."""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(HERE.parent / "symbolic_dag_dimacs_gate_20260925"))
from bounded import run_child, sha  # noqa: E402

PARENTS = (802, 804, 784)
PARENT_BLOBS = {
    "research/notes/ecc2k130/symbolic_dag_fullpoint_20260925/dag.py": "a277fb1b2f9b58ca8cc7cb33713950015f8c0446651b1d4e8657650819bcd419",
    "research/notes/ecc2k130/symbolic_dag_dimacs_gate_20260925/export.py": "77159083a897f65b91d0f762ef425eda2df73460a91b997e754351c7c331e8fe",
    "research/notes/ecc2k130/rotated_unequal_arity_20260925/evidence/raw/m10/producer/result.json": "65bc8e10ea1c058a0dddb3a6b8cc2c4de016b3b454b3850baf72552bf557820e",
}


def git(*args: str) -> str:
    return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()


def main_blob_sha(relative: str) -> str:
    data = subprocess.check_output(["git", "show", f"origin/main:{relative}"], cwd=ROOT)
    return hashlib.sha256(data).hexdigest()


def release_gate(frozen: dict) -> dict:
    if frozen["release_main_head"] is None:
        raise RuntimeError("NOT_ADMITTED: #802, #804 and #784 must merge; rebase and re-freeze")
    pr_states = {}
    for number in PARENTS:
        item = json.loads(subprocess.check_output(
            ["gh", "pr", "view", str(number), "--json", "state,mergeCommit"],
            cwd=ROOT, text=True))
        if item["state"] != "MERGED" or not item["mergeCommit"]:
            raise RuntimeError(f"NOT_ADMITTED: prerequisite #{number} unmerged")
        pr_states[str(number)] = item["mergeCommit"]["oid"]
    subprocess.run(["git", "fetch", "origin", "main"], cwd=ROOT, check=True)
    main_head = git("rev-parse", "origin/main")
    if main_head != frozen["release_main_head"]:
        raise RuntimeError("NOT_ADMITTED: main moved after capacity re-freeze")
    for number, merge_oid in pr_states.items():
        check = subprocess.run(["git", "merge-base", "--is-ancestor", merge_oid, main_head],
                               cwd=ROOT, check=False)
        if check.returncode != 0:
            raise RuntimeError(f"NOT_ADMITTED: prerequisite #{number} is not in frozen main")
    subprocess.run(["git", "merge-base", "--is-ancestor", main_head, "HEAD"],
                   cwd=ROOT, check=True)
    for relative, expected in PARENT_BLOBS.items():
        if main_blob_sha(relative) != expected:
            raise RuntimeError(f"NOT_ADMITTED: merged parent blob drift: {relative}")
    return {"main_head": main_head, "merged_parent_commits": pr_states}


def manifest(out: Path):
    rows = []
    for path in sorted(p for p in out.rglob("*") if p.is_file() and p.name != "MANIFEST.json"):
        rows.append({"path": path.relative_to(out).as_posix(),
                     "bytes": path.stat().st_size, "sha256": sha(path)})
    (out / "MANIFEST.json").write_text(json.dumps(rows, sort_keys=True, indent=2) + "\n")
    return sha(out / "MANIFEST.json")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    receipt = {"schema": "ecc2k130-m10-export-capacity-receipt-v1",
               "status": "STOP", "scope": "representation size only; no solver or PDP",
               "attempts": []}
    try:
        subprocess.run([sys.executable, str(HERE / "ci_replay.py")], cwd=ROOT, check=True)
        frozen = json.loads((HERE / "FROZEN.json").read_text())
        spec = json.loads((HERE / "INPUT.json").read_text())
        receipt["freeze_sha256"] = sha(HERE / "FROZEN.json")
        receipt["release"] = release_gate(frozen)
        caps = spec["caps"]
        for arm in ("balanced", "unequal"):
            command = [sys.executable, str(HERE / "produce.py"), "--arm", arm,
                       "--out", str(out / arm)]
            raw = run_child(
                command, cwd=ROOT, stdout=out / f"{arm}.stdout.txt",
                stderr=out / f"{arm}.stderr.txt",
                wall_cap=caps["per_arm_external_wall_seconds"],
                rss_cap=caps["per_arm_address_space_and_rss_bytes"])
            raw["arm"] = arm
            result_path = out / arm / "result.json"
            raw["result_sha256"] = sha(result_path)
            if result_path.is_file():
                try:
                    child = json.loads(result_path.read_text())
                    raw["reported_status"] = child.get("status")
                    raw["reported_peak_rss_bytes"] = child.get("resources", {}).get("peak_rss_bytes")
                except ValueError as exc:
                    raw["result_parse_error"] = str(exc)
            receipt["attempts"].append(raw)
        if len(receipt["attempts"]) == 2 and all(
            item["exit_code"] == 0 and item["stop_reason"] is None
            and item["reported_status"] in (
                "WITHIN_FROZEN_REPRESENTATION_CAPS", "CENSORED_CNF_BYTE_CAP",
                "CENSORED_DAG_NODE_CAP")
            for item in receipt["attempts"]):
            receipt["status"] = "PRODUCER_ARCHIVED"
        else:
            receipt["status"] = "STOP"
    except BaseException as exc:
        receipt["status"] = "STOP"
        receipt["error_type"] = type(exc).__name__
        receipt["error"] = str(exc)
    finally:
        receipt["manifest_sha256"] = manifest(out)
        (out / "receipt.json").write_text(json.dumps(receipt, sort_keys=True, indent=2) + "\n")
    return 0 if receipt["status"] == "PRODUCER_ARCHIVED" else 1


if __name__ == "__main__":
    sys.exit(main())
