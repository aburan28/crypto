#!/usr/bin/env python3
"""Fail-closed replay of every committed sparse-base and same-Q receipt."""

from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path
import subprocess

HERE = Path(__file__).resolve().parent
EVIDENCE = HERE / "evidence"
REPO = HERE.parents[3]


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def main() -> None:
    manifest = json.loads((EVIDENCE / "archive_manifest.json").read_bytes())
    sums = (EVIDENCE / "SHA256SUMS").read_text().splitlines()
    expected_paths = []
    for line in sums:
        expected, sep, relative = line.partition("  ")
        assert sep and len(expected) == 64
        path = EVIDENCE / relative
        assert path.is_file() and digest(path.read_bytes()) == expected, relative
        expected_paths.append(path.resolve())
    actual_paths = sorted(path.resolve() for path in EVIDENCE.rglob("*")
                          if path.is_file() and path.name != "SHA256SUMS")
    assert sorted(expected_paths) == actual_paths
    for source_name, expected in manifest["source_sha256"].items():
        with gzip.open(EVIDENCE / "source" / (source_name + ".gz"), "rb") as stream:
            assert digest(stream.read()) == expected
    assert digest((EVIDENCE / "source/pr737_independent_math.py").read_bytes()) == manifest["pr737_independent_math_sha256"]
    names = manifest["run_names"]
    assert len(names) == len(set(names))
    runs = {name: EVIDENCE / "runs" / name for name in names}
    assert sorted(runs.values()) == sorted(path for path in (EVIDENCE / "runs").iterdir() if path.is_dir())
    metadata = {name: json.loads((run / "manifest.json").read_bytes()) for name,run in runs.items()}
    # A source snapshot proves the bytes used by replay; the historical Git
    # object proves that each receipt names a real, reachable source commit.
    for name, meta in metadata.items():
        revision = f"{meta['source_commit']}:{meta['source_path']}"
        historical = subprocess.check_output(["git", "show", revision], cwd=REPO)
        assert digest(historical) == meta["source_sha256"], name
    failures = manifest["runner_failures"]
    assert len(failures) == len(set(failures))
    assert sorted(failures) == sorted(path.name for path in (EVIDENCE / "runner_failures").iterdir() if path.is_dir())
    for name in failures:
        failed = EVIDENCE / "runner_failures" / name
        if not (failed / "manifest.json").exists():
            continue  # Setup failure happened before source/producer receipt.
        meta = json.loads((failed / "manifest.json").read_bytes())
        revision = f"{meta['source_commit']}:{meta['source_path']}"
        historical = subprocess.check_output(["git", "show", revision], cwd=REPO)
        with gzip.open(failed / "producer_source.rs.gz", "rb") as stream:
            assert stream.read() == historical
        assert digest(historical) == meta["source_sha256"]
        receipt = json.loads((failed / "receipt.json").read_bytes())
        with gzip.open(failed / "producer.stdout.jsonl.gz", "rb") as stream:
            assert digest(stream.read()) == receipt["stdout_sha256"]
    final_training = {}
    for name, meta in metadata.items():
        if meta["mode"] == "train":
            key = (meta["n"], meta["R"])
            if key not in final_training or meta["count"] > metadata[final_training[key]]["count"]:
                final_training[key] = name
    import verify
    import verify_holdouts
    reports = {}
    for name, run in runs.items():
        meta = metadata[name]
        assert manifest["uncompressed_stdout_sha256"][name] == json.loads((run / "receipt.json").read_bytes())["stdout_sha256"]
        if meta["mode"] == "train":
            report = verify.verify_training(run)
        elif meta["mode"] == "holdout":
            key = (meta["n"], meta["R"])
            assert key in final_training
            report = verify_holdouts.check_holdout(runs[final_training[key]], run)
        elif meta["mode"] == "rho":
            report = verify_holdouts.check_rho(run)
        else:
            raise AssertionError(meta["mode"])
        validation = (json.dumps(report, indent=2, sort_keys=True) + "\n").encode()
        assert validation == (run / "independent_validation.json").read_bytes(), name
        reports[name] = report
    assert len([x for x in reports.values() if x["classification"] == "SPARSE_COMPACT_BASE_DIAGNOSTIC"]) >= 12
    print(json.dumps({"verdict": "PASS", "runs": len(reports),
                      "training_final": {f"n{n}_R{r}": name for (n,r),name in final_training.items()}},
                     sort_keys=True))


if __name__ == "__main__":
    main()
