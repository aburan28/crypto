#!/usr/bin/env python3
"""Fail-closed replay of every committed sparse-base and same-Q receipt."""

from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
EVIDENCE = HERE / "evidence"


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
