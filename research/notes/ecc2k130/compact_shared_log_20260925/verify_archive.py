#!/usr/bin/env python3
"""Offline integrity and independent replay of the committed shared-log archive."""
from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path

import verify

HERE = Path(__file__).resolve().parent
EVIDENCE = HERE / "evidence"


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def main() -> None:
    assert EVIDENCE.is_dir(), "archive not committed"
    ledger = (EVIDENCE / "SHA256SUMS").read_text().splitlines()
    listed = set()
    for row in ledger:
        digest, separator, name = row.partition("  ")
        assert separator and len(digest) == 64
        file = EVIDENCE / name
        assert file.is_file() and sha(file.read_bytes()) == digest, name
        listed.add(name)
    actual = {str(file.relative_to(EVIDENCE)) for file in EVIDENCE.rglob("*")
              if file.is_file() and file.name != "SHA256SUMS"}
    assert listed == actual
    manifest = json.loads((EVIDENCE / "archive_manifest.json").read_bytes())
    for own in ("make_inputs.py", "run.py", "run_panel.py", "verify.py", "verify_archive.py"):
        assert sha((EVIDENCE / "source" / own).read_bytes()) == sha((HERE / own).read_bytes())
    assert sha((EVIDENCE / "source/pr747_verify.py").read_bytes()) == verify.PR747_SHA
    assert sha((EVIDENCE / "source/pr737_verify.py").read_bytes()) == verify.PR737_SHA
    run_names = set()
    completed = 0
    censored = 0
    for n in manifest["panel_names"]:
        panel = json.loads((EVIDENCE / f"panel_n{n}.json").read_bytes())
        assert panel["n"] == n
        assert panel["input_spec_sha256"] == sha((HERE / "input_spec.json").read_bytes())
        training = EVIDENCE / "runs" / f"n{n}-train"
        train_report = None
        for attempt in panel["attempts"]:
            if "run_dir" not in attempt:
                assert attempt["status"] == "CENSORED_CURVE_WALL_BUDGET"
                censored += 1
                continue
            name = f"n{n}-{attempt['name']}"
            run_names.add(name)
            run = EVIDENCE / "runs" / name
            receipt = json.loads((run / "receipt.json").read_bytes())
            run_manifest = json.loads((run / "manifest.json").read_bytes())
            assert sha((run / "manifest.json").read_bytes()) == receipt["manifest_sha256"]
            assert sha((run / "producer.stderr.txt").read_bytes()) == receipt["stderr_sha256"]
            with gzip.open(run / "producer.stdout.jsonl.gz", "rb") as stream:
                stdout = stream.read()
            assert sha(stdout) == receipt["stdout_sha256"]
            assert sha(stdout) == manifest["uncompressed_stdout_sha256"][name]
            assert sha((run / run_manifest["input_file"]).read_bytes()) == run_manifest["input_sha256"]
            assert run_manifest["n"] == n and run_manifest["mode"] == attempt["mode"]
            assert run_manifest["source_sha256"] == manifest["source_sha256"][run_manifest["source_path"]]
            assert attempt["producer_receipt_sha256"] == sha((run / "receipt.json").read_bytes())
            if attempt["status"] == "COMPLETE" or attempt["status"] == "CENSORED_INCOMPLETE_LOGS":
                assert receipt["returncode"] == 0 and receipt["termination"] is None
                assert (run / "independent_validation.json").is_file()
                if attempt["mode"] == "train":
                    fresh = verify.training(run)
                    train_report = fresh
                elif attempt["mode"] == "compact":
                    assert train_report is not None
                    fresh = verify.compact(training, run, train_report)
                else:
                    fresh = verify.rho(run)
                saved = json.loads((run / "independent_validation.json").read_bytes())
                assert fresh == saved, name
                assert attempt["classification"] == saved["classification"]
                if attempt["status"] == "COMPLETE":
                    completed += 1
                else:
                    censored += 1
            else:
                assert attempt["status"].startswith("CENSORED_")
                censored += 1
    assert run_names == set(manifest["run_names"])
    assert set((EVIDENCE / "runs").iterdir()) == {EVIDENCE / "runs" / name for name in run_names}
    print(json.dumps({"classification": "ARCHIVE_REPLAY_PASS", "panels": manifest["panel_names"],
                      "runs": len(run_names), "complete_children": completed,
                      "censored_children": censored,
                      "sha256sums_sha256": sha((EVIDENCE / "SHA256SUMS").read_bytes())},
                     sort_keys=True))


if __name__ == "__main__":
    main()
