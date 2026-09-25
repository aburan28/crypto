#!/usr/bin/env python3
"""Fail-closed checksum and independent group-law replay of all archived arms."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path
import subprocess

import verify

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
EVIDENCE = HERE / "evidence"
ARMS = [(37, 3), (41, 8), (41, 12)]


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def main(evidence: Path) -> None:
    global EVIDENCE
    EVIDENCE = evidence.resolve()
    lines = (EVIDENCE / "SHA256SUMS").read_text().splitlines()
    expected = set()
    for line in lines:
        digest, sep, relative = line.partition("  ")
        assert sep and len(digest) == 64
        path = EVIDENCE / relative
        assert path.is_file() and sha(path.read_bytes()) == digest, relative
        assert relative not in expected
        expected.add(relative)
    actual = {str(path.relative_to(EVIDENCE)) for path in EVIDENCE.rglob("*")
              if path.is_file() and path.name != "SHA256SUMS"}
    assert expected == actual
    meta = json.loads((EVIDENCE / "archive_manifest.json").read_bytes())
    assert sorted(meta["run_names"]) == sorted(f"{mode}-n{n}-R{r}"
                                                    for n, r in ARMS for mode in ("oracle", "extractor"))
    for source_path, digest in meta["source_sha256"].items():
        snap = EVIDENCE / "source" / (Path(source_path).name + ".gz")
        with gzip.open(snap, "rb") as stream:
            data = stream.read()
        assert sha(data) == digest == sha((REPO / source_path).read_bytes())
    historical = subprocess.check_output(["git", "show",
        "fc27150df3238b6863ed5618c721e7fd8b6ce403:examples/koblitz_s5_sat_instance.rs"], cwd=REPO)
    assert sha(historical) == meta["source_sha256"]["examples/koblitz_s5_sat_instance.rs"]
    full_n41 = (HERE / "target_points_n41.jsonl").read_bytes()
    prefix_n41 = (HERE / "target_points_n41_R12.jsonl").read_bytes()
    assert prefix_n41 == b"".join(full_n41.splitlines(keepends=True)[:128])
    amendment = json.loads((HERE / "input_amendment.json").read_bytes())
    assert sha(full_n41) == amendment["full_n41_sha256"]
    assert sha(prefix_n41) == amendment["prefix_sha256"]
    assert len(prefix_n41.splitlines()) == amendment["prefix_count"] == 128
    expected_attempts = {"extractor-n37-R3-host-contended", "oracle-n41-R12-full-file-prefix",
                         "extractor-n41-R12-512-invalid"}
    assert set(meta["attempts"]) == expected_attempts
    for name in expected_attempts:
        run = EVIDENCE / "attempts" / name
        manifest = json.loads((run / "manifest.json").read_bytes())
        receipt = json.loads((run / "receipt.json").read_bytes())
        assert sha((run / "manifest.json").read_bytes()) == receipt["manifest_sha256"]
        with gzip.open(run / "producer.stdout.jsonl.gz", "rb") as stream:
            raw = stream.read()
        assert sha(raw) == receipt["stdout_sha256"]
        if name == "extractor-n41-R12-512-invalid":
            assert sha(raw) == amendment["failed_512_attempt_stdout_sha256"]
            assert sha((run / "manifest.json").read_bytes()) == amendment["failed_512_attempt_manifest_sha256"]
            assert sha((run / "receipt.json").read_bytes()) == amendment["failed_512_attempt_receipt_sha256"]
            assert json.loads(raw)["compact_orbit_point_batch"]["targets_requested"] == 512
        assert manifest["runner_sha256"] == meta["historical_source_sha256"]["run_before_prefix.py"]
    old_runner = (EVIDENCE / "source/run_before_prefix.py").read_bytes()
    old_protocol = (EVIDENCE / "source/PROTOCOL_before_prefix.md").read_bytes()
    assert sha(old_runner) == meta["historical_source_sha256"]["run_before_prefix.py"]
    assert sha(old_protocol) == meta["historical_source_sha256"]["PROTOCOL_before_prefix.md"]
    for filename, expected in meta["method_source_sha256"].items():
        assert sha((HERE / filename).read_bytes()) == expected
    for n, r in ARMS:
        dirs = {mode: EVIDENCE / "runs" / f"{mode}-n{n}-R{r}"
                for mode in ("oracle", "extractor")}
        for mode, run in dirs.items():
            manifest = json.loads((run / "manifest.json").read_bytes())
            receipt = json.loads((run / "receipt.json").read_bytes())
            assert (manifest["mode"], manifest["n"], manifest["R"]) == (mode, n, r)
            assert sha((run / "manifest.json").read_bytes()) == receipt["manifest_sha256"]
            assert sha((run / "producer.stderr.txt").read_bytes()) == receipt["stderr_sha256"]
            with gzip.open(run / "producer.stdout.jsonl.gz", "rb") as stream:
                stdout = stream.read()
            assert sha(stdout) == receipt["stdout_sha256"] == meta["stdout_sha256"][f"{mode}-n{n}-R{r}"]
            assert receipt["returncode"] == 0 and not receipt["timed_out"] and not receipt["sampled_rss_stop"]
            prior_run = ((n, r) == (41, 8) or ((n, r) == (37, 3) and mode == "oracle"))
            expected_protocol = old_protocol if prior_run else (HERE / "PROTOCOL.md").read_bytes()
            expected_runner = old_runner if prior_run else (HERE / "run.py").read_bytes()
            assert manifest["protocol_sha256"] == sha(expected_protocol)
            assert manifest["input_manifest_sha256"] == sha((HERE / "input_manifest.json").read_bytes())
            assert manifest["source_sha256"] == meta["source_sha256"][manifest["source_path"]]
            assert sha((REPO / manifest["base_file"]).read_bytes()) == manifest["base_file_sha256"]
            assert sha((REPO / manifest["targets_file"]).read_bytes()) == manifest["targets_sha256"]
            assert manifest["runner_sha256"] == sha(expected_runner)
            if (n, r) == (41, 12):
                assert manifest["targets_sha256"] == amendment["prefix_sha256"]
                assert manifest["input_amendment_sha256"] == sha((HERE / "input_amendment.json").read_bytes())
                assert manifest["count"] == 128
        report = verify.verify_arm(n, r,
                                  dirs["oracle"] / "producer.stdout.jsonl.gz",
                                  dirs["extractor"] / "producer.stdout.jsonl.gz")
        data = json.dumps(report, indent=2, sort_keys=True).encode() + b"\n"
        filename = f"independent-n{n}-R{r}.json"
        assert sha(data) == meta["validation_sha256"][filename]
        assert data == (EVIDENCE / "validations" / filename).read_bytes()
        print(json.dumps({"n": n, "R": r, "oracle_members": report["oracle_members"],
                          "extractor_hits": report["extractor_hits"],
                          "full_independent_miss_indices": report["full_independent_miss_indices"]},
                         sort_keys=True))
    print(json.dumps({"verdict": "PASS", "arms": len(ARMS)}, sort_keys=True))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence", type=Path, default=HERE / "evidence")
    main(parser.parse_args().evidence)
