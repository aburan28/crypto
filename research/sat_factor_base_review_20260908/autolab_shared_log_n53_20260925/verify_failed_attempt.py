#!/usr/bin/env python3
"""Rehash the retained first-run failure and replay its completed raw arms."""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
ORBIT = REPO / "research/sat_factor_base_review_20260908/autolab_orbit_extract_20260924"
EXPECTED_ARCHIVE = "3f525ea221922f6168b0cabb1d8d399c4d7f5509e923e7d87be838e375c9266d"
EXPECTED_OLD_CHECKER = "e339efe5c658d62f6f0533e923ee1dc3f07841006595939f27c3901ef3c79c8a"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bundle", type=Path, required=True)
    args = parser.parse_args()
    bundle = args.bundle.resolve()
    manifest = json.loads((bundle / "archive_manifest.json").read_text())
    archive = bundle / "evidence.tar.gz"
    assert manifest["archive_sha256"] == sha(archive) == EXPECTED_ARCHIVE
    assert manifest["archive_bytes"] == archive.stat().st_size
    assert manifest["classification"] == "INVALID_L32_REPLAY"
    with tempfile.TemporaryDirectory() as temporary:
        out = Path(temporary)
        with tarfile.open(archive, "r:gz") as tar:
            for member in tar:
                assert member.isfile()
                assert member.name == "SHA256SUMS" or (
                    member.name.startswith("panel/") and ".." not in Path(member.name).parts
                )
                target = out / member.name
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_bytes(tar.extractfile(member).read())
        sums = (out / "SHA256SUMS").read_text().splitlines()
        assert len(sums) == manifest["files"]
        listed = {}
        for line in sums:
            checksum, path = line.split("  ", 1)
            assert path not in listed
            listed[path] = checksum
            assert sha(out / path) == checksum
        assert set(listed) == {
            str(path.relative_to(out)) for path in (out / "panel").rglob("*") if path.is_file()
        }
        panel_path = out / "panel/panel.json"
        assert sha(panel_path) == manifest["panel_sha256"]
        assert panel_path.read_bytes() == (bundle / "panel.json").read_bytes()
        panel = json.loads(panel_path.read_text())
        assert panel["classification"] == "INVALID_L32_REPLAY"
        old_checker = bundle / "verify_pair_original.py"
        assert sha(old_checker) == panel["source_sha256"]["pair_verifier"] == EXPECTED_OLD_CHECKER
        assert 'item["s3_calls"] >= item["partner_roots"]' in old_checker.read_text()
        assert 'item["s3_calls"] >= item["partner_roots"]' not in (HERE / "verify_pair.py").read_text()
        source_paths = {
            "ic": REPO / "examples/koblitz_s5_sat_instance.rs",
            "rho": REPO / "examples/koblitz_rho_batch_ks.rs",
            "training_driver": ORBIT / "cold_batch_rank.py",
            "independent_replay": ORBIT / "independent_replay_20260924_codex/replay.py",
            "panel_runner": HERE / "run_panel.py",
            "archive_sealer": HERE / "archive.py",
        }
        for key, path in source_paths.items():
            assert sha(path) == panel["source_sha256"][key], key
        assert panel["source_sha256"]["ic"] == "c2bc8b05087df69bef9593363e9d7c112e843ef16da122da50eb29ab22115f09"
        assert sha(HERE / "PROTOCOL.md") == panel["protocol_sha256"]
        assert sha(HERE / "TARGET_MANIFEST.json") == panel["target_manifest_sha256"]
        assert sha(HERE / "generate_targets.py") == panel["generator_sha256"]
        training = out / "panel/training"
        prior = (training / "validation.json").read_bytes()
        check = subprocess.run(
            [sys.executable, str(ORBIT / "cold_batch_rank.py"), "--replay-only",
             "--out", str(training), "--targets", "512"],
            cwd=REPO, capture_output=True, text=True, timeout=600,
        )
        assert check.returncode == 0, check.stderr[-3000:]
        assert (training / "validation.json").read_bytes() == prior
        training_result = json.loads(prior)
        assert training_result["rank"] == training_result["columns"] == 220
        assert training_result["targets_extracted"] == 512
        assert training_result["factor_base_log_solution_verified"]
        assert len(panel["steps"]) == 1
        step = panel["steps"][0]
        assert (step["count"], step["block"], step["classification"]) == (
            32, 0, "INVALID_INDEPENDENT_REPLAY"
        )
        base = out / "panel/L32_b0"
        for arm in ("ic", "rho"):
            receipt = step["arms"][arm]
            assert receipt["returncode"] == 0 and not receipt["timed_out"] and not receipt["rss_gate"]
            assert sha(base / arm / f"{arm}.stdout.jsonl") == receipt["stdout_sha256"]
            assert sha(base / arm / f"{arm}.stderr.txt") == receipt["stderr_sha256"]
        assert step["verify_ic"]["returncode"] == 1
        assert sha(base / "verify_ic/verify_ic.stderr.txt") == step["verify_ic"]["stderr_sha256"]
        assert "AssertionError" in (base / "verify_ic/verify_ic.stderr.txt").read_text()
        recovered = []
        for arm in ("ic", "rho"):
            fresh = out / f"corrected_{arm}.json"
            check = subprocess.run(
                [sys.executable, str(HERE / "verify_pair.py"), "--arm", arm,
                 "--training", str(training), "--block", "0", "--count", "32",
                 "--raw", str(base / arm), "--out", str(fresh)],
                cwd=REPO, capture_output=True, text=True, timeout=600,
            )
            assert check.returncode == 0, (arm, check.stderr[-3000:])
            recovered.append(json.loads(fresh.read_text())["recovered_scalars"])
        assert len(recovered[0]) == len(recovered[1]) == 32
        assert recovered[0] == recovered[1]
        print("First failed attempt: hashes, old checker, 512 training relations and 32 same-Q logs PASS after checker correction", flush=True)


if __name__ == "__main__":
    main()
