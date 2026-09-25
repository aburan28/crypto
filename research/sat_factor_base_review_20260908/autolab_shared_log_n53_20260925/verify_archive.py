#!/usr/bin/env python3
"""Rehash and independently replay committed n53 shared-log panel evidence."""
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


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bundle", type=Path, required=True)
    args = parser.parse_args()
    bundle = args.bundle.resolve()
    manifest = json.loads((bundle / "archive_manifest.json").read_text())
    archive = bundle / "evidence.tar.gz"
    assert sha(archive) == manifest["archive_sha256"]
    assert archive.stat().st_size == manifest["archive_bytes"]
    with tempfile.TemporaryDirectory() as temporary:
        out = Path(temporary)
        with tarfile.open(archive, "r:gz") as tar:
            for member in tar:
                assert member.isfile()
                assert member.name == "SHA256SUMS" or (
                    member.name.startswith("panel/") and ".." not in Path(member.name).parts
                )
                destination = out / member.name
                destination.parent.mkdir(parents=True, exist_ok=True)
                destination.write_bytes(tar.extractfile(member).read())
        entries = (out / "SHA256SUMS").read_text().splitlines()
        assert len(entries) == manifest["files"]
        seen = set()
        for line in entries:
            checksum, name = line.split("  ", 1)
            assert name not in seen
            seen.add(name)
            assert sha(out / name) == checksum
        assert seen == {str(path.relative_to(out)) for path in (out / "panel").rglob("*") if path.is_file()}
        panel_path = out / "panel/panel.json"
        if not panel_path.exists():
            assert manifest["classification"] == "NO_PANEL_SUMMARY"
            raise AssertionError("archived producer infrastructure failure: no panel summary")
        assert sha(panel_path) == manifest["panel_sha256"]
        assert (bundle / "panel.json").read_bytes() == panel_path.read_bytes()
        panel = json.loads(panel_path.read_text())
        assert panel["classification"] == manifest["classification"]
        assert (panel["classification"] == "COMPLETE_ALL_SIX_PAIRS"
                or panel["classification"].startswith("CENSORED_")), panel["classification"]
        expected_sources = {
            "ic": REPO / "examples/koblitz_s5_sat_instance.rs",
            "rho": REPO / "examples/koblitz_rho_batch_ks.rs",
            "training_driver": ORBIT / "cold_batch_rank.py",
            "independent_replay": ORBIT / "independent_replay_20260924_codex/replay.py",
            "pair_verifier": HERE / "verify_pair.py",
            "panel_runner": HERE / "run_panel.py",
            "archive_sealer": HERE / "archive.py",
        }
        for key, path in expected_sources.items():
            assert sha(path) == panel["source_sha256"][key], key
        assert sha(HERE / "TARGET_MANIFEST.json") == panel["target_manifest_sha256"]
        assert sha(HERE / "generate_targets.py") == panel["generator_sha256"]
        # Protocol text can gain the result after a run without changing the
        # frozen target/producer source. Its measured hash remains in panel.json.
        training = out / "panel/training"
        if (training / "validation.json").exists():
            prior = (training / "validation.json").read_bytes()
            command = [sys.executable, str(ORBIT / "cold_batch_rank.py"),
                       "--replay-only", "--out", str(training), "--targets", "512"]
            result = subprocess.run(command, cwd=REPO, capture_output=True, text=True,
                                    timeout=600)
            assert result.returncode == 0, result.stderr[-3000:]
            assert (training / "validation.json").read_bytes() == prior
            print("training: fresh 512-relation rank/log replay PASS", flush=True)
        else:
            assert panel["classification"].startswith("CENSORED_TRAINING")
            print("training: censored, raw retained", flush=True)
            return
        for step in panel["steps"]:
            name = f"L{step['count']}_b{step['block']}"
            base = out / "panel" / name
            for arm in ("ic", "rho"):
                if arm not in step.get("arms", {}):
                    continue
                raw = base / arm / f"{arm}.stdout.jsonl"
                stderr = base / arm / f"{arm}.stderr.txt"
                receipt = step["arms"][arm]
                assert sha(raw) == receipt["stdout_sha256"]
                assert sha(stderr) == receipt["stderr_sha256"]
            if step.get("classification") != "PASS":
                assert panel["classification"].startswith("CENSORED")
                print(f"{name}: {step.get('classification')}, raw retained", flush=True)
                continue
            for arm in ("ic", "rho"):
                saved = json.loads((base / arm / "validation.json").read_text())
                fresh = out / f"fresh_{name}_{arm}.json"
                command = [sys.executable, str(HERE / "verify_pair.py"), "--arm", arm,
                           "--training", str(training), "--block", str(step["block"]),
                           "--count", str(step["count"]), "--raw", str(base / arm),
                           "--out", str(fresh)]
                result = subprocess.run(command, cwd=REPO, capture_output=True, text=True,
                                        timeout=600)
                assert result.returncode == 0, (name, arm, result.stderr[-3000:])
                assert json.loads(fresh.read_text()) == saved
            print(f"{name}: fresh same-Q IC/rho replay PASS", flush=True)
        if panel["classification"] == "COMPLETE_ALL_SIX_PAIRS":
            assert len(panel["steps"]) == 6
            assert all(step["classification"] == "PASS" for step in panel["steps"])
        print("Committed n53 panel archive PASS", flush=True)


if __name__ == "__main__":
    main()
