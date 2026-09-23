#!/usr/bin/env python3
"""Run the recorded confirmation only after the complete primary gates pass."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    primary = HERE / "run_01"
    manifest = json.loads((primary / "manifest.json").read_text())
    for name, digest in manifest["files"].items():
        assert sha(primary / name) == digest, name
    metadata = json.loads((primary / "metadata.json").read_text())
    assert metadata["complete"]
    result = json.loads((primary / "results.json").read_text())
    assert result["all_results_completed_and_verified"]
    assert any(v["retained_frontier"] == "PASS" for v in result["new_gates"].values()), "No primary treatment passed; retain the rejection instead."
    for name, digest in metadata["source_hashes"].items():
        assert sha(HERE / name) == digest, f"Primary source changed: {name}"
    protocol = json.loads((primary / "protocol.json").read_text())
    plan = json.loads((HERE / "confirmation_plan.json").read_text())
    assert plan["primary_protocol_sha256"] == sha(primary / "protocol.json")
    for key in ["variants", "repetitions", "limits", "new_candidates", "reference_arms", "matched_controls", "final_gate", "checksum_contract"]:
        assert plan[key] == protocol[key], key
    seen = set(protocol["discovery_seeds"] + protocol["regression_seeds"] + protocol["holdout_seeds"])
    assert not seen.intersection(plan["holdout_seeds"])
    assert plan["regression_seeds"] == protocol["holdout_seeds"]
    assert not (HERE / "run_02").exists(), "Use an additive successor; never overwrite a run."
    target = HERE / "confirmation_protocol.json"
    assert not target.exists(), "Confirmation protocol already exists; inspect its state."
    plan["primary_manifest_sha256"] = sha(primary / "manifest.json")
    plan["confirmation_plan_sha256"] = sha(HERE / "confirmation_plan.json")
    target.write_text(json.dumps(plan, indent=2) + "\n")
    subprocess.run([sys.executable, str(HERE / "run.py"), "--protocol", str(target), "--out", str(HERE / "run_02")], check=True)

if __name__ == "__main__":
    main()
