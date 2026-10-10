"""Combine the frozen N41 screen and independent confirmation target panels."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import random
import statistics


HERE = Path(__file__).resolve().parent
PILOT = HERE.parent / "koblitz-s3-c1800-n41-20261010"


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def geometric(values: list[float]) -> float:
    assert values and all(value > 0 for value in values)
    return math.exp(statistics.mean(map(math.log, values)))


def main() -> None:
    output = HERE / "local/combined_analysis.json"
    if output.exists():
        raise FileExistsError(output)
    first_analysis_path = PILOT / "local/analysis.json"
    second_analysis_path = HERE / "local/analysis.json"
    first_freeze_path = PILOT / "FREEZE.json"
    second_freeze_path = HERE / "FREEZE.json"
    first = json.loads(first_analysis_path.read_bytes())
    second = json.loads(second_analysis_path.read_bytes())
    first_freeze = json.loads(first_freeze_path.read_bytes())
    second_freeze = json.loads(second_freeze_path.read_bytes())
    assert first["complete"] and second["complete"]
    assert first["verified_runs"] == 16 and second["verified_runs"] == 48
    assert first["freeze_sha256"] == sha(first_freeze_path.read_bytes())
    assert second["freeze_sha256"] == sha(second_freeze_path.read_bytes())
    assert second_freeze["pilot_freeze_sha256"] == sha(first_freeze_path.read_bytes())
    assert second_freeze["pilot_analysis_sha256"] == sha(first_analysis_path.read_bytes())
    assert first_freeze["source_tree_git_oid"] == second_freeze["source_tree_git_oid"]
    assert {arm: entry["candidate_id"] for arm, entry in first_freeze["arms"].items()} == {
        arm: entry["candidate_id"] for arm, entry in second_freeze["arms"].items()}
    first_labels = {entry["label"] for entry in first["target_results"]}
    second_labels = {entry["label"] for entry in second["target_results"]}
    assert len(first_labels) == 8 and len(second_labels) == 24
    assert not first_labels & second_labels
    ratios = [entry["online_ratio"] for entry in first["target_results"]]
    ratios += [entry["online_ratio"] for entry in second["target_results"]]
    rng = random.Random(20261010)
    boot = sorted(geometric([ratios[rng.randrange(len(ratios))]
                             for _ in ratios]) for _ in range(20000))
    result = {
        "kind": "s3_n41_c1200_c1800_combined32_analysis_v1",
        "pilot_analysis_sha256": sha(first_analysis_path.read_bytes()),
        "confirmation_analysis_sha256": sha(second_analysis_path.read_bytes()),
        "pilot_freeze_sha256": sha(first_freeze_path.read_bytes()),
        "confirmation_freeze_sha256": sha(second_freeze_path.read_bytes()),
        "script_sha256": sha(Path(__file__).read_bytes()),
        "pilot_target_count": 8,
        "confirmation_target_count": 24,
        "combined_target_count": 32,
        "online_geometric_baseline_over_expanded": geometric(ratios),
        "online_target_bootstrap_95pct": [boot[499], boot[19499]],
        "log_ratio_sample_sd": statistics.stdev([math.log(value) for value in ratios]),
        "combined_lower_endpoint_at_least_twofold": boot[499] >= 2.0,
        "confirmation_only_lower_endpoint_at_least_twofold":
            second["confirmation_gate_bootstrap_lower_at_least_twofold"],
        "host_exclusive_partition_receipt": None,
        "controlled_speedup": None,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"target_count": 32,
                      "combined_ratio": result["online_geometric_baseline_over_expanded"],
                      "combined_interval": result["online_target_bootstrap_95pct"],
                      "confirmation_gate": result["confirmation_only_lower_endpoint_at_least_twofold"]}))


if __name__ == "__main__":
    main()
