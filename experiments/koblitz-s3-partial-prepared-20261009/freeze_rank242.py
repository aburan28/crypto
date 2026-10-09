"""Freeze fresh N53 rank-242 targets and full-rank stress controls before runs."""

from __future__ import annotations

import json
from pathlib import Path

from run_paired_holdout import HERE, ROOT, sha


def main() -> None:
    parent_path = HERE / "FREEZE.json"
    parent = json.loads(parent_path.read_text())
    fresh_path = HERE / "rank242/holdout/n53/fixtures.json"
    fresh = json.loads(fresh_path.read_text())
    prior_path = HERE / "holdout/n53/fixtures.json"
    prior = json.loads(prior_path.read_text())
    assert fresh["count"] == 24 and fresh["seed_label"] == "s3-partial-rank242-fresh-20261009-v1"
    assert prior["count"] == 100
    assert sha(ROOT / "examples/s3_partial_prepared_targetspan.rs") == parent["partial_source_sha256"]
    assert sha(ROOT / "experiments/koblitz-s3-pair-query-20261007-v4/candidate.rs") == parent["reference_source_sha256"]

    def fresh_target(index: int) -> dict:
        fixture = fresh["fixtures"][index - 1]
        assert fixture["label"] == f"T{index:03}"
        relative = f"rank242/holdout/n53/T{index:03}/public_target.json"
        return {"label": f"F{index:03}", "family": "fresh",
                "point_file_relative": relative,
                "point_sha256": sha(HERE / relative),
                "public_point": fixture["public_point"],
                "fixture_scalar": int(fixture["fixture_scalar"]),
                "scalar_derivation_sha256": fixture["scalar_derivation_sha256"]}

    def stress_target(index: int) -> dict:
        fixture = prior["fixtures"][index - 1]
        assert fixture["label"] == f"T{index:03}"
        relative = f"holdout/n53/T{index:03}/public_target.json"
        return {"label": f"S{index:03}", "family": "full_rank_stress_control",
                "point_file_relative": relative,
                "point_sha256": sha(HERE / relative),
                "public_point": fixture["public_point"],
                "fixture_scalar": int(fixture["fixture_scalar"]),
                "scalar_derivation_sha256": fixture["scalar_derivation_sha256"]}

    targets = ([fresh_target(i) for i in range(1, 13)] +
               [stress_target(i) for i in (67, 85, 96)] +
               [fresh_target(i) for i in range(13, 25)])
    plan = {"kind": "s3_rank242_fresh_n53_panel_v1",
            "parent_freeze_sha256": sha(parent_path),
            "fresh_fixture_sha256": sha(fresh_path),
            "prior_fixture_sha256": sha(prior_path),
            "fresh_seed_label": fresh["seed_label"],
            "fresh_target_count": 24, "stress_target_count": 3,
            "stress_target_policy": "previous_extension_full_rank_cases_67_85_96",
            "one_target_per_process": True,
            "n": 53, "relation_seed": parent["relation_seed"],
            "workers": parent["workers"],
            "timeout_seconds_per_run": parent["timeout_seconds_per_run"],
            "rank_policies": {"shallow": 237, "deep": 242},
            "arm_order": "rotate_reference_shallow_deep_by_sequence_index_mod_3",
            "partial_source_sha256": parent["partial_source_sha256"],
            "reference_source_sha256": parent["reference_source_sha256"],
            "partial_binary_sha256": parent["partial_binary_sha256"],
            "reference_binary_sha256": parent["reference_binary_sha256"],
            "targets": targets, "controlled_speedup": None}
    output = HERE / "rank242/FREEZE.json"
    if output.exists():
        raise FileExistsError(output)
    output.write_text(json.dumps(plan, sort_keys=True, separators=(",", ":")) + "\n")
    print(json.dumps({"plan": str(output), "sha256": sha(output),
                      "target_count": len(targets)}))


if __name__ == "__main__":
    main()
