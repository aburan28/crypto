"""Freeze a new N53 rank-243 panel and the same full-rank stress controls."""

from __future__ import annotations

import json

from run_paired_holdout import HERE, sha


def main() -> None:
    prior_path = HERE / "rank242/FREEZE.json"
    prior = json.loads(prior_path.read_text())
    fresh_path = HERE / "rank243/holdout/n53/fixtures.json"
    fresh = json.loads(fresh_path.read_text())
    assert fresh["count"] == 24 and fresh["seed_label"] == "s3-partial-rank243-fresh-20261009-v1"
    assert prior["rank_policies"] == {"shallow": 237, "deep": 242}
    targets = []
    for index in range(1, 25):
        fixture = fresh["fixtures"][index - 1]
        assert fixture["label"] == f"T{index:03}"
        relative = f"rank243/holdout/n53/T{index:03}/public_target.json"
        targets.append({"label": f"F{index:03}", "family": "fresh",
                        "point_file_relative": relative,
                        "point_sha256": sha(HERE / relative),
                        "public_point": fixture["public_point"],
                        "fixture_scalar": int(fixture["fixture_scalar"]),
                        "scalar_derivation_sha256": fixture["scalar_derivation_sha256"]})
    stress = [target for target in prior["targets"]
              if target["family"] == "full_rank_stress_control"]
    assert [target["label"] for target in stress] == ["S067", "S085", "S096"]
    sequence = targets[:12] + stress + targets[12:]
    plan = {**prior,
            "kind": "s3_rank243_fresh_n53_panel_v1",
            "preceding_rank242_plan_sha256": sha(prior_path),
            "fresh_fixture_sha256": sha(fresh_path),
            "fresh_seed_label": fresh["seed_label"],
            "rank_policies": {"shallow": 242, "deep": 243},
            "targets": sequence}
    output = HERE / "rank243/FREEZE.json"
    if output.exists():
        raise FileExistsError(output)
    output.write_text(json.dumps(plan, sort_keys=True, separators=(",", ":")) + "\n")
    print(json.dumps({"plan": str(output), "sha256": sha(output),
                      "target_count": len(sequence)}))


if __name__ == "__main__":
    main()
