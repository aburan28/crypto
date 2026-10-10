"""Freeze fresh N53 target identities and ABBA/BAAB order before a host pilot."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
STUDY = ROOT / "experiments/koblitz-s3-c1800-combined-20261010"
PARENT = STUDY / "CONFIRM_FREEZE.json"
ROBUSTNESS = STUDY / "robustness.json"
GENERATED = HERE / "GENERATED_TARGETS.json"
VERIFICATION = HERE / "local/verification.json"
OUTPUT = HERE / "TARGET_FREEZE.json"


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def canonical(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False).encode()


def main() -> None:
    if OUTPUT.exists():
        raise FileExistsError(OUTPUT)
    parent_raw = PARENT.read_bytes()
    parent = json.loads(parent_raw)
    robustness_raw = ROBUSTNESS.read_bytes()
    robustness = json.loads(robustness_raw)
    generated_raw = GENERATED.read_bytes()
    generated = json.loads(generated_raw)
    verification_raw = VERIFICATION.read_bytes()
    verification = json.loads(verification_raw)
    assert verification["all_scalar_replays_verified"] is True
    assert verification["target_count"] == 34
    assert verification["generated_targets_sha256"] == sha(generated_raw)
    assert verification["verifier_sha256"] == sha((HERE / "verify_targets_sage.py").read_bytes())
    assert generated["fixture_generator_sha256"] == sha((HERE / "make_targets_sage.py").read_bytes())
    assert robustness["sample_size_planning"]["estimated_distinct_targets_for_80_percent_power"] == 34
    assert robustness["sample_size_planning"]["process_runs_per_target"] == 4
    assert parent["kind"] == "s3_n53_c1800_combined_confirmation_freeze_v1"
    assert parent["repetitions_per_arm_per_target"] == 2
    assert parent["target_workers"] == 1
    assert parent["relation_workers"] == 14
    assert generated["parent_freeze_sha256"] == sha((ROOT / "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json").read_bytes())
    targets = generated["targets"]
    assert [row["label"] for row in targets] == [f"T{i:03}" for i in range(247, 281)]
    assert len({row["workload_id"] for row in targets}) == len(targets) == 34
    for row in targets:
        point_raw = (ROOT / row["point_file_relative"]).read_bytes()
        workload_raw = (ROOT / row["workload_manifest_relative"]).read_bytes()
        assert sha(point_raw) == row["point_sha256"]
        assert sha(workload_raw) == row["workload_manifest_sha256"]
        assert json.loads(point_raw) == row["public_point"]
        workload = json.loads(workload_raw)
        assert workload["record"]["target_point"] == row["public_point"]
        assert sha(canonical(workload["record"])) == row["workload_identity_sha256"]
        assert row["workload_id"] == row["workload_identity_sha256"][:12]
    arms = {}
    for arm in ("baseline", "combined"):
        frozen = parent["arms"][arm]
        manifest = ROOT / frozen["candidate_manifest_relative"]
        source = ROOT / frozen["source_relative"]
        assert sha(manifest.read_bytes()) == frozen["candidate_manifest_sha256"]
        assert sha(source.read_bytes()) == frozen["source_sha256"]
        assert json.loads(manifest.read_bytes())["candidate_id"] == frozen["candidate_id"]
        arms[arm] = {key: frozen[key] for key in (
            "candidate_id", "candidate_manifest_relative", "candidate_manifest_sha256",
            "source_relative", "source_sha256", "columns", "factor_base_points",
            "factor_base_digest", "prepare_rank", "root_table_slots",
            "target_search_schedule", "target_search_prefix_columns")}
    order = []
    for index, target in enumerate(targets):
        names = (("baseline", "combined", "combined", "baseline")
                 if index % 2 == 0 else
                 ("combined", "baseline", "baseline", "combined"))
        repetitions = {"baseline": 0, "combined": 0}
        for arm in names:
            repetitions[arm] += 1
            order.append({"label": target["label"], "arm": arm,
                          "repetition": repetitions[arm]})
    assert len(order) == 136
    freeze = {
        "kind": "s3_n53_c1800_combined_confirmation34_target_freeze_v1",
        "comparison": "N53 C1200 packed versus C1800 prefix-prefetch-projective S3 online DLP",
        "primary_metric": "paired_one_target_online_ms_baseline_over_combined",
        "online_ratio_gate": "2.0",
        "target_bootstrap_95pct_lower_gate": "2.0",
        "development_status": "fresh_targets_unrun",
        "parent_confirmation_freeze_sha256": sha(parent_raw),
        "prior_robustness_sha256": sha(robustness_raw),
        "generated_targets_sha256": sha(generated_raw),
        "independent_verification_sha256": sha(verification_raw),
        "fixture_generator_sha256": generated["fixture_generator_sha256"],
        "sage_runtime_info_sha256": generated["sage_runtime_info_sha256"],
        "target_count": len(targets),
        "repetitions_per_arm_per_target": 2,
        "planned_process_rows": len(order),
        "targets": [{key: row[key] for key in (
            "label", "public_point", "point_file_relative", "point_sha256",
            "workload_id", "workload_identity_sha256",
            "workload_manifest_relative", "workload_manifest_sha256")}
                    for row in targets],
        "arms": arms,
        "run_order": order,
        "relation_seed": parent["relation_seed"],
        "relation_workers": parent["relation_workers"],
        "target_workers": parent["target_workers"],
        "host_pilot_receipt": None,
        "execution_resource_cap": None,
        "host_exclusive_partition_receipt": None,
        "online_speedup": None,
        "preparer_sha256": sha(Path(__file__).read_bytes()),
    }
    raw = (json.dumps(freeze, indent=2, sort_keys=True) + "\n").encode()
    OUTPUT.write_bytes(raw)
    print(json.dumps({"target_freeze_sha256": sha(raw),
                      "target_count": len(targets),
                      "planned_process_rows": len(order)}, sort_keys=True))


if __name__ == "__main__":
    main()
