"""Bind the frozen 34-target panel to the reviewed Linux executables and runners."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PANEL = HERE / "PANEL_FREEZE.json"
GENERATED = HERE / "GENERATED_TARGETS.json"
BUILD = ROOT / "experiments/koblitz-s3-c1800-combined-20261010/modal/build/receipt.json"
OUTPUT = HERE / "EXECUTION_FREEZE.json"
EXAMPLES = {
    "baseline": "s3_partial_prepared_fullrank_c2800mixed_targetspan",
    "combined": "s3_partial_prepared_fullrank_prefetch_projective_c1800prefix_targetspan",
}


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def canonical(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False).encode()


def main() -> None:
    if OUTPUT.exists():
        raise FileExistsError(OUTPUT)
    panel_raw = PANEL.read_bytes()
    panel = json.loads(panel_raw)
    generated_raw = GENERATED.read_bytes()
    generated = json.loads(generated_raw)
    build_raw = BUILD.read_bytes()
    build = json.loads(build_raw)
    assert panel["kind"] == "s3_n53_c1800_confirmation34_panel_freeze_v1"
    assert panel["generated_targets_sha256"] == sha(generated_raw)
    assert build["kind"] == "s3_modal_c1800_combined_cross_build_v1"
    assert build["exit_code"] == 0 and build["target"] == "x86_64-unknown-linux-gnu"
    assert panel["repetitions_per_arm_per_target"] == 2
    assert len(panel["targets"]) == 34 and len(panel["run_order"]) == 136
    assert panel["relation_workers"] == 14 and panel["target_workers"] == 1
    assert panel["timeout_seconds_per_process"] == 600
    assert panel["observed_peak_rss_envelope_bytes"] == 34359738368
    generated_targets = {row["label"]: row for row in generated["targets"]}
    assert set(generated_targets) == set(panel["target_labels"])
    for target in panel["targets"]:
        assert target == generated_targets[target["label"]]
        point_raw = (ROOT / target["point_file_relative"]).read_bytes()
        assert sha(point_raw) == target["point_sha256"]
        assert json.loads(point_raw) == target["public_point"]
        workload_raw = (ROOT / target["workload_manifest_relative"]).read_bytes()
        assert sha(workload_raw) == target["workload_manifest_sha256"]
        workload = json.loads(workload_raw)
        assert workload["identity_sha256"] == target["workload_identity_sha256"]
        assert workload["workload_id"] == target["workload_id"]
        assert workload["record"]["target_point"] == target["public_point"]
    binaries = {}
    for arm, example in EXAMPLES.items():
        frozen = panel["arms"][arm]
        source_raw = (ROOT / frozen["source_relative"]).read_bytes()
        assert sha(source_raw) == frozen["source_sha256"] == build["source_sha256"][example]
        manifest_raw = (ROOT / frozen["candidate_manifest_relative"]).read_bytes()
        assert sha(manifest_raw) == frozen["candidate_manifest_sha256"]
        manifest = json.loads(manifest_raw)
        candidate_id = manifest.pop("candidate_id")
        assert candidate_id == frozen["candidate_id"]
        assert candidate_id.endswith("h" + sha(canonical(manifest))[:12])
        built = build["binaries"][example]
        binary_raw = Path(built["path"]).read_bytes()
        assert binary_raw[:5] == b"\x7fELF\x02"
        assert int.from_bytes(binary_raw[18:20], "little") == 62
        assert sha(binary_raw) == built["sha256"]
        assert len(binary_raw) == built["bytes"]
        binaries[arm] = {"example": example, "sha256": built["sha256"],
                         "bytes": built["bytes"], "architecture": "x86_64-linux-elf"}
    execution = {
        "kind": "s3_n53_c1800_confirmation34_execution_freeze_v1",
        "panel_freeze_sha256": sha(panel_raw),
        "generated_targets_sha256": sha(generated_raw),
        "linux_build_receipt_sha256": sha(build_raw),
        "linux_build_source_tree_git_oid": build["library_source_tree_git_oid"],
        "binaries": binaries,
        "single_runner_sha256": sha((HERE / "run_one.py").read_bytes()),
        "panel_runner_sha256": sha((HERE / "run_panel.py").read_bytes()),
        "analyzer_sha256": sha((HERE / "analyze_panel.py").read_bytes()),
        "sage_replay_sha256": sha((HERE / "replay_sage.py").read_bytes()),
        "run_order_sha256": sha(canonical(panel["run_order"])),
        "planned_run_count": 136,
        "memory_envelope_bytes": panel["observed_peak_rss_envelope_bytes"],
        "per_process_timeout_seconds": panel["timeout_seconds_per_process"],
        "host_exclusive_partition_receipt": None,
        "controlled_speedup": None,
    }
    OUTPUT.write_text(json.dumps(execution, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"execution_freeze_sha256": sha(OUTPUT.read_bytes()),
                      "planned_run_count": execution["planned_run_count"],
                      "binaries": binaries}, sort_keys=True))


if __name__ == "__main__":
    main()
