"""Replay the bounded N83 pair-arm receipt contract and frozen input hashes."""

from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys


HERE = Path(__file__).resolve().parent
STUDY = HERE.parents[1]
REPO = HERE.parents[3]
sys.path.insert(0, str(STUDY))
import chain_search_supervisor as search  # noqa: E402


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def read(path: Path) -> dict:
    return json.loads(path.read_text())


def absent_container(name: str) -> bool:
    result = subprocess.run(
        ["docker", "ps", "-aq", "--filter", f"name=^{name}$"],
        capture_output=True, text=True, timeout=15, check=True,
    )
    return not result.stdout.strip()


def main() -> None:
    build = read(HERE / "evidence-01/build/build.json")
    binary = Path(build["binary"])
    assert build["schema"] == "n83.ordered-first-pair-linux-build/v1"
    assert build["status"] == "PASS_build_only" and build["exit_code"] == 0
    assert build["cleanup_verified"] and build["source_clean"] and build["source_same"]
    assert not build["solver_search_executed"] and not build["timed_out"]
    assert build["conservative_charge_seconds"] == 280
    assert digest(binary) == build["binary_sha256"]
    assert absent_container("n83-ordered-pair-build-20261011")
    panel_receipts = HERE / "evidence-01/panel"
    panel = read(panel_receipts / "manifest.json")
    panel_replay = read(panel_receipts / "replay.json")
    storage = read(panel_receipts / "upload-receipt.json")
    s3_head = read(panel_receipts / "s3-head.json")
    assert panel_replay["status"] == storage["status"] == "PASS"
    assert panel_replay["checks"][0]["points_checked"] == 339968
    assert panel_replay["checks"][0]["representatives_checked"] == 2048

    rows = []
    targets = []
    for arm, encoding in (
        ("ordered", "frobenius_orbit_phase_ordered_first_pair"),
        ("unrestricted", "frobenius_orbit_phase"),
    ):
        root = HERE / "evidence-01" / arm
        config, outer = read(root / "config.json"), read(root / "outer.json")
        worker_config = read(root / "run/config.json")
        worker = read(root / "run/summary.json")
        assert search.receipt_error(worker, worker_config, config) is None
        assert config["domain_encoding"] == encoding
        assert config["source_commit"] == build["source_commit"] == worker["source_commit"]
        assert config["binary_sha256"] == build["binary_sha256"]
        assert config["selected_object_sha256"] == digest(
            Path(config["panel_dir"]) / config["selected_object"]
        )
        assert all(
            digest(panel_receipts / name) == expected == digest(Path(config["panel_dir"]) / name)
            for name, expected in config["panel_sha256"].items()
        )
        assert all(
            digest(Path(config["fixture_dir"]) / name) == expected
            for name, expected in config["fixture_sha256"].items()
        )
        for name, expected in config["source_sha256"].items():
            if name == "attestation.json":
                actual = digest(root / "source" / name)
            else:
                actual = sha256(subprocess.check_output(
                    ["git", "-C", str(REPO), "show", f"{build['source_commit']}:{name}"],
                    timeout=15,
                )).hexdigest()
            assert actual == expected
        assert outer["config_sha256"] == digest(root / "config.json")
        assert outer["worker_config_sha256"] == digest(root / "run/config.json")
        assert outer["worker_sha256"] == digest(root / "run/summary.json")
        assert outer["stdout_sha256"] == digest(root / "stdout.log")
        assert outer["stderr_sha256"] == digest(root / "stderr.log")
        assert outer["status"] == worker["status"] == "UNKNOWN_solver_cap"
        replay = search.replay_group(
            worker,
            Path(config["fixture_dir"]) / "probe-corpus.json",
            Path(config["panel_dir"]) / config["selected_object"],
        )
        assert replay == outer["group_replay"] == {"status": "PASS", "queries": 1, "relations": 0}
        assert all(
            outer[key] is False for key in
            ("binary_changed", "config_changed", "input_changed", "source_changed")
        )
        assert all(
            outer[key] is None for key in
            ("cleanup_error", "launch_error", "parse_error", "worker_error", "group_replay_error")
        )
        assert outer["worker_exit_code"] == 0
        assert outer["memory_cgroup_limit_bytes"] == 4096 * 1024 * 1024
        assert outer["memory_cgroup_swap_limit_bytes"] == 0
        assert absent_container(config["container_name"])
        report = worker["report"]
        assert report["trials"] == report["sat_calls"] == report["sat_unknowns"] == 1
        assert report["sat_conflicts"] == 1000
        assert all(
            report[key] == 0 for key in
            ("sat_models", "relations", "independent_relations", "matrix_rows", "terminal_matrix_rank")
        )
        assert worker["verified_log"] is None
        assert worker["total_index_calculus_runtime_ms"] is None
        assert worker["selected_best_total_runtime"] is None
        assert len(worker["attempt_records"]) == 1
        targets.append(worker["attempt_records"][0]["target"])
        rows.append({
            "arm": arm, "status": outer["status"], "domain_encoding": encoding,
            "outer_wall_ms_L0": outer["process_wall_ms"],
            "worker_solver_ms_L0": worker["solver_ms"],
            "base_import_ms_L0": worker["base_import_ms"],
            "cgroup_peak_bytes": outer["worker_cgroup_peak_bytes"],
            "conflicts": report["sat_conflicts"], "models": report["sat_models"],
            "relations": report["relations"], "rank": report["terminal_matrix_rank"],
            "container_cleanup_verified": True,
            "selected_object_sha256": config["selected_object_sha256"],
        })
    assert targets[0] == targets[1]
    assert rows[0]["selected_object_sha256"] == rows[1]["selected_object_sha256"]
    s3_uri = panel["bases"][0]["s3_uri"]
    assert s3_uri.startswith("s3://crypto-autoresearcher/factor-bases/icv1/etc/")
    assert s3_head["ContentLength"] == panel["bases"][0]["bytes"]
    assert s3_head["Metadata"]["blake3"] == panel["bases"][0]["compressed_blake3"]
    assert storage["objects"] == [{
        "compressed_blake3": panel["bases"][0]["compressed_blake3"],
        "downloaded_hash_matches": True, "s3_uri": s3_uri,
    }]
    audit = {
        "schema": "n83.ordered-first-pair-audit/v1", "status": "PASS_bounded_receipts_only",
        "source_commit": build["source_commit"], "binary_sha256": build["binary_sha256"],
        "selected_object_s3_uri": s3_uri, "target": targets[0],
        "build_wall_ms_L0": build["wall_elapsed_seconds_L0"] * 1000,
        "arms": rows, "conservative_stage_charges_seconds": [280, 280, 280],
        "cleanup_and_receipt_reserve_charged_seconds": 60,
        "additional_charged_seconds": 900, "prior_charged_seconds": 9900,
        "cumulative_charged_seconds": 10800,
        "solver_search_executed": True, "selected_best_total_runtime": None,
        "total_index_calculus_runtime_ms": None,
    }
    (HERE / "evidence-01/audit.json").write_text(json.dumps(audit, sort_keys=True, indent=2) + "\n")
    print(json.dumps({"status": audit["status"], "audit": str(HERE / "evidence-01/audit.json")}))


if __name__ == "__main__":
    main()
