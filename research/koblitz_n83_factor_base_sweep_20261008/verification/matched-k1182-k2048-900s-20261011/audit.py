"""Audit the frozen matched-size N83 search receipts and same-point input."""

from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys
import tarfile


HERE = Path(__file__).resolve().parent
STUDY = HERE.parents[1]
REPO = HERE.parents[3]
sys.path.insert(0, str(STUDY))
import chain_search_supervisor as search  # noqa: E402
import v2_prefix_lp_verify as prefix  # noqa: E402


def digest(path: Path) -> str:
    h = sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read(path: Path) -> dict:
    return json.loads(path.read_text())


def container_absent(name: str) -> bool:
    result = subprocess.run(
        ["docker", "ps", "-aq", "--filter", f"name=^{name}$"],
        capture_output=True, text=True, timeout=15, check=True,
    )
    return not result.stdout.strip()


def main() -> None:
    evidence = HERE / "evidence-01"
    build = read(evidence / "build/build.json")
    assert build["schema"] == "n83.ordered-first-pair-linux-build/v1"
    assert build["status"] == "PASS_build_only" and build["exit_code"] == 0
    assert build["cleanup_verified"] and build["source_clean"] and build["source_same"]
    assert not build["timed_out"] and not build["solver_search_executed"]
    assert build["wall_cap_seconds"] == build["conservative_charge_seconds"] == 280
    assert digest(Path(build["binary"])) == build["binary_sha256"]
    assert container_absent("n83-ordered-pair-build-20261011")

    frozen = {}
    for columns in (1182, 2048):
        panel = evidence / f"panel-k{columns}"
        manifest = read(panel / "manifest.json")
        row = manifest["bases"][0]
        object_path = (Path(read(evidence / f"k{columns}/config.json")["panel_dir"])
                       / row["object"])
        prefix.manifest_row(panel / "manifest.json", columns, object_path)
        replay, storage = read(panel / "replay.json"), read(panel / "upload-receipt.json")
        assert replay["status"] == storage["status"] == "PASS"
        assert replay["checks"] == [{
            "object": row["object"], "points_checked": columns * 166,
            "representatives_checked": columns, "status": "PASS",
        }]
        assert storage["objects"] == [{
            "compressed_blake3": row["compressed_blake3"],
            "downloaded_hash_matches": True, "s3_uri": row["s3_uri"],
        }]
        frozen[columns] = (row, object_path)

    geometry = prefix.check_prefix(frozen[1182][1], frozen[2048][1], 1182, 2048, 166)
    assert geometry["representative_prefix"] is True
    assert geometry["point_row_prefix"] is True
    git_source_present = subprocess.run(
        ["git", "-C", str(REPO), "cat-file", "-e", f"{build['source_commit']}^{{commit}}"],
        capture_output=True, timeout=15,
    ).returncode == 0
    with tarfile.open(evidence / "source-snapshot.tar.gz", "r:gz") as archive:
        members = archive.getmembers()
        assert len(members) == len(search.SOURCES)
        assert {member.name for member in members} == set(search.SOURCES)
        source_bytes = {}
        for member in members:
            assert member.isfile() and 0 < member.size <= 5_000_000
            source = archive.extractfile(member)
            assert source is not None
            data = source.read()
            assert len(data) == member.size
            source_bytes[member.name] = data
    arms = []
    targets = []
    for columns in (1182, 2048):
        arm = evidence / f"k{columns}"
        config, outer = read(arm / "config.json"), read(arm / "outer.json")
        worker_config = read(arm / "run/config.json")
        worker = read(arm / "run/summary.json")
        assert search.receipt_error(worker, worker_config, config) is None
        assert config["source_commit"] == worker["source_commit"] == build["source_commit"]
        assert config["binary_sha256"] == build["binary_sha256"]
        assert config["domain_encoding"] == "frobenius_orbit_phase"
        assert config["columns"] == columns and config["summands"] == 5
        assert config["policy"] == "public_x_hash" and config["seed"] == 2026100801
        assert config["max_trials"] == 1 and config["max_variables"] == 150000
        assert config["max_domain_clauses"] == 1000000
        assert config["max_models"] == 32 and config["conflict_budget"] == 5000
        assert config["wall_seconds"] == 280 and config["cpu_limit"] == 1.0
        assert config["memory_cgroup_limit_bytes"] == 4096 * 1024 * 1024
        assert config["memory_cgroup_swap_limit_bytes"] == 0
        assert config["selected_object_sha256"] == digest(frozen[columns][1])
        assert config["selected_point_set_blake3"] == frozen[columns][0]["point_set_blake3"]
        for name, expected in config["panel_sha256"].items():
            assert digest(arm.parent / f"panel-k{columns}" / name) == expected
            assert digest(Path(config["panel_dir"]) / name) == expected
        for name, expected in config["fixture_sha256"].items():
            assert digest(Path(config["fixture_dir"]) / name) == expected
        assert set(config["source_sha256"]) == set(search.SOURCES) | {"attestation.json"}
        for name, expected in config["source_sha256"].items():
            if name == "attestation.json":
                actual = digest(arm / "source" / name)
            else:
                actual = sha256(source_bytes[name]).hexdigest()
                if git_source_present:
                    committed = sha256(subprocess.check_output(
                        ["git", "-C", str(REPO), "show", f"{build['source_commit']}:{name}"],
                        timeout=15,
                    )).hexdigest()
                    assert committed == actual
            assert actual == expected
        for name, path in (("config_sha256", arm / "config.json"),
                           ("worker_config_sha256", arm / "run/config.json"),
                           ("worker_sha256", arm / "run/summary.json"),
                           ("stdout_sha256", arm / "stdout.log"),
                           ("stderr_sha256", arm / "stderr.log")):
            assert outer[name] == digest(path)
        assert outer["status"] == worker["status"] == "UNKNOWN_solver_cap"
        assert outer["worker_exit_code"] == 0
        assert outer["group_replay"] == search.replay_group(
            worker, Path(config["fixture_dir"]) / "probe-corpus.json",
            frozen[columns][1],
        ) == {"status": "PASS", "queries": 1, "relations": 0}
        assert all(outer[name] is False for name in
                   ("binary_changed", "config_changed", "input_changed", "source_changed"))
        assert all(outer[name] is None for name in
                   ("cleanup_error", "launch_error", "parse_error", "worker_error",
                    "group_replay_error"))
        assert container_absent(config["container_name"])
        report = worker["report"]
        assert report["factor_base_size"] == columns * 166
        assert report["trials"] == report["sat_calls"] == report["sat_unknowns"] == 1
        assert report["sat_conflicts"] == 5000
        assert all(report[name] == 0 for name in
                   ("sat_models", "relations", "independent_relations", "matrix_rows",
                    "terminal_matrix_rank", "verification_failures", "sat_group_rejected_models"))
        assert worker["verified_log"] is None
        assert worker["total_index_calculus_runtime_ms"] is None
        assert worker["selected_best_total_runtime"] is None
        assert len(worker["attempt_records"]) == 1
        targets.append(worker["attempt_records"][0]["target"])
        arms.append({
            "columns": columns, "status": outer["status"],
            "selected_object_s3_uri": frozen[columns][0]["s3_uri"],
            "selected_object_sha256": config["selected_object_sha256"],
            "outer_wall_ms_L0": outer["process_wall_ms"],
            "base_import_ms_L0": worker["base_import_ms"],
            "solver_ms_L0": worker["solver_ms"],
            "cgroup_peak_bytes": outer["worker_cgroup_peak_bytes"],
            "conflicts": report["sat_conflicts"], "models": report["sat_models"],
            "relations": report["relations"], "rank": report["terminal_matrix_rank"],
            "container_cleanup_verified": True,
        })
    assert targets[0] == targets[1]
    audit = {
        "schema": "n83.matched-retained-size-audit/v1",
        "status": "PASS_bounded_receipts_only",
        "source_commit": build["source_commit"],
        "source_git_object_and_snapshot_checked": git_source_present,
        "binary_sha256": build["binary_sha256"],
        "target": targets[0], "prefix_geometry": geometry,
        "build_wall_ms_L0": build["wall_elapsed_seconds_L0"] * 1000,
        "arms": arms,
        "conservative_stage_charges_seconds": [280, 280, 280],
        "cleanup_and_receipt_reserve_charged_seconds": 60,
        "additional_charged_seconds": 900,
        "prior_charged_seconds": 10800,
        "cumulative_charged_seconds": 11700,
        "selected_best_total_runtime": None,
        "total_index_calculus_runtime_ms": None,
    }
    (evidence / "audit.json").write_text(json.dumps(audit, sort_keys=True, indent=2) + "\n")
    print(json.dumps({"status": audit["status"], "audit": str(evidence / "audit.json")}))


if __name__ == "__main__":
    main()
