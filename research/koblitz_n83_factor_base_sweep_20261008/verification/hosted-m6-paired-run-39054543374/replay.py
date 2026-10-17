#!/usr/bin/env python3
"""Independently replay the frozen hosted K600 4+2 Gray/hash artifact."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
STUDY = HERE.parents[1]
ROOT = HERE.parents[3]
PROTOCOL = STUDY / "verification/hosted-m6-paired-protocol-20261017"
PRIOR = STUDY / "verification/k600-gray-hash-paired-protocol-20261014"
sys.path[:0] = [str(PROTOCOL), str(PRIOR)]
import plan  # noqa: E402
import receipt_audit as prior_audit  # noqa: E402

SOURCE_COMMIT = "17f03c518377468eeb49331ecfd8c3dedc2d72d0"
PLAN_SHA256 = "7788514127328a78cf2b0f634003ddae3ff18a9ff0c0c8cd75b50fcc006152cd"
IMAGE = "rust:1.98.1-bookworm@sha256:93ce27a88655056a51dbdd8f5f2d7ddc071c7b0070fb288a37b5a285fc83971e"
ARCHIVE_DIR = "n83-archives"
SETUP_DIR = "n83-m6-setup"
RUN_DIR = "n83-m6-pair"
PILOT_MANIFEST = STUDY / "pilot-01/manifest.json"
PUBLIC_CORPUS = STUDY / "pilot-01/probe-corpus.json"


def require(ok: bool, reason: str) -> None:
    if not ok:
        raise ValueError(reason)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_json(path: Path) -> dict:
    row = json.loads(path.read_text())
    require(isinstance(row, dict), f"JSON object missing: {path}")
    return row


def same(row: dict, expected: dict, label: str) -> None:
    for key, value in expected.items():
        require(type(row.get(key)) is type(value) and row[key] == value,
                f"{label} {key} differs")


def audit_mounts(command: list, inspected: dict, arm: dict, stage: str) -> None:
    require(isinstance(command, list) and command[:2] == ["docker", "run"]
            and command.count("--mount") == 6,
            f"{stage} Docker invocation differs")
    mounted = {}
    for index, token in enumerate(command[:-1]):
        if token != "--mount":
            continue
        parts = command[index + 1].split(",")
        fields = dict(part.split("=", 1) for part in parts if "=" in part)
        fields["readonly"] = "readonly" in parts
        require(fields.get("type") == "bind" and "dst" in fields
                and "src" in fields and fields["dst"] not in mounted,
                f"{stage} bind syntax differs")
        mounted[fields["dst"]] = fields
    expected = {
        "/worker": ("/n83-m6-setup/worker", True),
        "/base.jsonl.gz": (f"/n83-archives/{arm['base_object_blake3']}.jsonl.gz", True),
        "/pool.jsonl.gz": (f"/n83-archives/{arm['pool_object_blake3']}.jsonl.gz", True),
        "/pilot-manifest.json": (
            "/research/koblitz_n83_factor_base_sweep_20261008/pilot-01/manifest.json", True),
        "/corpus.json": (
            "/research/koblitz_n83_factor_base_sweep_20261008/pilot-01/probe-corpus.json", True),
        "/out": (f"/n83-m6-pair/{arm['role']}", False),
    }
    require(set(mounted) == set(expected), f"{stage} bind targets differ")
    actual_mounts = {mount["Destination"]: mount for mount in inspected["Mounts"]
                     if mount["Type"] == "bind"}
    require(set(actual_mounts) == set(expected), f"{stage} inspected binds differ")
    for target, (suffix, readonly) in expected.items():
        row = mounted[target]
        require(row["src"].endswith(suffix), f"{stage} bind source differs: {target}")
        require(row["readonly"] is readonly,
                f"{stage} bind mode differs: {target}")
        require(actual_mounts[target]["Source"] == row["src"]
                and actual_mounts[target]["RW"] is not readonly,
                f"{stage} inspected bind differs: {target}")


def audit_s3_bases(directory: Path, expected: dict) -> None:
    import blake3

    receipt = read_json(HERE / "factor-bases-s3-check.json")
    same(receipt, {"schema": "n83.hosted-m6-s3-base-replay/v1",
                   "status": "PASS_3_FROZEN_S3_OBJECTS",
                   "source_plan_sha256": PLAN_SHA256,
                   "selected_best_total_runtime": None}, "S3 base replay")
    rows = receipt.get("objects")
    require(isinstance(rows, list) and len(rows) == 3,
            "S3 base replay object count differs")
    by_digest = {row.get("compressed_blake3"): row for row in rows}
    require(len(by_digest) == 3, "S3 base replay object identity differs")
    for arm in expected["arms"]:
        for role in ("base", "pool"):
            digest = arm[f"{role}_object_blake3"]
            row = by_digest[digest]
            path = directory / f"{digest}.jsonl.gz"
            data = path.read_bytes()
            same(row, {"s3_uri": arm[f"{role}_s3_uri"],
                       "compressed_blake3": digest,
                       "compressed_sha256": arm[f"{role}_compressed_sha256"],
                       "bytes": len(data),
                       "etag": '"' + hashlib.md5(data).hexdigest() + '"'},
                 "S3 base")
            require(hashlib.sha256(data).hexdigest() == row["compressed_sha256"]
                    and blake3.blake3(data).hexdigest() == digest,
                    "fresh S3 base bytes differ")


def audit_search(config: dict, result: dict, arm: dict,
                 base_points: list, residual: set) -> str:
    import blake3

    same(config, {
        "schema": "n83.cross-stream-pilot-partial-config/v1",
        "study": "koblitz_n83_factor_base_sweep_20261008",
        "base_object_blake3": arm["base_object_blake3"],
        "pool_object_blake3": arm["pool_object_blake3"],
        "pilot_manifest_blake3": blake3.blake3(PILOT_MANIFEST.read_bytes()).hexdigest(),
        "source_blake3": blake3.blake3(
            (ROOT / "examples/koblitz_n83_mixed_orbit_worker.rs").read_bytes()).hexdigest(),
        "public_corpus_path": "/corpus.json",
        "public_corpus_canonical_json_blake3": prior_audit.PRIOR.CORPUS_BLAKE3,
        "budget_seconds": 2300, "conflict_budget": 20000,
        "max_models": 64, "base_summands": 4, "residual_summands": 2,
        "order_first_pair": True, "max_variables": 150000,
        "max_domain_clauses": 12000000,
        "memory_cgroup_limit_bytes": 4294967296,
        "memory_cgroup_swap_limit_bytes": 0,
    }, "search config")
    same(result, {
        "schema": "n83.cross-stream-pilot-partial-result/v1",
        "study": "koblitz_n83_factor_base_sweep_20261008",
        "base_object_blake3": arm["base_object_blake3"],
        "pool_object_blake3": arm["pool_object_blake3"],
        "base_policy": "public_x_sequential", "base_seed": 2026100802,
        "pool_policy": "public_x_gray_prefix" if arm["role"] == "gray_pool" else "public_x_hash",
        "pool_seed": 2026100801,
        "public_corpus_canonical_json_blake3": prior_audit.PRIOR.CORPUS_BLAKE3,
        "base_signed_orbits": 600, "pool_signed_orbits": 600,
        "shared_signed_orbits": arm["shared_signed_orbits"],
        "residual_signed_orbits": arm["residual_signed_orbits"],
        "trial": 0, "a": "1833078656608906180176350",
        "b": "163089659398503559150141",
        "graph_rows": None, "matrix_rank": None,
        "total_index_calculus_runtime_ms": None,
        "selected_best_total_runtime": None,
    }, "search result")
    require(prior_audit.PRIOR.finite_nonnegative(result.get("base_import_elapsed_ms_L0"))
            and prior_audit.PRIOR.finite_nonnegative(result.get("worker_wall_ms_L0"))
            and result["base_import_elapsed_ms_L0"] <= result["worker_wall_ms_L0"] <= 2300000,
            "worker time counters differ")
    peak = result.get("memory_cgroup_peak_bytes")
    require(type(peak) is int and 0 < peak <= 4294967296, "memory peak differs")
    stats = result.get("search")
    require(isinstance(stats, dict), "search stats missing")
    capacity = arm["source_capacity"]
    same(stats, {"domain_clauses": capacity["domain_clauses"]}, "search stats")
    for key, cap in (("conflicts", 20000), ("models", 64),
                     ("invalid_models", 64), ("group_rejected_models", 64),
                     ("solver_calls", 80), ("patterns_attempted", 16),
                     ("patterns_refuted", 16),
                     ("peak_variables", capacity["total_variables"]),
                     ("peak_xors", capacity["xor_rows"])):
        value = stats.get(key)
        require(type(value) is int and 0 <= value <= cap, f"search {key} exceeds cap")
    require(type(stats.get("peak_clauses")) is int and stats["peak_clauses"] >= 0,
            "peak clauses missing")
    require(type(stats.get("exhausted")) is bool
            and type(stats.get("refuted")) is bool, "solver disposition missing")
    status = result.get("status")
    require(status in {"PASS_one_bounded_partial", "PASS_one_trial_refuted",
                       "UNKNOWN_search_cap", "UNKNOWN_preflight_cap",
                       "PRODUCER_FAILURE_invalid_model"}, "search status unclassified")
    candidate = result.get("candidate")
    require((candidate is not None) == (status == "PASS_one_bounded_partial"),
            "candidate/status differs")
    require((stats["invalid_models"] > 0) == (status == "PRODUCER_FAILURE_invalid_model"),
            "invalid-model disposition differs")
    if status == "UNKNOWN_search_cap":
        require(stats["exhausted"], "search cap lacks exhausted counter")
    elif status == "UNKNOWN_preflight_cap":
        require(stats["exhausted"] and stats["patterns_attempted"] == 0
                and stats["solver_calls"] == 0, "preflight cap made solver calls")
    elif status == "PASS_one_trial_refuted":
        require(stats["refuted"] and not stats["exhausted"], "refutation counters differ")
    if candidate is not None:
        require(isinstance(candidate, dict) and candidate.get("source_id") == 1
                and (candidate.get("a"), candidate.get("b")) ==
                (result["a"], result["b"]), "candidate source differs")
        indices = candidate.get("base_point_indices")
        points = [prior_audit.PRIOR.parse_point(p)
                  for p in candidate.get("residual_points", [])]
        require(isinstance(indices, list) and len(indices) == 4
                and all(type(i) is int and 0 <= i < len(base_points) for i in indices)
                and len(points) == 2 and all(p in residual for p in points),
                "candidate lies outside frozen domains")
        require(prior_audit.PRIOR.replay_equation(
            candidate["a"], candidate["b"],
            [base_points[i] for i in indices], points, prior_audit.PUBLIC_TARGET),
            "candidate fails independent group equation")
    return status


def audit_artifact(root: Path, s3_bases: Path) -> dict:
    require(root.is_dir(), "artifact directory missing")
    expected = plan.build()
    require(sha256(plan.SAVED) == PLAN_SHA256
            and json.loads(plan.SAVED.read_text()) == expected, "frozen plan differs")
    require(tuple(map(int, expected["public_target"])) == prior_audit.PUBLIC_TARGET
            and expected["first_trial_a"] == "1833078656608906180176350"
            and expected["first_trial_b"] == "163089659398503559150141",
            "public target or first-trial coefficients differ")
    audit_s3_bases(s3_bases, expected)
    require((root / SETUP_DIR / "source-commit.txt").read_text().strip() == SOURCE_COMMIT,
            "artifact source commit differs")
    require((root / SETUP_DIR / "source-status.txt").read_text() == "",
            "hosted checkout was dirty")
    require((root / SETUP_DIR / "plan-check.json").is_file(), "plan replay missing")
    same(read_json(root / SETUP_DIR / "plan-check.json"),
         {"status": "PASS", "sha256": PLAN_SHA256}, "plan replay")
    sources = (root / SETUP_DIR / "sources.sha256").read_text().splitlines()
    require(len(sources) == 6, "hosted source list differs")
    for line in sources:
        digest, relative = line.split("  ", 1)
        require(not Path(relative).is_absolute() and sha256(ROOT / relative) == digest,
                f"hosted source hash differs: {relative}")
    image_identity = (root / SETUP_DIR / "image-identity.txt").read_text().strip()
    image_id, architecture, system, repo_digests = image_identity.split(" ", 3)
    require(image_id.startswith("sha256:") and len(image_id) == 71
            and architecture == "arm64" and system == "linux"
            and any(value.endswith("@" + IMAGE.split("@", 1)[1])
                    for value in json.loads(repo_digests)),
            "pinned Linux ARM image identity differs")
    regeneration = read_json(root / SETUP_DIR / "regeneration.json")
    same(regeneration, {"schema": "n83.hosted-regenerated-archives/v1",
                        "status": "PASS_3_FROZEN_OBJECTS",
                        "source_commit": SOURCE_COMMIT,
                        "workflow_dispatch_head_sha": SOURCE_COMMIT,
                        "pilot_completed_base_count": 54}, "regeneration")
    produced_manifest = root / "n83-pilot-parent/panel/manifest.json"
    same(read_json(produced_manifest), {
        "schema": "n83.factor-base-panel/v1",
        "status": "completed_factor_base_panel",
        "source_commit": SOURCE_COMMIT,
        "completed_base_count": 54,
        "expected_base_count": 54,
    }, "regenerated pilot manifest")
    require(sha256(produced_manifest) == regeneration["pilot_manifest_sha256"],
            "pilot manifest digest differs")
    archive_hashes = {}
    for arm in expected["arms"]:
        for role in ("base", "pool"):
            key = arm[f"{role}_object_blake3"]
            archive = root / ARCHIVE_DIR / f"{key}.jsonl.gz"
            digest = sha256(archive)
            require(digest == arm[f"{role}_compressed_sha256"],
                    f"archive {key} differs")
            archive_hashes[f"{key}.jsonl.gz"] = digest
    require(regeneration.get("archive_sha256") == archive_hashes,
            "regenerated archive receipt differs")
    config = read_json(root / RUN_DIR / "config.json")
    same(config, {"schema": "n83.hosted-m6-paired-config/v1",
                  "source_commit": SOURCE_COMMIT,
                  "plan_sha256": PLAN_SHA256,
                  "worker_sha256": sha256(root / SETUP_DIR / "worker"),
                  "manifest_sha256": sha256(PILOT_MANIFEST),
                  "corpus_sha256": sha256(PUBLIC_CORPUS),
                  "image": IMAGE,
                  "run_id": "39054543374",
                  "conservative_full_charge_seconds": 7500,
                  "selected_best_total_runtime": None}, "outer config")
    require((root / SETUP_DIR / "worker.sha256").read_text().startswith(
        config["worker_sha256"] + "  "), "worker binary hash receipt differs")
    outcomes = []
    for arm in expected["arms"]:
        role = arm["role"]
        arm_dir = root / RUN_DIR / role
        imported = read_json(arm_dir / "import.stage.json")
        same(imported, {"schema": "n83.hosted-m6-paired-stage/v1",
                        "stage": "import", "role": role,
                        "wall_cap_seconds": 300, "exit_code": 0,
                        "timed_out": False, "inspect_exit_code": 0,
                        "cleanup_exit_code": 0}, "import stage")
        require(imported["command"][-len(arm["worker_verify_argv"]):] ==
                arm["worker_verify_argv"], "import worker argv differs")
        for stage in ("import", "search"):
            receipt = read_json(arm_dir / f"{stage}.stage.json")
            require(type(receipt.get("elapsed_wall_seconds_L0")) in (int, float)
                    and 0 <= receipt["elapsed_wall_seconds_L0"]
                    <= receipt["wall_cap_seconds"] + 30,
                    f"{stage} wall counter differs")
            require(sha256(arm_dir / f"{stage}.stdout") == receipt["stdout_sha256"]
                    and sha256(arm_dir / f"{stage}.stderr") == receipt["stderr_sha256"],
                    f"{stage} output digest differs")
            inspected = json.loads((arm_dir / f"{stage}.inspect.json").read_text())
            require(isinstance(inspected, list) and len(inspected) == 1,
                    f"{stage} container inspection missing")
            state, host = inspected[0]["State"], inspected[0]["HostConfig"]
            audit_mounts(receipt["command"], inspected[0], arm, stage)
            require(state["OOMKilled"] is False and host["Memory"] == 4294967296
                    and host["MemorySwap"] == 4294967296
                    and host["NanoCpus"] == 1000000000
                    and host["PidsLimit"] == 64
                    and host["NetworkMode"] == "none"
                    and host["ReadonlyRootfs"] is True
                    and [cap.upper() for cap in host["CapDrop"]] == ["ALL"]
                    and any(option.startswith("no-new-privileges")
                            for option in host["SecurityOpt"]),
                    f"{stage} container limits differ")
        prior_audit.audit_verify(read_json(arm_dir / "import.stdout"), arm)
        base = root / ARCHIVE_DIR / f"{arm['base_object_blake3']}.jsonl.gz"
        pool = root / ARCHIVE_DIR / f"{arm['pool_object_blake3']}.jsonl.gz"
        base_points, residual = prior_audit.domains(base, pool, arm["shared_signed_orbits"])
        searched = read_json(arm_dir / "search.stage.json")
        same(searched, {"schema": "n83.hosted-m6-paired-stage/v1",
                        "stage": "search", "role": role,
                        "wall_cap_seconds": 2400, "inspect_exit_code": 0,
                        "cleanup_exit_code": 0}, "search stage")
        require(type(searched.get("timed_out")) is bool,
                "search wall-cap flag missing")
        require(searched["command"][-len(arm["worker_search_argv"]):] ==
                arm["worker_search_argv"], "search worker argv differs")
        if searched["timed_out"]:
            require(searched["exit_code"] is None,
                    "outer wall cap has process exit code")
            status = "UNKNOWN_outer_wall_cap"
        elif searched["exit_code"] == 124:
            cap = read_json(arm_dir / "search/cap.json")
            same(cap, {"schema": "n83.cross-stream-pilot-partial-cap/v1",
                       "status": "UNKNOWN_budget", "budget_seconds": 2300}, "worker cap")
            status = "UNKNOWN_worker_wall_cap"
        else:
            require(searched["exit_code"] == 0, "search container failed")
            status = audit_search(read_json(arm_dir / "search/config.json"),
                                  read_json(arm_dir / "search/result.json"),
                                  arm, base_points, residual)
        arm_outer = read_json(arm_dir / "outer.json")
        require(arm_outer["status"] == status and arm_outer["role"] == role,
                "arm outer status differs")
        outcomes.append({"role": role, "status": status,
                         "search_exit_code": searched["exit_code"],
                         "elapsed_search_wall_seconds_L0": searched["elapsed_wall_seconds_L0"]})
    outer = read_json(root / RUN_DIR / "outer.json")
    same(outer, {"schema": "n83.hosted-m6-paired-outer/v1",
                 "status": "BOUNDED_PARTIAL_one_query_pair",
                 "source_commit": SOURCE_COMMIT, "run_id": "39054543374",
                 "conservative_full_charge_seconds": 7500,
                 "natural_partial_yield": None, "matrix_rank": None,
                 "total_index_calculus_runtime_ms": None,
                 "selected_best_total_runtime": None}, "outer")
    require([r["role"] for r in outer["arms"]] == expected["arm_order"]
            and [r["status"] for r in outer["arms"]] == [r["status"] for r in outcomes],
            "outer arms differ")
    return {"schema": "n83.hosted-m6-paired-audit/v1",
            "status": "PASS_MATCHED_4_PLUS_2_ONE_QUERY_REPLAY",
            "source_commit": SOURCE_COMMIT, "plan_sha256": PLAN_SHA256,
            "archives_sha256": archive_hashes, "outcomes": outcomes,
            "s3_factor_bases_status": "PASS_3_FROZEN_S3_OBJECTS",
            "conservative_full_charge_seconds": 7500,
            "selected_best_total_runtime": None}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", required=True, type=Path)
    parser.add_argument("--s3-bases", required=True, type=Path)
    args = parser.parse_args()
    print(json.dumps(audit_artifact(args.artifact, args.s3_bases), sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
