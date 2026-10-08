#!/usr/bin/env python3
"""Bounded static interface and inherited-corpus admission audit; never runs a solver."""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import re
import resource
import shutil
import subprocess
import sys
import tarfile
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rss_bytes() -> int:
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return value if sys.platform == "darwin" else value * 1024


def binary_inventory(names: list[str], supplements: list[dict]) -> dict:
    found = {}
    for name in names:
        executable = shutil.which(name)
        if executable is None:
            found[name] = {"available_on_PATH": False}
            continue
        path = Path(executable).resolve()
        try:
            output = subprocess.run([str(path), "--version"], text=True,
                                    capture_output=True, timeout=3, check=False)
            version = (output.stdout + output.stderr).splitlines()[:3]
            version_status = output.returncode
        except (OSError, subprocess.TimeoutExpired) as error:
            version = [repr(error)]
            version_status = "unavailable"
        found[name] = {"available_on_PATH": True, "path": str(path),
                       "bytes": path.stat().st_size, "sha256": sha(path),
                       "version_first_lines": version, "version_exit": version_status}
    for item in supplements:
        path = Path(item["path"])
        extra = {"on_PATH": False, "path": str(path), "exists_on_host": path.is_file(),
                 "expected_sha256": item["expected_sha256"]}
        if path.is_file():
            extra.update({"bytes": path.stat().st_size, "sha256": sha(path)})
            assert extra["sha256"] == item["expected_sha256"], item["name"]
        found[item["name"]] = extra
    return found


def audit(data: dict, include_binaries: bool) -> dict:
    source_hashes = {}
    for label, item in data["reference_files"].items():
        path = ROOT / item["path"]
        actual = sha(path)
        assert actual == item["sha256"], f"source or reference drift: {label}"
        source_hashes[label] = actual
    builder = (ROOT / data["reference_files"]["koblitz_builder"]["path"]).read_text()
    binary_semaev = (ROOT / data["reference_files"]["binary_semaev"]["path"]).read_text()
    example = (ROOT / data["reference_files"]["sat_example"]["path"]).read_text()
    export_note = (ROOT / data["reference_files"]["export_gate_result"]["path"]).read_text()
    theorem_note = (ROOT / data["reference_files"]["fibre_theorem_result"]["path"]).read_text()
    assert re.search(r"pub const MAX_VARS:\s*usize\s*=\s*64;", builder)
    assert "SymElement::from_subspace_vars(basis, i * ell, n, n_vars)" in builder
    assert "if n_vars > MAX_VARS" in builder
    functions = sorted(set(re.findall(r"pub fn (binary_semaev_s\d+)\s*\(", binary_semaev)))
    assert functions == ["binary_semaev_s3", "binary_semaev_s4"], functions
    assert "Four factor-base points" in example.splitlines()[2]
    assert "No direct resultant, recursive polynomial" in export_note
    assert "An exporter must encode the infinity/inverse branches" in theorem_note

    archive = ROOT / data["reference_files"]["corpus_archive"]["path"]
    facts = {}
    with tarfile.open(archive, "r:gz") as stream:
        for arm in data["arms"]:
            label = arm["name"]
            member = stream.extractfile(f"raw/{label}/targets.json")
            assert member is not None
            targets = json.load(io.TextIOWrapper(member, encoding="utf-8"))
            assert len(targets) == 8
            assert [row["class"] for row in targets] == ["planted"] * 4 + ["negative"] * 4
            assert all(sum(row["coset_multiplicities"]) > 0 for row in targets[:4])
            assert all(sum(row["coset_multiplicities"]) == 0 for row in targets[4:])
            s3 = json.loads((ROOT / data["reference_files"][
                "n13_s3_result" if label == "n13-m5" else "n19_s3_result"]["path"]).read_text())
            assert s3["arm"] == label and len(s3["targets"]) == 32
            assert s3["rational_point_tuple_count"] == arm["labelled_tuples"]
            assert s3["factor_sizes"] == arm["factor_sizes"]
            for q_index, point in enumerate(targets):
                for t_index in range(4):
                    row = s3["targets"][4 * q_index + t_index]
                    assert row["Q_index"] == q_index and row["T_index"] == t_index
                    assert row["target_class"] == point["class"]
                    assert row["true_point_tuple_count"] == point["coset_multiplicities"][t_index]
            # This is a target-level observation for the already frozen corpus,
            # not a proof that an affine-only encoding is complete in general.
            affine_per_q = [sum(s3["targets"][4*i+j]["candidate_path_count"]
                                for j in range(4)) for i in range(8)]
            assert all(v > 0 for v in affine_per_q[:4])
            assert affine_per_q[4:] == [0] * 4
            chain_vars = arm["m"] * arm["d"] + (arm["m"] - 2) * arm["n"]
            facts[label] = {"point_targets": 8, "planted_Q": 4,
                            "oracle_certified_negative_Q": 4,
                            "exact_Q_plus_T_branches": 32,
                            "point_tuple_count": arm["labelled_tuples"],
                            "affine_S3_paths_per_Q": affine_per_q,
                            "raw_chain_bits_before_auxiliaries": chain_vars,
                            "generic_u64_builder_accepts_layout": chain_vars <= 64}
    exceptional = json.loads((ROOT / data["reference_files"]["n13_s3_result"]["path"]).read_text())["targets"][12]
    assert exceptional["target"] == [7256, 3272]
    assert exceptional["true_exceptional_only_masks"] == [[0, 0, 0, 2, 1]]
    assert exceptional["true_mask_count"] == 3 and exceptional["candidate_mask_count"] == 2
    facts["n13-m5"]["frozen_O_prefix_control"] = {
        "Q_index": 3, "T_index": 0, "target": [7256, 3272],
        "exceptional_only_mask": [0, 0, 0, 2, 1],
        "full_point_masks": 3, "affine_candidate_masks": 2}
    facts["audited_interface"] = {
        "generic_boolean_monomial_limit": 64,
        "generic_builder_one_shared_basis_for_all_slots": True,
        "binary_semaev_functions": functions,
        "S5_example_factor_points": 4,
        "complete_rotated_S6_S7_or_branch_complete_chain_export_in_audited_paths": False}
    missing = {
        "n13-m5": ["frozen_rotated-slot polynomial/circuit exporter",
                     "complete O/inverse and rational-lift branches",
                     "four-coset projected model-to-point equivalence"],
        "n19-m6": ["frozen rotated-slot polynomial/circuit exporter",
                     "multiword handling for 88 chain bits before auxiliaries",
                     "complete O/inverse and rational-lift branches",
                     "four-coset projected model-to-point equivalence"]}
    result = {"domain": data["domain"], "base_commit": data["base_commit"],
              "reference_sha256": source_hashes, "facts": facts,
              "admitted_solver_arms": [], "missing_prerequisites": missing,
              "decision": "BLOCKED_BEFORE_SOLVER_TIMING"}
    if include_binaries:
        result["binary_inventory"] = binary_inventory(data["binary_names"], data["supplementary_binaries"])
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--skip-binaries", action="store_true")
    args = parser.parse_args()
    data = json.loads((HERE / "INPUT.json").read_text())
    wall, cpu = time.perf_counter(), time.process_time()
    try:
        result = audit(data, not args.skip_binaries)
        result["wall_seconds"] = time.perf_counter() - wall
        result["cpu_seconds"] = time.process_time() - cpu
        result["peak_rss_bytes"] = rss_bytes()
        assert result["wall_seconds"] <= data["caps"]["preflight_wall_seconds"]
        assert result["peak_rss_bytes"] <= data["caps"]["preflight_rss_bytes"]
        args.out.write_text(json.dumps(result, sort_keys=True, separators=(",", ":")) + "\n")
    except BaseException as error:
        args.out.with_suffix(".failure.json").write_text(json.dumps({
            "error": repr(error), "wall_seconds": time.perf_counter() - wall,
            "peak_rss_bytes": rss_bytes()}, sort_keys=True) + "\n")
        raise


if __name__ == "__main__":
    main()
