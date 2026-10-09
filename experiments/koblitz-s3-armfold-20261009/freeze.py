"""Freeze matched ARM table and PMULL-fold IC1 candidates."""

from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
import subprocess


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PARENT = ROOT / "experiments/koblitz-s3-crosscurve-20261009/FREEZE.json"
BATCHED = ROOT / "experiments/koblitz-s3-batched-target-20261009/FREEZE.json"
SOURCE = ROOT / "examples/s3_partial_prepared_batched_group_targetspan.rs"
SHARED = ROOT / "src/cryptanalysis/semaev_decomp.rs"
BINARY = Path("/private/tmp/s3-prepared-target-20261009/release/examples/s3_partial_prepared_batched_group_targetspan")


def canonical(value: dict) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False).encode()


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def git(*args: str) -> str:
    return subprocess.check_output(["git", "-C", str(ROOT), *args], text=True).strip()


def write_exact(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        assert path.read_bytes() == data, f"frozen artifact changed: {path}"
    else:
        path.write_bytes(data)


def main() -> None:
    assert subprocess.check_output(["uname", "-m"], text=True).strip() == "arm64"
    parent_bytes = PARENT.read_bytes()
    parent = json.loads(parent_bytes)
    batched_bytes = BATCHED.read_bytes()
    batched = json.loads(batched_bytes)
    assert batched["parent_plan_sha256"] == sha(parent_bytes)
    assert batched["source_sha256"] == sha(SOURCE.read_bytes())
    source_tree = git("rev-parse", "HEAD:src")
    source_commit = git("rev-parse", "HEAD")
    shared_sha = sha(SHARED.read_bytes())
    curves = {}
    for n in (41, 53):
        old_id = batched["curves"][str(n)]["batched_candidate_id"]
        old_path = BATCHED.parent / "identities/candidates" / f"{old_id}.json"
        base = json.loads(old_path.read_bytes())
        assert base.pop("candidate_id") == old_id
        assert sha(canonical(base))[:12] == old_id[-12:]
        modes = {}
        for mode, enabled in (("table", False), ("fold", True)):
            record = copy.deepcopy(base)
            record["implementation"]["rust_library_source_tree_git_oid"] = source_tree
            record["implementation"]["field_arithmetic_source_sha256"] = shared_sha
            record["implementation"]["frozen_source_commit"] = source_commit
            record["implementation"]["platform_constraint"] = "physical_aarch64_with_aes_pmull"
            record["implementation"]["algorithm_flags"]["GF2_ARM_PMULL_FOLD"] = (
                "1" if enabled else "unset")
            record["point_decomposition"]["field_reduction_kernel"] = (
                "aarch64_pmull_sparse_tail_fold" if enabled else "aarch64_pmull_table")
            ident = (f"IC1N{n}Ckb1fb{parent['curves'][str(n)]['factor_base_points']}"
                     f"PDP4rootRCsampleLAgaussTDdirectISO0h{sha(canonical(record))[:12]}")
            write_exact(HERE / "identities/candidates" / f"{ident}.json",
                        canonical({**record, "candidate_id": ident}) + b"\n")
            modes[f"{mode}_candidate_id"] = ident
        curves[str(n)] = {
            "curve_id": parent["curves"][str(n)]["curve_id"],
            "target_query_workers": batched["curves"][str(n)]["target_query_workers"],
            **modes}
    selected = [{"n": target["n"], "label": target["label"],
                 "workload_id": target["workload_id"],
                 "point_sha256": target["point_sha256"]}
                for target in parent["targets"]
                if 25 <= int(target["label"][1:]) <= 32]
    assert len(selected) == 16
    freeze = {"kind": "s3_arm_pmull_fold_paired_v1",
              "parent_plan_sha256": sha(parent_bytes),
              "batched_plan_sha256": sha(batched_bytes),
              "example_source_sha256": sha(SOURCE.read_bytes()),
              "field_arithmetic_source_sha256": shared_sha,
              "local_binary_sha256": sha(BINARY.read_bytes()),
              "source_tree_git_oid": source_tree,
              "source_commit": source_commit,
              "architecture": "aarch64",
              "kernel_modes": {"table": "GF2_ARM_PMULL_FOLD unset",
                               "fold": "GF2_ARM_PMULL_FOLD=1"},
              "relation_seed": parent["relation_seed"],
              "workers": parent["workers"],
              "prepare_rank": 243,
              "timeout_seconds_per_run": parent["timeout_seconds_per_run"],
              "repetitions_per_arm_per_target": 2,
              "curves": curves,
              "targets": selected,
              "controlled_speedup": None}
    write_exact(HERE / "FREEZE.json", canonical(freeze) + b"\n")
    print(json.dumps({"field_source_sha256": shared_sha,
                      "binary_sha256": freeze["local_binary_sha256"],
                      "curves": curves,
                      "plan_sha256": sha((HERE / "FREEZE.json").read_bytes())}, sort_keys=True))


if __name__ == "__main__":
    main()
