#!/usr/bin/env python3
"""Create the commit-bound Stage 19 execution manifest after the first clean commit."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import stat
import subprocess
from typing import Any


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
ARTIFACT = HERE / "stage-19-magma-calculator-panel-20260909"
OUTPUT = ARTIFACT / "execution-manifest.json"
SCHEMA = "koblitz_magma_calculator_stage19_execution_manifest.v1"


class PreparationError(RuntimeError):
    """The first pre-execution commit is not clean or fully bound."""


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise PreparationError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


RENDERER = load_module("stage19_renderer_for_manifest", HERE / "render_stage19_magma_calculator_panel.py")
VERIFIER = load_module("stage19_verifier_for_manifest", HERE / "verify_stage19_magma_calculator_panel.py")
CHILD = load_module("stage19_child_for_manifest", HERE / "post_stage19_magma_calculator_request.py")


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def canonical_bytes(value: Any) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True) + "\n").encode()


def git(*args: str) -> bytes:
    result = subprocess.run(
        ["git", *args], cwd=REPO, capture_output=True, check=False
    )
    if result.returncode != 0:
        raise PreparationError(
            f"git {' '.join(args)} failed: {result.stderr.decode(errors='replace').strip()}"
        )
    return result.stdout


def relevant_paths() -> list[str]:
    fixed = [
        ".github/workflows/koblitz-sota-reproduction.yml",
        "scripts/process_meter.py",
        str((HERE / "STAGE19_RESULTS.md").relative_to(REPO)),
        str((HERE / "stage-19-magma-calculator-panel-protocol.json").relative_to(REPO)),
        str((HERE / "render_stage19_magma_calculator_panel.py").relative_to(REPO)),
        str((HERE / "post_stage19_magma_calculator_request.py").relative_to(REPO)),
        str((HERE / "run_stage19_magma_calculator_panel.py").relative_to(REPO)),
        str((HERE / "verify_stage19_magma_calculator_panel.py").relative_to(REPO)),
        str((HERE / "test_stage19_magma_calculator_panel.py").relative_to(REPO)),
        str(Path(__file__).resolve().relative_to(REPO)),
        str((ARTIFACT / "plan.json").relative_to(REPO)),
        str((ARTIFACT / "prepared-summary.json").relative_to(REPO)),
    ]
    inputs = sorted(
        str(path.relative_to(REPO)) for path in (ARTIFACT / "inputs").glob("*.magma")
    )
    if len(inputs) != 10:
        raise PreparationError("execution manifest requires exactly ten prepared inputs")
    return sorted(fixed + inputs)


def blob_record(revision: str, relative: str) -> dict:
    path = REPO / relative
    info = path.lstat()
    if not stat.S_ISREG(info.st_mode) or info.st_nlink != 1:
        raise PreparationError(f"relevant path is not a single-link regular file: {relative}")
    line = git("ls-tree", revision, "--", relative).decode().strip()
    if not line or "\t" not in line:
        raise PreparationError(f"relevant path is not tracked at {revision}: {relative}")
    metadata, observed_path = line.split("\t", 1)
    fields = metadata.split()
    if len(fields) != 3 or fields[1] != "blob" or observed_path != relative:
        raise PreparationError(f"unexpected Git tree entry for {relative}")
    mode, _, blob_oid = fields
    current = path.read_bytes()
    committed = git("show", f"{revision}:{relative}")
    if current != committed:
        raise PreparationError(f"working bytes differ from {revision}:{relative}")
    return {
        "path": relative,
        "git_mode": mode,
        "git_blob_oid": blob_oid,
        "bytes": len(current),
        "sha256": sha256_bytes(current),
    }


def build_manifest() -> dict:
    if OUTPUT.exists() or OUTPUT.is_symlink():
        raise PreparationError("execution-manifest.json already exists")
    status = git("status", "--porcelain=v1", "--untracked-files=all").decode()
    if status:
        raise PreparationError("first pre-execution commit must have a completely clean checkout")
    prepared = VERIFIER.verify(ARTIFACT)
    if prepared.get("artifact_status") != "prepared_not_executed":
        raise PreparationError("prepared Stage 19 artifact is not execution-free")
    plan, _ = RENDERER.build_plan()
    if (ARTIFACT / "plan.json").read_bytes() != canonical_bytes(plan):
        raise PreparationError("prepared plan differs from deterministic rendering")
    revision = git("rev-parse", "HEAD").decode().strip()
    tree = git("rev-parse", "HEAD^{tree}").decode().strip()
    if len(revision) != 40 or len(tree) != 40:
        raise PreparationError("Git did not return full commit and tree identities")
    paths = relevant_paths()
    records = [blob_record(revision, path) for path in paths]
    return {
        "schema": SCHEMA,
        "first_preexecution_commit": revision,
        "first_preexecution_tree": tree,
        "second_preexecution_commit_rule": {
            "must_descend_from_first_commit": True,
            "only_allowed_tree_delta": [
                {
                    "status": "A",
                    "path": str(OUTPUT.relative_to(REPO)),
                }
            ],
            "execution_manifest_must_be_committed": True,
        },
        "relevant_blobs": records,
        "relevant_path_count": len(records),
        "relevant_path_list_sha256": sha256_bytes(
            canonical_bytes([record["path"] for record in records])
        ),
        "prepared_plan_sha256": sha256_bytes((ARTIFACT / "plan.json").read_bytes()),
        "prepared_summary_sha256": sha256_bytes(
            (ARTIFACT / "prepared-summary.json").read_bytes()
        ),
        "execution_requirements": {
            "checkout_clean_outside_runtime_artifacts": True,
            "current_bytes_must_match_relevant_blobs": True,
            "record_exact_second_preexecution_commit_and_tree_in_run_and_attempts": True,
            "network_execution_before_second_commit": False,
        },
        "claim_boundary": plan["claim_boundary"],
    }


def self_test() -> dict:
    paths = relevant_paths()
    required_suffixes = {
        "plan.json",
        "prepared-summary.json",
        "render_stage19_magma_calculator_panel.py",
        "post_stage19_magma_calculator_request.py",
        "run_stage19_magma_calculator_panel.py",
        "verify_stage19_magma_calculator_panel.py",
        "test_stage19_magma_calculator_panel.py",
        "koblitz-sota-reproduction.yml",
    }
    names = {Path(path).name for path in paths}
    if not required_suffixes.issubset(names) or len([p for p in paths if p.endswith(".magma")]) != 10:
        raise AssertionError("execution-manifest relevant path set is incomplete")
    return {"self_test": "pass", "relevant_paths": len(paths), "prepared_inputs": 10}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        print(json.dumps(self_test(), indent=2, sort_keys=True))
        return
    if not args.write:
        parser.error("manifest creation requires --write after the first clean commit")
    try:
        manifest = build_manifest()
        CHILD.atomic_write(OUTPUT, canonical_bytes(manifest))
    except (PreparationError, VERIFIER.VerificationError, RENDERER.RenderError, CHILD.ChildError) as error:
        parser.error(str(error))
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
