#!/usr/bin/env python3
"""Adversarial synthetic tests for compact Stage-23 terminal evidence."""

from __future__ import annotations

from contextlib import ExitStack
from copy import deepcopy
import importlib.util
import json
import os
from pathlib import Path
import shutil
import stat
import subprocess
import tempfile
import types
import unittest
from unittest import mock

import koblitz_stage23_terminal_evidence as core

HERE = Path(__file__).resolve().parent
FIXTURE_SPEC = importlib.util.spec_from_file_location(
    "stage23_fixture_helpers", HERE / "test_koblitz_unknown_scalar_panel.py"
)
assert FIXTURE_SPEC and FIXTURE_SPEC.loader
fixtures = importlib.util.module_from_spec(FIXTURE_SPEC)
FIXTURE_SPEC.loader.exec_module(fixtures)


REPO = HERE.parent
PROTOCOL_RELATIVE = (
    "research/sat_factor_base_review_20260908/continuation-05-sota-gates/"
    "stage-23-unknown-scalar-protocol.json"
)


def write(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(data)


def run(command: list[str], cwd: Path) -> str:
    completed = subprocess.run(command, cwd=cwd, text=True, capture_output=True, check=False)
    if completed.returncode:
        raise AssertionError(f"command failed: {command!r}\n{completed.stderr}")
    return completed.stdout.strip()


def build_source_repo(root: Path) -> Path:
    source = root / "source"
    source.mkdir()
    write(
        source / ".gitignore",
        b"Cargo.lock\n!research/sat_factor_base_review_20260908/continuation-05-sota-gates/stage-20-rust-build/Cargo.lock\n",
    )
    actual_sources = {
        "scripts/run_koblitz_unknown_scalar_panel.py": HERE / "run_koblitz_unknown_scalar_panel.py",
        "scripts/run_koblitz_relation_yield_bridge.py": HERE / "run_koblitz_relation_yield_bridge.py",
        "scripts/process_meter.py": HERE / "process_meter.py",
    }
    for relative in core.DIRECT_SOURCE_PATHS:
        if relative == "Cargo.lock":
            continue
        destination = source / relative
        if relative in actual_sources:
            write(destination, actual_sources[relative].read_bytes())
        elif relative == core.FROZEN_LOCK_RELATIVE:
            write(destination, b"# synthetic frozen lock\nversion = 3\n")
        elif relative == "Cargo.toml":
            write(destination, b"[package]\nname='stage23-synthetic'\nversion='0.0.0'\n")
        else:
            write(destination, f"// synthetic closure member {relative}\n".encode())
    write(source / "Cargo.lock", (source / core.FROZEN_LOCK_RELATIVE).read_bytes())
    write(source / PROTOCOL_RELATIVE, (REPO / PROTOCOL_RELATIVE).read_bytes())
    run(["git", "init", "-q"], source)
    run(["git", "config", "user.name", "Stage23 Test"], source)
    run(["git", "config", "user.email", "stage23@example.invalid"], source)
    run(["git", "add", "-A"], source)
    run(["git", "commit", "-q", "-m", "synthetic source closure"], source)
    if run(["git", "status", "--porcelain=v1", "--untracked-files=all"], source):
        raise AssertionError("synthetic source repository is dirty")
    return source.resolve()


def patched_stage23(source: Path) -> ExitStack:
    stage23 = fixtures.stage23
    stack = ExitStack()
    replacements = {
        "REPO": source,
        "PROTOCOL": source / PROTOCOL_RELATIVE,
        "METER": source / "scripts/process_meter.py",
        "DISCOVERY_SOURCE": source / "examples/koblitz_public_factor_base_discovery.rs",
        "PANEL_SOURCE": source / "examples/koblitz_unknown_scalar_panel.rs",
        "LOCK": source / core.FROZEN_LOCK_RELATIVE,
        "RUNNER_SOURCE": source / "scripts/run_koblitz_unknown_scalar_panel.py",
    }
    for name, value in replacements.items():
        stack.enter_context(mock.patch.object(stage23, name, value))
    stack.enter_context(mock.patch.object(stage23.custody, "REPO", source))
    return stack


def make_bundle(
    root: Path, *, hardlink_build_output: bool = False
) -> tuple[Path, dict[str, Path]]:
    source = build_source_repo(root)
    evidence = root / "evidence"
    with patched_stage23(source):
        run_root, initial_outer = fixtures.make_complete_rank_deficient_run(evidence)
        outer_root = evidence / "outer"
        outer_root.mkdir()
        outer_metrics = outer_root / "driver.metrics.json"
        shutil.move(initial_outer, outer_metrics)
        write(outer_root / "driver.stdout", b"synthetic driver stdout\n")
        write(outer_root / "driver.stderr", b"")
        project = evidence / "project-verification"
        fixtures.stage23.verify(
            types.SimpleNamespace(
                run_root=run_root,
                outer_metrics=outer_metrics,
                output=project,
            )
        )
    if hardlink_build_output:
        os.link(
            run_root / "build-target/release/examples/koblitz_unknown_scalar_panel",
            root / "external-cargo-build-hardlink",
        )
    bundle = root / "bundle"
    core.package_bundle(
        run_root=run_root,
        outer_root=outer_root,
        outer_metrics=outer_metrics,
        project_verification_root=project,
        output=bundle,
    )
    return bundle, {
        "run": run_root,
        "outer": outer_root,
        "outer_metrics": outer_metrics,
        "project": project,
        "source": source,
    }


def write_json(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def reseal_bundle(bundle: Path) -> None:
    seal_path = bundle / "bundle-seal.json"
    seal = json.loads(seal_path.read_text())
    seal["manifest"] = core.identity(
        bundle / "bundle-manifest.json", recorded_path="bundle-manifest.json"
    )
    inventory = core.tree_inventory(bundle, excluded={"bundle-seal.json"})
    seal["inventory"] = inventory
    seal["inventory_sha256"] = core.canonical_sha256(inventory)
    seal.pop("seal_payload_sha256", None)
    seal["seal_payload_sha256"] = core.canonical_sha256(seal)
    write_json(seal_path, seal)


def coherently_reseal_original_layers(bundle: Path) -> None:
    manifest_path = bundle / "bundle-manifest.json"
    manifest = json.loads(manifest_path.read_text())
    roots = {entry["role"]: entry for entry in manifest["path_map"]["roots"]}
    original_run = Path(roots["run"]["original_absolute"])
    original_project = Path(roots["project_verification"]["original_absolute"])
    run_root = bundle / "original/run"
    summary_path = run_root / "run-summary.json"
    run_seal_path = run_root / "run-seal.json"
    run_seal = json.loads(run_seal_path.read_text())
    omitted = [
        record for record in run_seal["inventory"]
        if core.is_build_target_path(record["path"])
    ]
    run_seal["summary"] = core.identity(
        summary_path, recorded_path=str(original_run / "run-summary.json")
    )
    run_inventory = sorted(
        [*core.tree_inventory(run_root, excluded={"run-seal.json"}), *omitted],
        key=lambda record: record["path"],
    )
    run_seal["inventory"] = run_inventory
    run_seal["inventory_sha256"] = core.canonical_sha256(run_inventory)
    run_seal.pop("seal_payload_sha256", None)
    run_seal["seal_payload_sha256"] = core.canonical_sha256(run_seal)
    write_json(run_seal_path, run_seal)

    manifest["run_partition"]["source_inventory_sha256"] = run_seal["inventory_sha256"]
    manifest["run_partition"]["retained"] = [
        record for record in run_inventory if not core.is_build_target_path(record["path"])
    ]
    manifest["run_partition"]["omitted"] = [
        record for record in run_inventory if core.is_build_target_path(record["path"])
    ]
    manifest["run_partition"]["run_seal"] = core.identity(
        run_seal_path, recorded_path="original/run/run-seal.json"
    )

    project_root = bundle / "original/project-verification"
    verification_path = project_root / "verification.json"
    verification = json.loads(verification_path.read_text())
    verification["run_seal"] = core.identity(
        run_seal_path, recorded_path=str(original_run / "run-seal.json")
    )
    verification["run_inventory_sha256"] = run_seal["inventory_sha256"]
    write_json(verification_path, verification)
    verification_seal_path = project_root / "verification-seal.json"
    verification_seal = json.loads(verification_seal_path.read_text())
    verification_seal["verification"] = core.identity(
        verification_path, recorded_path=str(original_project / "verification.json")
    )
    verification_inventory = core.tree_inventory(
        project_root, excluded={"verification-seal.json"}
    )
    verification_seal["inventory"] = verification_inventory
    verification_seal["inventory_sha256"] = core.canonical_sha256(verification_inventory)
    verification_seal.pop("seal_payload_sha256", None)
    verification_seal["seal_payload_sha256"] = core.canonical_sha256(verification_seal)
    write_json(verification_seal_path, verification_seal)
    manifest["project_verification"]["inventory"] = core.tree_inventory(project_root)
    write_json(manifest_path, manifest)

    bundle_seal = json.loads((bundle / "bundle-seal.json").read_text())
    bundle_seal["source_run_inventory_sha256"] = run_seal["inventory_sha256"]
    write_json(bundle / "bundle-seal.json", bundle_seal)
    reseal_bundle(bundle)


def coherently_forge_summary_ratio(bundle: Path) -> None:
    summary_path = bundle / "original/run/run-summary.json"
    summary = json.loads(summary_path.read_text())
    summary["ratios"]["online_ic_over_rho_core"] = 999.0
    write_json(summary_path, summary)
    coherently_reseal_original_layers(bundle)


class Stage23TerminalEvidenceTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory(prefix="stage23-terminal-evidence-")
        self.root = Path(self.temporary.name).resolve()
        self.bundle, self.original = make_bundle(self.root)

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def clone_bundle(self, name: str) -> Path:
        destination = self.root / name
        shutil.copytree(self.bundle, destination)
        return destination

    def test_package_preserves_bytes_and_reconstructs_complete_graph(self) -> None:
        result = core.verify_bundle(self.bundle)
        self.assertEqual(result["status"], "compact_terminal_evidence_structurally_verified")
        self.assertEqual(result["task_count"], 8)
        self.assertEqual(result["row_count"], 2)
        self.assertEqual(result["completed_rows"], 2)
        self.assertEqual(result["binary_build_cross_bindings"], 2)
        self.assertTrue(result["exact_task_graph_commands_resources_rows_and_ratios_reconstructed"])
        self.assertFalse(result["scientific_measurement_admitted"])
        self.assertFalse(result["external_portable_verification_satisfied"])
        self.assertFalse(result["full_cost_gate_passed"])
        self.assertFalse(result["koblitz_index_calculus_sota"])

        preserved = (
            (self.original["run"] / "tasks/04-row-01-ic/stdout", self.bundle / "original/run/tasks/04-row-01-ic/stdout"),
            (self.original["run"] / "inputs/source.json", self.bundle / "original/run/inputs/source.json"),
            (self.original["run"] / "binaries/koblitz_unknown_scalar_panel", self.bundle / "original/run/binaries/koblitz_unknown_scalar_panel"),
            (self.original["outer"] / "driver.stdout", self.bundle / "original/outer/driver.stdout"),
            (self.original["project"] / "verification.json", self.bundle / "original/project-verification/verification.json"),
        )
        for original, archived in preserved:
            self.assertEqual(original.read_bytes(), archived.read_bytes())
        self.assertFalse((self.bundle / "original/run/build-target").exists())
        manifest = json.loads((self.bundle / "bundle-manifest.json").read_text())
        workspace_lock = manifest["source_closure"]["workspace_cargo_lock"]
        self.assertFalse(workspace_lock["tracked_at_bound_commit"])
        self.assertTrue(workspace_lock["ignored_workspace_file"])
        self.assertEqual(
            (self.bundle / "source-extra/Cargo.lock").read_bytes(),
            (self.bundle / "source-tree" / core.FROZEN_LOCK_RELATIVE).read_bytes(),
        )

    def test_packager_allows_hardlinks_only_in_omitted_cargo_build_tree(self) -> None:
        with tempfile.TemporaryDirectory(prefix="stage23-cargo-hardlink-") as temporary:
            bundle, _ = make_bundle(
                Path(temporary).resolve(), hardlink_build_output=True
            )
            result = core.verify_bundle(bundle)
            self.assertEqual(
                result["status"], "compact_terminal_evidence_structurally_verified"
            )
            self.assertEqual(result["completed_rows"], 2)

    def test_rejects_signed_extra_file_and_path_traversal(self) -> None:
        extra = self.clone_bundle("extra")
        write(extra / "signed-but-unexpected.txt", b"extra")
        reseal_bundle(extra)
        with self.assertRaisesRegex(core.EvidenceError, "missing or extra"):
            core.verify_bundle(extra)

        traversal = self.clone_bundle("traversal")
        manifest_path = traversal / "bundle-manifest.json"
        manifest = json.loads(manifest_path.read_text())
        manifest["run_partition"]["retained"][0]["path"] = "../escape"
        write_json(manifest_path, manifest)
        reseal_bundle(traversal)
        with self.assertRaisesRegex(core.EvidenceError, "archive-relative|normalized"):
            core.verify_bundle(traversal)

    def test_rejects_duplicate_json_keys_and_links(self) -> None:
        duplicate = self.clone_bundle("duplicate")
        manifest_path = duplicate / "bundle-manifest.json"
        raw = manifest_path.read_bytes()
        self.assertTrue(raw.startswith(b"{\n"))
        manifest_path.write_bytes(b'{\n  "schema": "duplicate",\n' + raw[2:])
        reseal_bundle(duplicate)
        with self.assertRaisesRegex(core.EvidenceError, "duplicate JSON key"):
            core.verify_bundle(duplicate)

        linked = self.clone_bundle("linked")
        target = linked / "original/outer/driver.stdout"
        target.unlink()
        target.symlink_to("driver.stderr")
        with self.assertRaisesRegex(core.EvidenceError, "symlink"):
            core.verify_bundle(linked)

        hardlinked = self.clone_bundle("hardlinked")
        target = hardlinked / "original/outer/driver.stdout"
        extra = hardlinked / "original/outer/hardlink"
        os.link(target, extra)
        with self.assertRaisesRegex(core.EvidenceError, "inadmissible|hard-linked"):
            core.verify_bundle(hardlinked)

    def test_rejects_cross_binding_and_coherently_resealed_ratio_forgery(self) -> None:
        cross_binding = self.clone_bundle("cross-binding")
        manifest_path = cross_binding / "bundle-manifest.json"
        manifest = json.loads(manifest_path.read_text())
        manifest["binary_build_cross_bindings"][0]["omitted_build_path"] = (
            "build-target/release/examples/other"
        )
        write_json(manifest_path, manifest)
        reseal_bundle(cross_binding)
        with self.assertRaisesRegex(core.EvidenceError, "cross-bindings"):
            core.verify_bundle(cross_binding)

        forged = self.clone_bundle("ratio-forgery")
        coherently_forge_summary_ratio(forged)
        with self.assertRaisesRegex(core.EvidenceError, "resources or ratios"):
            core.verify_bundle(forged)

    def test_rejects_coherently_resealed_extra_files_in_each_original_tree(self) -> None:
        run_extra = self.clone_bundle("run-extra")
        write(run_extra / "original/run/signed-extra.json", b"{}\n")
        coherently_reseal_original_layers(run_extra)
        with self.assertRaisesRegex(core.EvidenceError, "exact file grammar"):
            core.verify_bundle(run_extra)

        outer_extra = self.clone_bundle("outer-extra")
        write(outer_extra / "original/outer/signed-extra", b"extra")
        manifest_path = outer_extra / "bundle-manifest.json"
        manifest = json.loads(manifest_path.read_text())
        manifest["outer"]["inventory"] = core.tree_inventory(outer_extra / "original/outer")
        write_json(manifest_path, manifest)
        reseal_bundle(outer_extra)
        with self.assertRaisesRegex(core.EvidenceError, "outer tree violates"):
            core.verify_bundle(outer_extra)

        project_extra = self.clone_bundle("project-extra")
        write(project_extra / "original/project-verification/signed-extra.json", b"{}\n")
        coherently_reseal_original_layers(project_extra)
        with self.assertRaisesRegex(core.EvidenceError, "two-file grammar"):
            core.verify_bundle(project_extra)

    def test_rejects_reserved_build_target_file_and_profile_or_blocker_forgery(self) -> None:
        reserved = self.clone_bundle("reserved-build-target")
        seal_path = reserved / "original/run/run-seal.json"
        seal = json.loads(seal_path.read_text())
        seal["inventory"].append(
            {"path": "build-target", "bytes": 1, "sha256": core.sha256(b"x")}
        )
        seal["inventory"] = sorted(seal["inventory"], key=lambda record: record["path"])
        write_json(seal_path, seal)
        coherently_reseal_original_layers(reserved)
        with self.assertRaisesRegex(core.EvidenceError, "reserved build-target"):
            core.verify_bundle(reserved)

        profile = self.clone_bundle("profile-forgery")
        manifest_path = profile / "bundle-manifest.json"
        manifest = json.loads(manifest_path.read_text())
        manifest["profile"] = "production"
        write_json(manifest_path, manifest)
        reseal_bundle(profile)
        with self.assertRaisesRegex(core.EvidenceError, "manifest profile"):
            core.verify_bundle(profile)

        blocker = self.clone_bundle("blocker-forgery")
        verification_path = blocker / "original/project-verification/verification.json"
        verification = json.loads(verification_path.read_text())
        verification["external_portable_verification_blocker"] = "forged blocker"
        write_json(verification_path, verification)
        coherently_reseal_original_layers(blocker)
        with self.assertRaisesRegex(core.EvidenceError, "project verification cannot be reconstructed"):
            core.verify_bundle(blocker)

    def test_external_tool_map_is_bound_and_bundled_code_is_never_invoked(self) -> None:
        forged_tool = self.clone_bundle("forged-tool")
        manifest_path = forged_tool / "bundle-manifest.json"
        manifest = json.loads(manifest_path.read_text())
        python = next(
            entry for entry in manifest["path_map"]["external_tools"]
            if entry["role"] == "python"
        )
        python["original_absolute"] = "/usr/bin/forged-python"
        python["identity"] = {
            "path": "/usr/bin/forged-python", "bytes": 1, "sha256": "0" * 64,
        }
        write_json(manifest_path, manifest)
        reseal_bundle(forged_tool)
        with self.assertRaisesRegex(core.EvidenceError, "external tool map differs"):
            core.verify_bundle(forged_tool)

        inert = self.clone_bundle("inert-bundled-code")
        sentinel = self.root / "bundled-code-was-executed"
        bundled_verifier = inert / "verifier/verify_koblitz_stage23_terminal_evidence.py"
        bundled_verifier.write_text(
            "from pathlib import Path\n"
            f"Path({str(sentinel)!r}).write_text('executed')\n"
            "raise RuntimeError('bundled evidence code executed')\n"
        )
        manifest_path = inert / "bundle-manifest.json"
        manifest = json.loads(manifest_path.read_text())
        manifest["portable_tools"]["verifier"] = core.identity(
            bundled_verifier,
            recorded_path="verifier/verify_koblitz_stage23_terminal_evidence.py",
        )
        write_json(manifest_path, manifest)
        reseal_bundle(inert)
        result = core.verify_bundle(inert)
        self.assertFalse(result["bundled_verifier_sources_executed"])
        self.assertFalse(sentinel.exists())

    def test_descriptor_walk_rejects_symlinked_parent_component(self) -> None:
        alias = self.root / "bundle-parent-alias"
        alias.symlink_to(self.bundle.parent, target_is_directory=True)
        with self.assertRaisesRegex(core.EvidenceError, "cannot open|Too many levels"):
            core.verify_bundle(alias / self.bundle.name)

    def test_descriptor_read_rejects_same_inode_same_size_in_place_mutation(self) -> None:
        target = self.root / "same-size-race.bin"
        target.write_bytes(b"before")
        original = target.stat()
        real_fstat = os.fstat
        triggered = False

        def mutate_after_first_file_fstat(descriptor: int):
            nonlocal triggered
            metadata = real_fstat(descriptor)
            if (
                not triggered
                and stat.S_ISREG(metadata.st_mode)
                and metadata.st_dev == original.st_dev
                and metadata.st_ino == original.st_ino
            ):
                triggered = True
                writer = os.open(target, os.O_WRONLY)
                try:
                    os.write(writer, b"after!")
                    os.fsync(writer)
                finally:
                    os.close(writer)
                os.utime(
                    target,
                    ns=(original.st_atime_ns, original.st_mtime_ns + 1_000_000_000),
                )
            return metadata

        with mock.patch.object(core.os, "fstat", side_effect=mutate_after_first_file_fstat):
            with self.assertRaisesRegex(core.EvidenceError, "changed while it was read"):
                core.regular_bytes(target, "same-size race fixture")
        final = target.stat()
        self.assertTrue(triggered)
        self.assertEqual((final.st_dev, final.st_ino, final.st_size), (original.st_dev, original.st_ino, original.st_size))

    def test_directory_snapshot_rejects_late_insertion_and_final_inventory_change(self) -> None:
        directory = self.root / "late-insertion"
        directory.mkdir()
        (directory / "first").write_bytes(b"first")
        expected = directory.stat()
        real_listdir = os.listdir
        real_fstat = os.fstat
        triggered = False

        def insert_after_list(descriptor: int):
            nonlocal triggered
            names = real_listdir(descriptor)
            metadata = real_fstat(descriptor)
            if (
                not triggered
                and metadata.st_dev == expected.st_dev
                and metadata.st_ino == expected.st_ino
            ):
                triggered = True
                late = os.open(
                    "late", os.O_WRONLY | os.O_CREAT | os.O_EXCL,
                    0o644, dir_fd=descriptor,
                )
                try:
                    os.write(late, b"late")
                finally:
                    os.close(late)
            return names

        with mock.patch.object(core.os, "listdir", side_effect=insert_after_list):
            with self.assertRaisesRegex(core.EvidenceError, "directory changed during traversal"):
                core.tree_inventory(directory)
        self.assertTrue(triggered)

        original_check = core._verify_original_project_verification

        def mutate_before_return(**kwargs):
            result = original_check(**kwargs)
            (self.bundle / "late-before-return").write_bytes(b"late")
            return result

        with mock.patch.object(
            core, "_verify_original_project_verification", side_effect=mutate_before_return
        ):
            with self.assertRaisesRegex(core.EvidenceError, "changed between authenticated inventory"):
                core.verify_bundle(self.bundle)


if __name__ == "__main__":
    unittest.main()
