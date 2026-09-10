#!/usr/bin/env python3
"""Compose additive Stage 18 finalization records without changing the panel."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import secrets
import stat
import tempfile

import verify_stage18_finalization_correction as verify


class ComposerError(RuntimeError):
    pass


def atomic_bytes(path: Path, data: bytes, *, candidate_names: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    parent = path.parent
    parent_stat = os.lstat(parent)
    if not stat.S_ISDIR(parent_stat.st_mode) or parent.is_symlink():
        raise ComposerError(f"additive output parent is not a regular directory: {parent}")
    if not hasattr(os, "O_NOFOLLOW") or not hasattr(os, "O_DIRECTORY"):
        raise ComposerError("platform lacks required no-follow directory-safe creation flags")
    directory_flags = os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW
    directory_fd = os.open(parent, directory_flags)
    temporary: Path | None = None
    temporary_fd: int | None = None
    created = False
    supplied = iter(candidate_names) if candidate_names is not None else None
    try:
        for _ in range(128):
            if supplied is None:
                name = f".{path.name}.{secrets.token_hex(16)}.tmp"
            else:
                try:
                    name = next(supplied)
                except StopIteration as error:
                    raise ComposerError("exclusive temporary-name candidates were exhausted") from error
            if Path(name).name != name or name in {".", ".."}:
                raise ComposerError("unsafe exclusive temporary name")
            candidate = parent / name
            flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW
            try:
                temporary_fd = os.open(candidate, flags, 0o600)
            except FileExistsError:
                continue
            temporary = candidate
            created = True
            break
        if temporary_fd is None or temporary is None:
            raise ComposerError("could not create an exclusive fresh temporary file")
        before = os.lstat(temporary)
        opened = os.fstat(temporary_fd)
        if (
            not stat.S_ISREG(before.st_mode) or not stat.S_ISREG(opened.st_mode)
            or before.st_nlink != 1 or opened.st_nlink != 1
            or (before.st_dev, before.st_ino) != (opened.st_dev, opened.st_ino)
        ):
            raise ComposerError("exclusive temporary file identity is unsafe")
        view = memoryview(data)
        while view:
            written = os.write(temporary_fd, view)
            if written <= 0:
                raise ComposerError("short write to exclusive temporary file")
            view = view[written:]
        os.fsync(temporary_fd)
        os.close(temporary_fd)
        temporary_fd = None
        after = os.lstat(temporary)
        if (
            not stat.S_ISREG(after.st_mode) or after.st_nlink != 1
            or (after.st_dev, after.st_ino) != (before.st_dev, before.st_ino)
        ):
            raise ComposerError("exclusive temporary file changed before publication")
        if path.is_symlink():
            raise ComposerError(f"atomic destination is a symlink: {path}")
        if path.exists():
            destination = os.lstat(path)
            if not stat.S_ISREG(destination.st_mode) or destination.st_nlink != 1:
                raise ComposerError(f"atomic destination is not a single-link regular file: {path}")
        os.replace(temporary, path)
        created = False
        published = os.lstat(path)
        if not stat.S_ISREG(published.st_mode) or published.st_nlink != 1:
            raise ComposerError("published additive output is not a single-link regular file")
        os.fsync(directory_fd)
    finally:
        if temporary_fd is not None:
            os.close(temporary_fd)
        if created and temporary is not None:
            try:
                temporary.unlink()
            except FileNotFoundError:
                pass
        os.close(directory_fd)


def atomic_json(path: Path, value: dict, *, candidate_names: list[str] | None = None) -> None:
    atomic_bytes(path, verify.output_bytes(value), candidate_names=candidate_names)


def safe_additive_output(path: Path, panel: Path) -> Path:
    resolved = path.resolve()
    try:
        resolved.relative_to(panel.resolve(strict=True))
    except ValueError:
        pass
    else:
        raise ComposerError(f"additive output cannot be inside the immutable panel: {resolved}")
    if resolved in {Path("/").resolve(), verify.REPO.resolve(), panel.resolve()}:
        raise ComposerError(f"unsafe additive output: {resolved}")
    if path.is_symlink() or (path.exists() and not path.is_file()):
        raise ComposerError(f"additive output is not a regular file: {resolved}")
    if path.exists() and os.lstat(path).st_nlink != 1:
        raise ComposerError(f"additive output has multiple hard links: {resolved}")
    return resolved


def compose(panel: Path, amendment: Path, summary_path: Path,
            certificate_path: Path, results_path: Path,
            manifest_path: Path) -> dict:
    panel = panel.resolve(strict=True)
    summary_path = safe_additive_output(summary_path, panel)
    certificate_path = safe_additive_output(certificate_path, panel)
    results_path = safe_additive_output(results_path, panel)
    manifest_path = safe_additive_output(manifest_path, panel)
    outputs = (summary_path, certificate_path, results_path, manifest_path)
    if len(set(outputs)) != len(outputs):
        raise ComposerError("additive output paths collide")
    before_rows, before_bytes, before_hash = verify.panel_inventory(panel)
    summary, certificate, results, manifest = verify.compose_bundle_expected(
        panel, amendment.resolve(strict=True),
    )
    for path, data in (
        (summary_path, verify.output_bytes(summary)),
        (certificate_path, verify.output_bytes(certificate)),
        (results_path, results),
        (manifest_path, verify.output_bytes(manifest)),
    ):
        if not path.exists() or path.read_bytes() != data:
            atomic_bytes(path, data)
    after_rows, after_bytes, after_hash = verify.panel_inventory(panel)
    if (after_rows, after_bytes, after_hash) != (before_rows, before_bytes, before_hash):
        raise ComposerError("immutable panel changed during additive composition")
    verified = verify.verify_outputs(
        panel, amendment.resolve(), summary_path, certificate_path,
        results_path, manifest_path,
    )
    return {
        "composition": "pass",
        "status": summary["status"],
        "summary": str(summary_path),
        "summary_sha256": verify.sha256_file(summary_path),
        "certificate": str(certificate_path),
        "certificate_sha256": verify.sha256_file(certificate_path),
        "results": str(results_path),
        "results_sha256": verify.sha256_file(results_path),
        "manifest": str(manifest_path),
        "manifest_sha256": verify.sha256_file(manifest_path),
        "panel_tree_sha256": after_hash,
        "panel_files_modified": verified["existing_panel_files_modified_by_correction"],
        "scientific_tasks_added": verified["scientific_tasks_added_by_correction"],
    }


def self_test() -> dict:
    summary, certificate, results, manifest = verify.compose_bundle_expected()
    with tempfile.TemporaryDirectory(prefix="stage18-finalization-composer-") as temporary:
        root = Path(temporary)
        summary_path = root / "summary.json"
        certificate_path = root / "certificate.json"
        results_path = root / "STAGE18_RESULTS.md"
        manifest_path = root / "manifest.json"
        result = compose(
            verify.PANEL, verify.AMENDMENT, summary_path, certificate_path,
            results_path, manifest_path,
        )
        if (
            verify.read_json(summary_path) != summary
            or verify.read_json(certificate_path) != certificate
            or results_path.read_bytes() != results
            or verify.read_json(manifest_path) != manifest
        ):
            raise AssertionError("temporary additive composition changed its payload")
        symlink_victim = root / "symlink-victim"
        symlink_victim.write_bytes(b"symlink-victim")
        symlink_name = ".preexisting-symlink.tmp"
        (root / symlink_name).symlink_to(symlink_victim)
        atomic_bytes(
            root / "symlink-output", b"safe",
            candidate_names=[symlink_name, ".fresh-after-symlink.tmp"],
        )
        if symlink_victim.read_bytes() != b"symlink-victim":
            raise AssertionError("preexisting temporary symlink target was modified")
        if not (root / symlink_name).is_symlink():
            raise AssertionError("preexisting temporary symlink was replaced")
        hardlink_victim = root / "hardlink-victim"
        hardlink_victim.write_bytes(b"hardlink-victim")
        hardlink_name = ".preexisting-hardlink.tmp"
        os.link(hardlink_victim, root / hardlink_name)
        atomic_bytes(
            root / "hardlink-output", b"safe",
            candidate_names=[hardlink_name, ".fresh-after-hardlink.tmp"],
        )
        if hardlink_victim.read_bytes() != b"hardlink-victim":
            raise AssertionError("preexisting temporary hardlink target was modified")
        if not os.path.samefile(hardlink_victim, root / hardlink_name):
            raise AssertionError("preexisting temporary hardlink was replaced")
        if verify.panel_inventory(verify.PANEL)[2] != verify.PANEL_TREE_SHA256:
            raise AssertionError("atomic-writer controls changed the immutable panel")
        try:
            safe_additive_output(verify.PANEL / "forbidden.json", verify.PANEL)
        except ComposerError:
            pass
        else:
            raise AssertionError("composer accepted an output inside the immutable panel")
        return {
            "self_test": "pass", "checks": 9,
            "status": result["status"], "panel_tree_sha256": result["panel_tree_sha256"],
            "panel_files_modified": result["panel_files_modified"],
            "scientific_tasks_added": result["scientific_tasks_added"],
        }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel", type=Path, default=verify.PANEL)
    parser.add_argument("--amendment", type=Path, default=verify.AMENDMENT)
    parser.add_argument("--summary", type=Path, default=verify.CORRECTED_SUMMARY)
    parser.add_argument("--certificate", type=Path, default=verify.CERTIFICATE)
    parser.add_argument("--results", type=Path, default=verify.RESULTS)
    parser.add_argument("--manifest", type=Path, default=verify.BUNDLE_MANIFEST)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    try:
        result = self_test() if args.self_test else compose(
            args.panel, args.amendment, args.summary, args.certificate,
            args.results, args.manifest,
        )
        print(json.dumps(result, indent=2, sort_keys=True))
    except (ComposerError, verify.CorrectionError, OSError) as error:
        parser.exit(1, f"Stage 18 finalization composition failed: {error}\n")


if __name__ == "__main__":
    main()
