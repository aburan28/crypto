#!/usr/bin/env python3
"""Compact, package, and verify Stage-23 terminal evidence without path rewriting."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import math
import os
from pathlib import Path, PurePosixPath
import stat
import subprocess
from typing import Any, Iterable


MANIFEST_SCHEMA = "koblitz_stage23_terminal_evidence_manifest.v1"
SEAL_SCHEMA = "koblitz_stage23_terminal_evidence_seal.v1"
PATH_MAP_SCHEMA = "koblitz_stage23_archive_path_map.v1"
VERIFY_SCHEMA = "koblitz_stage23_terminal_evidence_verification.v1"
RUN_SEAL_SCHEMA = "koblitz_unknown_scalar_run_seal.v1"
RUN_SUMMARY_SCHEMA = "koblitz_unknown_scalar_run_summary.v1"
PROJECT_VERIFICATION_SCHEMA = "koblitz_unknown_scalar_verification.v1"
PROJECT_VERIFICATION_SEAL_SCHEMA = "koblitz_unknown_scalar_verification_seal.v1"
HEX40 = set("0123456789abcdef")
DIRECT_SOURCE_PATHS = (
    "Cargo.toml",
    "Cargo.lock",
    "research/sat_factor_base_review_20260908/continuation-05-sota-gates/stage-20-rust-build/Cargo.lock",
    "examples/koblitz_public_factor_base_discovery.rs",
    "examples/koblitz_unknown_scalar_panel.rs",
    "scripts/run_koblitz_unknown_scalar_panel.py",
    "scripts/process_meter.py",
    "scripts/run_koblitz_relation_yield_bridge.py",
    "src/cryptanalysis/koblitz_index_calculus.rs",
    "src/cryptanalysis/koblitz_pdp_phase_a.rs",
    "src/cryptanalysis/sat.rs",
    "src/cryptanalysis/semaev_sat.rs",
)
FROZEN_LOCK_RELATIVE = (
    "research/sat_factor_base_review_20260908/continuation-05-sota-gates/"
    "stage-20-rust-build/Cargo.lock"
)
BINARY_NAMES = (
    "koblitz_public_factor_base_discovery",
    "koblitz_unknown_scalar_panel",
)
ROOT_LAYOUT = {
    "run": ("sealed_run", "original/run"),
    "outer": ("outer_process_receipts", "original/outer"),
    "project_verification": ("project_verification", "original/project-verification"),
    "source_repository": ("committed_source_tree", "source-tree"),
}
EXPECTED_ENVIRONMENT = {
    "CARGO_BUILD_JOBS": "1",
    "LANG": "C",
    "LC_ALL": "C",
    "MKL_NUM_THREADS": "1",
    "OMP_NUM_THREADS": "1",
    "OPENBLAS_NUM_THREADS": "1",
    "PATH": "/usr/bin:/bin:/usr/sbin:/sbin:/opt/homebrew/bin",
    "RAYON_NUM_THREADS": "1",
    "TMPDIR": "/tmp",
    "TZ": "UTC",
    "VECLIB_MAXIMUM_THREADS": "1",
}
PENDING_MATH_REPLAY = [
    "domain-separated target generation and finite-field subgroup membership",
    "public factor-base discovery and materialization",
    "every relation-attempt target and decomposition witness over the curve",
    "curve-validity of every retained relation row",
    "IC recovered-scalar elliptic-curve point verification",
    "signed-Frobenius rho walk and recovered-scalar point verification",
]
ORIGINAL_PORTABLE_BLOCKER = (
    "run identities retain absolute checkout paths; no archive-relative source and executable replay exists"
)


class EvidenceError(RuntimeError):
    """Fail-closed terminal-evidence error."""


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode()


def canonical_sha256(value: Any) -> str:
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def git_object_oid(kind: str, data: bytes) -> str:
    framed = f"{kind} {len(data)}\0".encode() + data
    return hashlib.sha1(framed).hexdigest()


def require_exact(value: Any, keys: set[str], context: str) -> dict[str, Any]:
    if not isinstance(value, dict) or set(value) != keys:
        raise EvidenceError(f"{context} uses an unexpected schema")
    return value


def require_hex(value: Any, length: int, context: str) -> str:
    if (
        not isinstance(value, str)
        or len(value) != length
        or any(character not in HEX40 for character in value)
    ):
        raise EvidenceError(f"{context} must be {length}-character lowercase hexadecimal")
    return value


def safe_relative(value: Any, context: str) -> str:
    if not isinstance(value, str) or not value or "\\" in value or "\0" in value:
        raise EvidenceError(f"{context} is not a portable archive-relative path")
    path = PurePosixPath(value)
    if path.is_absolute() or value != path.as_posix() or any(part in {"", ".", ".."} for part in path.parts):
        raise EvidenceError(f"{context} is not a normalized archive-relative path")
    if any("\n" in part or "\r" in part or "\t" in part for part in path.parts):
        raise EvidenceError(f"{context} contains a control separator")
    return value


def is_build_target_path(relative: str) -> bool:
    """Classify the reserved build target itself and every descendant."""
    return relative == "build-target" or relative.startswith("build-target/")


def normalized_absolute(value: Any, context: str) -> str:
    if not isinstance(value, str) or not value.startswith("/") or "\0" in value:
        raise EvidenceError(f"{context} must be an absolute path")
    if os.path.normpath(value) != value or any(part == ".." for part in Path(value).parts):
        raise EvidenceError(f"{context} must be normalized")
    return value


def safe_join(root: Path, relative: str, context: str) -> Path:
    relative = safe_relative(relative, context)
    candidate = root.joinpath(*PurePosixPath(relative).parts)
    try:
        candidate.relative_to(root)
    except ValueError as error:
        raise EvidenceError(f"{context} escapes its archive root") from error
    return candidate


def _directory_flags() -> int:
    flags = os.O_RDONLY
    flags |= getattr(os, "O_DIRECTORY", 0)
    flags |= getattr(os, "O_CLOEXEC", 0)
    flags |= getattr(os, "O_NOFOLLOW", 0)
    return flags


def _stable_stat_signature(metadata: os.stat_result) -> tuple[int, int, int, int, int, int, int]:
    return (
        metadata.st_mode,
        metadata.st_nlink,
        metadata.st_dev,
        metadata.st_ino,
        metadata.st_size,
        metadata.st_mtime_ns,
        metadata.st_ctime_ns,
    )


def _open_absolute_directory(path: Path, context: str, *, create: bool = False) -> int:
    value = normalized_absolute(str(path), context)
    parts = Path(value).parts[1:]
    descriptor = os.open("/", _directory_flags())
    try:
        for part in parts:
            try:
                child = os.open(part, _directory_flags(), dir_fd=descriptor)
            except FileNotFoundError:
                if not create:
                    raise
                os.mkdir(part, 0o755, dir_fd=descriptor)
                child = os.open(part, _directory_flags(), dir_fd=descriptor)
            metadata = os.fstat(child)
            if not stat.S_ISDIR(metadata.st_mode):
                os.close(child)
                raise EvidenceError(f"{context} component is not a directory: {part}")
            os.close(descriptor)
            descriptor = child
        return descriptor
    except BaseException:
        os.close(descriptor)
        raise


def _read_regular_at(parent_fd: int, name: str, context: str) -> bytes:
    flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0)
    descriptor = os.open(name, flags, dir_fd=parent_fd)
    try:
        metadata = os.fstat(descriptor)
        if not stat.S_ISREG(metadata.st_mode):
            raise EvidenceError(f"{context} must be a regular non-symlink file")
        if metadata.st_nlink != 1:
            raise EvidenceError(f"{context} must not be hard-linked")
        chunks: list[bytes] = []
        while True:
            chunk = os.read(descriptor, 1024 * 1024)
            if not chunk:
                break
            chunks.append(chunk)
        data = b"".join(chunks)
        final = os.fstat(descriptor)
        if (
            _stable_stat_signature(final) != _stable_stat_signature(metadata)
            or len(data) != metadata.st_size
        ):
            raise EvidenceError(f"{context} changed while it was read")
        return data
    finally:
        os.close(descriptor)


def regular_bytes(path: Path, context: str) -> bytes:
    path = Path(path)
    if not path.is_absolute():
        path = Path.cwd() / path
    normalized_absolute(str(path), context)
    parent_fd = _open_absolute_directory(path.parent, f"{context} parent")
    try:
        return _read_regular_at(parent_fd, path.name, f"{context} {path}")
    except OSError as error:
        raise EvidenceError(f"cannot read {context} {path}: {error}") from error
    finally:
        os.close(parent_fd)


def real_directory(path: Path, context: str) -> Path:
    path = Path(path)
    if not path.is_absolute():
        path = Path.cwd() / path
    normalized_absolute(str(path), context)
    try:
        descriptor = _open_absolute_directory(path, context)
    except OSError as error:
        raise EvidenceError(f"cannot open {context} {path}: {error}") from error
    try:
        if not stat.S_ISDIR(os.fstat(descriptor).st_mode):
            raise EvidenceError(f"{context} must be a real directory")
    finally:
        os.close(descriptor)
    return path


def read_json(path: Path, context: str) -> tuple[dict[str, Any], bytes]:
    data = regular_bytes(path, context)

    def reject_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        value: dict[str, Any] = {}
        for key, child in pairs:
            if key in value:
                raise EvidenceError(f"duplicate JSON key {key!r} in {context}")
            value[key] = child
        return value

    try:
        value = json.loads(data, object_pairs_hook=reject_duplicates)
    except EvidenceError:
        raise
    except json.JSONDecodeError as error:
        raise EvidenceError(f"invalid JSON in {context}: {error}") from error
    if not isinstance(value, dict):
        raise EvidenceError(f"{context} must contain one JSON object")
    return value, data


def identity(path: Path, *, recorded_path: str | None = None, context: str = "file") -> dict[str, Any]:
    data = regular_bytes(path, context)
    return {
        "path": recorded_path if recorded_path is not None else str(path),
        "bytes": len(data),
        "sha256": sha256(data),
    }


def validate_identity_shape(value: Any, context: str) -> dict[str, Any]:
    value = require_exact(value, {"path", "bytes", "sha256"}, context)
    if not isinstance(value["path"], str):
        raise EvidenceError(f"{context}.path must be a string")
    if isinstance(value["bytes"], bool) or not isinstance(value["bytes"], int) or value["bytes"] < 0:
        raise EvidenceError(f"{context}.bytes must be a nonnegative integer")
    require_hex(value["sha256"], 64, f"{context}.sha256")
    return value


def validate_inventory_records(value: Any, context: str) -> list[dict[str, Any]]:
    if not isinstance(value, list):
        raise EvidenceError(f"{context} must be a list")
    records: list[dict[str, Any]] = []
    paths: set[str] = set()
    for index, raw in enumerate(value):
        record = validate_identity_shape(raw, f"{context}[{index}]")
        path = safe_relative(record["path"], f"{context}[{index}].path")
        if path in paths:
            raise EvidenceError(f"{context} contains duplicate path {path}")
        paths.add(path)
        records.append(record)
    if records != sorted(records, key=lambda item: item["path"]):
        raise EvidenceError(f"{context} is not sorted by path")
    return records


def tree_inventory(root: Path, *, excluded: Iterable[str] = ()) -> list[dict[str, Any]]:
    excluded_set = set(excluded)
    records: list[dict[str, Any]] = []
    root_fd = _open_absolute_directory(root, "artifact-tree root")

    def walk(directory_fd: int, prefix: PurePosixPath | None = None) -> None:
        directory_before = os.fstat(directory_fd)
        names_before = sorted(os.listdir(directory_fd))
        for name in names_before:
            if name in {".", ".."} or "/" in name or "\0" in name:
                raise EvidenceError("artifact tree contains an invalid directory entry")
            relative = name if prefix is None else (prefix / name).as_posix()
            metadata = os.stat(name, dir_fd=directory_fd, follow_symlinks=False)
            if stat.S_ISLNK(metadata.st_mode):
                raise EvidenceError(f"artifact tree contains symlink {relative}")
            if stat.S_ISDIR(metadata.st_mode):
                child = os.open(name, _directory_flags(), dir_fd=directory_fd)
                try:
                    opened = os.fstat(child)
                    if opened.st_dev != metadata.st_dev or opened.st_ino != metadata.st_ino:
                        raise EvidenceError(f"artifact directory changed during traversal: {relative}")
                    walk(child, PurePosixPath(relative))
                finally:
                    os.close(child)
                continue
            if not stat.S_ISREG(metadata.st_mode):
                raise EvidenceError(f"artifact tree contains inadmissible file {relative}")
            data = _read_regular_at(directory_fd, name, f"artifact file {relative}")
            if relative not in excluded_set:
                records.append({"path": relative, "bytes": len(data), "sha256": sha256(data)})
        names_after = sorted(os.listdir(directory_fd))
        directory_after = os.fstat(directory_fd)
        if (
            names_after != names_before
            or _stable_stat_signature(directory_after) != _stable_stat_signature(directory_before)
        ):
            label = prefix.as_posix() if prefix is not None else "."
            raise EvidenceError(f"artifact directory changed during traversal: {label}")

    try:
        walk(root_fd)
    finally:
        os.close(root_fd)
    return sorted(records, key=lambda record: record["path"])


def write_new(path: Path, data: bytes, mode: int = 0o644) -> None:
    path = Path(path)
    if not path.is_absolute():
        path = Path.cwd() / path
    normalized_absolute(str(path), "new output file")
    parent_fd = _open_absolute_directory(path.parent, "new output parent", create=True)
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    descriptor = -1
    try:
        try:
            descriptor = os.open(path.name, flags, mode, dir_fd=parent_fd)
        except FileExistsError as error:
            raise EvidenceError(f"refusing to overwrite {path}") from error
        metadata = os.fstat(descriptor)
        if not stat.S_ISREG(metadata.st_mode) or metadata.st_nlink != 1:
            raise EvidenceError(f"new output is not a private regular file: {path}")
        os.fchmod(descriptor, mode)
        stream = os.fdopen(descriptor, "wb")
        descriptor = -1
        with stream:
            stream.write(data)
            stream.flush()
            os.fsync(stream.fileno())
    finally:
        if descriptor >= 0:
            os.close(descriptor)
        os.close(parent_fd)


def write_json_new(path: Path, value: Any) -> None:
    write_new(path, json.dumps(value, indent=2, sort_keys=True).encode() + b"\n")


def safe_new_output(path: Path) -> Path:
    path = Path(path)
    if not path.is_absolute():
        path = Path.cwd() / path
    normalized_absolute(str(path), "output")
    if path == Path("/"):
        raise EvidenceError(f"output must be a new, non-root path: {path}")
    parent_fd = _open_absolute_directory(path.parent, "output parent")
    try:
        try:
            os.stat(path.name, dir_fd=parent_fd, follow_symlinks=False)
        except FileNotFoundError:
            pass
        else:
            raise EvidenceError(f"output must not already exist: {path}")
    finally:
        os.close(parent_fd)
    return path


def create_new_directory(path: Path, context: str) -> None:
    parent_fd = _open_absolute_directory(path.parent, f"{context} parent")
    try:
        os.mkdir(path.name, 0o755, dir_fd=parent_fd)
        descriptor = os.open(path.name, _directory_flags(), dir_fd=parent_fd)
        try:
            if not stat.S_ISDIR(os.fstat(descriptor).st_mode):
                raise EvidenceError(f"{context} is not a directory")
        finally:
            os.close(descriptor)
    except FileExistsError as error:
        raise EvidenceError(f"{context} already exists: {path}") from error
    finally:
        os.close(parent_fd)


def copy_tree_bytes(source: Path, destination: Path, records: list[dict[str, Any]]) -> None:
    for record in records:
        source_path = safe_join(source, record["path"], "source inventory path")
        data = regular_bytes(source_path, f"source file {record['path']}")
        if len(data) != record["bytes"] or sha256(data) != record["sha256"]:
            raise EvidenceError(f"source inventory digest changed for {record['path']}")
        mode = 0o755 if record["path"].startswith("binaries/") else 0o644
        write_new(safe_join(destination, record["path"], "destination inventory path"), data, mode)


def validate_self_hash(value: dict[str, Any], field: str, context: str) -> None:
    payload = dict(value)
    claimed = payload.pop(field, None)
    if claimed != canonical_sha256(payload):
        raise EvidenceError(f"{context} self-hash is invalid")


def run_command(command: list[str], context: str, *, cwd: Path | None = None) -> bytes:
    completed = subprocess.run(command, cwd=cwd, capture_output=True, check=False)
    if completed.returncode != 0:
        stderr = completed.stderr.decode(errors="replace").strip()
        raise EvidenceError(f"{context} failed: {stderr}")
    return completed.stdout


def parse_git_tree(repo: Path, commit: str) -> tuple[list[dict[str, Any]], bytes]:
    raw = run_command(
        ["git", "-C", str(repo), "ls-tree", "-r", "-z", "--full-tree", "-l", commit],
        "Git tree inventory",
    )
    records: list[dict[str, Any]] = []
    for item in raw.split(b"\0"):
        if not item:
            continue
        try:
            header, raw_path = item.split(b"\t", 1)
            mode, kind, oid, size = header.decode("ascii").split(" ", 3)
            size = size.strip()
            path = raw_path.decode("utf-8")
        except (ValueError, UnicodeDecodeError) as error:
            raise EvidenceError("Git tree contains an unsupported entry") from error
        safe_relative(path, "Git tree path")
        if kind != "blob" or mode not in {"100644", "100755"} or not size.isdigit():
            raise EvidenceError(f"Git tree entry is not a regular file: {path}")
        require_hex(oid, 40, f"Git blob {path}")
        records.append({"path": path, "mode": mode, "bytes": int(size), "blob_oid": oid})
    if records != sorted(records, key=lambda row: row["path"]):
        raise EvidenceError("Git tree inventory is not sorted")
    if len({record["path"] for record in records}) != len(records):
        raise EvidenceError("Git tree inventory contains duplicate paths")
    listing = run_command(
        ["git", "-C", str(repo), "ls-tree", "-r", "--full-tree", commit],
        "Git tree listing",
    )
    return records, listing


def git_blob_bytes(repo: Path, oid: str, expected_size: int, context: str) -> bytes:
    data = run_command(["git", "-C", str(repo), "cat-file", "blob", oid], context)
    if len(data) != expected_size or git_object_oid("blob", data) != oid:
        raise EvidenceError(f"{context} does not match its Git object identity")
    return data


def compute_tree_oid(records: list[dict[str, Any]]) -> str:
    root: dict[str, Any] = {}
    for record in records:
        node = root
        parts = PurePosixPath(record["path"]).parts
        for part in parts[:-1]:
            child = node.setdefault(part, {})
            if not isinstance(child, dict):
                raise EvidenceError("source tree contains a file/directory collision")
            node = child
        if parts[-1] in node:
            raise EvidenceError("source tree contains a duplicate entry")
        node[parts[-1]] = (record["mode"], record["blob_oid"])

    def digest(node: dict[str, Any]) -> str:
        rendered: list[tuple[bytes, bytes]] = []
        for name, child in node.items():
            encoded = name.encode("utf-8")
            if isinstance(child, dict):
                oid = digest(child)
                rendered.append((encoded + b"/", b"40000 " + encoded + b"\0" + bytes.fromhex(oid)))
            else:
                mode, oid = child
                rendered.append((encoded, mode.encode() + b" " + encoded + b"\0" + bytes.fromhex(oid)))
        body = b"".join(entry for _, entry in sorted(rendered, key=lambda pair: pair[0]))
        return git_object_oid("tree", body)

    return digest(root)


def materialize_git_tree(repo: Path, records: list[dict[str, Any]], destination: Path) -> list[dict[str, Any]]:
    completed: list[dict[str, Any]] = []
    for record in records:
        data = git_blob_bytes(repo, record["blob_oid"], record["bytes"], f"Git blob {record['path']}")
        target = safe_join(destination, record["path"], "source-tree destination")
        write_new(target, data, 0o755 if record["mode"] == "100755" else 0o644)
        completed.append({**record, "sha256": sha256(data)})
    return completed


@dataclass(frozen=True)
class ArchiveMap:
    bundle: Path
    roots: dict[str, dict[str, Any]]
    files: dict[str, dict[str, Any]]
    external_tools: dict[str, dict[str, Any]]

    @classmethod
    def from_manifest(cls, bundle: Path, value: Any) -> "ArchiveMap":
        value = require_exact(value, {"schema", "roots", "files", "external_tools"}, "path map")
        if value["schema"] != PATH_MAP_SCHEMA:
            raise EvidenceError("path map schema changed")
        roots_value = value["roots"]
        if not isinstance(roots_value, list):
            raise EvidenceError("path-map roots must be a list")
        roots: dict[str, dict[str, Any]] = {}
        for entry in roots_value:
            entry = require_exact(
                entry,
                {"role", "kind", "original_absolute", "archive_relative"},
                "path-map root",
            )
            role = entry["role"]
            if role not in ROOT_LAYOUT or role in roots:
                raise EvidenceError(f"unexpected or duplicate path-map root {role!r}")
            expected_kind, expected_archive = ROOT_LAYOUT[role]
            if entry["kind"] != expected_kind or entry["archive_relative"] != expected_archive:
                raise EvidenceError(f"path-map root {role} changed its type or archive location")
            normalized_absolute(entry["original_absolute"], f"path-map root {role}")
            safe_relative(entry["archive_relative"], f"path-map archive root {role}")
            roots[role] = entry
        if set(roots) != set(ROOT_LAYOUT):
            raise EvidenceError("path map omits a required typed root")

        files_value = value["files"]
        if not isinstance(files_value, list):
            raise EvidenceError("path-map files must be a list")
        files: dict[str, dict[str, Any]] = {}
        expected_file_roles = {
            "ignored_workspace_lock": ("ignored_workspace_dependency_lock", "source-extra/Cargo.lock"),
            "source_commit_object": ("git_commit_object", "source-metadata/commit.object"),
            "packager_source": ("portable_tool_source", "verifier/package_koblitz_stage23_terminal_evidence.py"),
            "core_source": ("portable_tool_source", "verifier/koblitz_stage23_terminal_evidence.py"),
            "verifier_source": ("portable_tool_source", "verifier/verify_koblitz_stage23_terminal_evidence.py"),
        }
        for entry in files_value:
            entry = require_exact(
                entry,
                {"role", "kind", "original_absolute", "archive_relative"},
                "path-map file",
            )
            role = entry["role"]
            if role not in expected_file_roles or role in files:
                raise EvidenceError(f"unexpected or duplicate path-map file {role!r}")
            expected_kind, expected_archive = expected_file_roles[role]
            if entry["kind"] != expected_kind or entry["archive_relative"] != expected_archive:
                raise EvidenceError(f"path-map file {role} changed its type or archive location")
            normalized_absolute(entry["original_absolute"], f"path-map file {role}")
            safe_relative(entry["archive_relative"], f"path-map archive file {role}")
            files[role] = entry
        if set(files) != set(expected_file_roles):
            raise EvidenceError("path map omits a required typed file")

        external = value["external_tools"]
        if not isinstance(external, list) or len(external) != 3:
            raise EvidenceError("path map must retain three external tool identities")
        roles: set[str] = set()
        external_tools: dict[str, dict[str, Any]] = {}
        for entry in external:
            entry = require_exact(
                entry,
                {"role", "kind", "original_absolute", "archived", "identity"},
                "external tool path",
            )
            role = entry["role"]
            if role not in {"python", "cargo", "rustc"} or role in roles:
                raise EvidenceError("external tool path roles changed")
            roles.add(role)
            external_tools[role] = entry
            if entry["kind"] != "external_tool_identity" or entry["archived"] is not False:
                raise EvidenceError("external tools must remain identity-only and unarchived")
            normalized_absolute(entry["original_absolute"], f"external tool {role}")
            identity_value = validate_identity_shape(entry["identity"], f"external tool {role} identity")
            if identity_value["path"] != entry["original_absolute"]:
                raise EvidenceError(f"external tool {role} path and identity differ")
        if roles != {"python", "cargo", "rustc"}:
            raise EvidenceError("path map omits an external tool role")
        return cls(bundle=bundle, roots=roots, files=files, external_tools=external_tools)

    def original(self, role: str, relative: str | None = None) -> str:
        root = self.roots[role]["original_absolute"]
        if relative is None:
            return root
        safe_relative(relative, f"{role} original child")
        return str(Path(root).joinpath(*PurePosixPath(relative).parts))

    def archived(self, role: str, relative: str | None = None) -> Path:
        root = safe_join(self.bundle, self.roots[role]["archive_relative"], f"{role} archive root")
        if relative is None:
            return root
        return safe_join(root, relative, f"{role} archive child")

    def archived_file(self, role: str) -> Path:
        return safe_join(self.bundle, self.files[role]["archive_relative"], f"{role} archive file")

    def validate_bound_file(
        self, value: Any, role: str, relative: str, context: str
    ) -> dict[str, Any]:
        value = validate_identity_shape(value, context)
        if value["path"] != self.original(role, relative):
            raise EvidenceError(f"{context} points to the wrong typed original path")
        actual = identity(
            self.archived(role, relative), recorded_path=value["path"], context=context
        )
        if actual != value:
            raise EvidenceError(f"{context} changed after its identity was recorded")
        return value


def validate_metrics(value: Any, context: str) -> dict[str, Any]:
    value = require_exact(
        value,
        {"command", "returncode", "watchdog_seconds", "timed_out", "orphan_group_terminated", "metrics"},
        context,
    )
    if not isinstance(value["command"], list) or not all(isinstance(item, str) for item in value["command"]):
        raise EvidenceError(f"{context}.command must be a string list")
    watchdog = value["watchdog_seconds"]
    if isinstance(watchdog, bool) or not isinstance(watchdog, (int, float)) or not math.isfinite(float(watchdog)) or watchdog <= 0:
        raise EvidenceError(f"{context}.watchdog_seconds is invalid")
    if isinstance(value["returncode"], bool) or not isinstance(value["returncode"], int):
        raise EvidenceError(f"{context}.returncode must be an integer")
    if not isinstance(value["timed_out"], bool) or not isinstance(value["orphan_group_terminated"], bool):
        raise EvidenceError(f"{context} termination flags must be Boolean")
    metrics = require_exact(
        value["metrics"],
        {"wall_seconds", "user_seconds", "system_seconds", "total_core_seconds", "single_core_seconds", "peak_rss_bytes", "meter"},
        f"{context}.metrics",
    )
    numbers = [metrics[field] for field in ("wall_seconds", "user_seconds", "system_seconds", "total_core_seconds", "single_core_seconds")]
    if any(isinstance(number, bool) or not isinstance(number, (int, float)) or not math.isfinite(float(number)) for number in numbers):
        raise EvidenceError(f"{context} metric is not finite numeric data")
    if metrics["wall_seconds"] <= 0 or metrics["user_seconds"] < 0 or metrics["system_seconds"] < 0:
        raise EvidenceError(f"{context} time metric is outside its range")
    if not math.isclose(metrics["total_core_seconds"], metrics["user_seconds"] + metrics["system_seconds"], rel_tol=0, abs_tol=1e-9):
        raise EvidenceError(f"{context} CPU ledger does not balance")
    if metrics["single_core_seconds"] != metrics["total_core_seconds"]:
        raise EvidenceError(f"{context} legacy CPU alias changed")
    if isinstance(metrics["peak_rss_bytes"], bool) or not isinstance(metrics["peak_rss_bytes"], int) or metrics["peak_rss_bytes"] < 0:
        raise EvidenceError(f"{context} RSS is invalid")
    if metrics["meter"] != "fresh-process getrusage(RUSAGE_CHILDREN)":
        raise EvidenceError(f"{context} meter identity changed")
    return value


def process_clean(value: dict[str, Any]) -> bool:
    return value["returncode"] == 0 and not value["timed_out"] and not value["orphan_group_terminated"]


# Minimal, project-independent unkeyed BLAKE3 replay.  The verifier must treat
# every archived source file as data; it never imports or executes bundle code.
_B3_IV = (
    0x6A09E667, 0xBB67AE85, 0x3C6EF372, 0xA54FF53A,
    0x510E527F, 0x9B05688C, 0x1F83D9AB, 0x5BE0CD19,
)
_B3_PERMUTATION = (2, 6, 3, 10, 7, 0, 4, 13, 1, 11, 12, 5, 9, 14, 15, 8)
_B3_CHUNK_START = 1
_B3_CHUNK_END = 2
_B3_PARENT = 4
_B3_ROOT = 8


def _b3_rotr(value: int, shift: int) -> int:
    return ((value >> shift) | (value << (32 - shift))) & 0xFFFFFFFF


def _b3_g(state: list[int], a: int, b: int, c: int, d: int, x: int, y: int) -> None:
    state[a] = (state[a] + state[b] + x) & 0xFFFFFFFF
    state[d] = _b3_rotr(state[d] ^ state[a], 16)
    state[c] = (state[c] + state[d]) & 0xFFFFFFFF
    state[b] = _b3_rotr(state[b] ^ state[c], 12)
    state[a] = (state[a] + state[b] + y) & 0xFFFFFFFF
    state[d] = _b3_rotr(state[d] ^ state[a], 8)
    state[c] = (state[c] + state[d]) & 0xFFFFFFFF
    state[b] = _b3_rotr(state[b] ^ state[c], 7)


def _b3_round(state: list[int], message: list[int]) -> None:
    _b3_g(state, 0, 4, 8, 12, message[0], message[1])
    _b3_g(state, 1, 5, 9, 13, message[2], message[3])
    _b3_g(state, 2, 6, 10, 14, message[4], message[5])
    _b3_g(state, 3, 7, 11, 15, message[6], message[7])
    _b3_g(state, 0, 5, 10, 15, message[8], message[9])
    _b3_g(state, 1, 6, 11, 12, message[10], message[11])
    _b3_g(state, 2, 7, 8, 13, message[12], message[13])
    _b3_g(state, 3, 4, 9, 14, message[14], message[15])


def _b3_words(block: bytes) -> list[int]:
    padded = block + bytes(64 - len(block))
    return [int.from_bytes(padded[index:index + 4], "little") for index in range(0, 64, 4)]


def _b3_compress(
    chaining_value: tuple[int, ...] | list[int],
    block_words: list[int],
    counter: int,
    block_len: int,
    flags: int,
) -> list[int]:
    state = list(chaining_value) + list(_B3_IV[:4]) + [
        counter & 0xFFFFFFFF, (counter >> 32) & 0xFFFFFFFF, block_len, flags,
    ]
    message = list(block_words)
    for _ in range(7):
        _b3_round(state, message)
        message = [message[index] for index in _B3_PERMUTATION]
    return [state[index] ^ state[index + 8] for index in range(8)] + [
        state[index + 8] ^ chaining_value[index] for index in range(8)
    ]


def _b3_chunk_output(chunk: bytes, counter: int) -> tuple[list[int], list[int], int, int, int]:
    blocks = [chunk[index:index + 64] for index in range(0, len(chunk), 64)] or [b""]
    chaining_value = list(_B3_IV)
    for index, block in enumerate(blocks[:-1]):
        flags = _B3_CHUNK_START if index == 0 else 0
        chaining_value = _b3_compress(
            chaining_value, _b3_words(block), counter, len(block), flags
        )[:8]
    final_index = len(blocks) - 1
    final_flags = _B3_CHUNK_END | (_B3_CHUNK_START if final_index == 0 else 0)
    final = blocks[-1]
    return chaining_value, _b3_words(final), counter, len(final), final_flags


def _b3_output_cv(output: tuple[list[int], list[int], int, int, int]) -> list[int]:
    return _b3_compress(*output)[:8]


def _b3_parent_output(left: list[int], right: list[int]) -> tuple[list[int], list[int], int, int, int]:
    return list(_B3_IV), left + right, 0, 64, _B3_PARENT


def blake3_bytes(data: bytes) -> bytes:
    chunks = [data[index:index + 1024] for index in range(0, len(data), 1024)] or [b""]
    stack: list[list[int]] = []
    for chunk_index, chunk in enumerate(chunks[:-1]):
        chaining_value = _b3_output_cv(_b3_chunk_output(chunk, chunk_index))
        total_chunks = chunk_index + 1
        while total_chunks & 1 == 0:
            chaining_value = _b3_output_cv(_b3_parent_output(stack.pop(), chaining_value))
            total_chunks >>= 1
        stack.append(chaining_value)
    output = _b3_chunk_output(chunks[-1], len(chunks) - 1)
    while stack:
        output = _b3_parent_output(stack.pop(), _b3_output_cv(output))
    words = _b3_compress(output[0], output[1], 0, output[3], output[4] | _B3_ROOT)
    return b"".join(word.to_bytes(4, "little") for word in words)[:32]


def json_blake3(value: Any) -> str:
    return blake3_bytes(canonical_bytes(value)).hex()


def stage23_protocol(value: dict[str, Any], profile: str) -> dict[str, Any]:
    if value.get("schema") != "koblitz_unknown_scalar_panel_protocol.v1":
        raise EvidenceError("wrong Stage-23 protocol schema")
    if value.get("status") != "frozen_before_production_execution":
        raise EvidenceError("Stage-23 protocol is not frozen")
    if profile not in {"production", "smoke"}:
        raise EvidenceError("Stage-23 profile is invalid")
    if value.get("targets", {}).get("count") != 5 or value.get("factor_base", {}).get("divisor_indices") != [0, 2]:
        raise EvidenceError("Stage-23 target panel or factor base changed")
    if value.get("index_calculus", {}).get("parallel_workers") != 1 or value.get("execution", {}).get("parallel_workers") != 1:
        raise EvidenceError("Stage-23 execution is not single-worker")
    return value


def stage23_profile_values(profile: str, frozen: dict[str, Any]) -> dict[str, Any]:
    if profile == "production":
        return {
            "targets": frozen["targets"]["count"],
            "n": frozen["curve"]["n"],
            "a": frozen["curve"]["a"],
            "dimension": frozen["factor_base"]["dimension"],
            "evidence_class": "public_synthetic_candidate_pending_independent_replay",
        }
    smoke = frozen["smoke"]
    return {
        "targets": smoke["targets"], "n": smoke["curve_n"], "a": smoke["curve_a"],
        "dimension": 4, "evidence_class": "operational_smoke",
    }


def validate_discovery_result(
    result: dict[str, Any], curve_a: int, values: dict[str, Any], frozen: dict[str, Any]
) -> None:
    if (
        result.get("schema") != "koblitz_public_factor_base_discovery.v1"
        or result.get("n") != values["n"] or result.get("a") != curve_a
        or result.get("m") != 2 or result.get("requested_dimension") != values["dimension"]
    ):
        raise EvidenceError(f"curve-a={curve_a} discovery parameters changed")
    forbidden = result.get("forbidden_inputs")
    expected_forbidden = {
        "target_constructed", "target_subgroup_enumerated", "discrete_log_labels_constructed",
        "relation_yield_used", "solver_timing_used",
    }
    if not isinstance(forbidden, dict) or set(forbidden) != expected_forbidden or any(child is not False for child in forbidden.values()):
        raise EvidenceError(f"curve-a={curve_a} discovery used a forbidden input")
    candidates = result.get("candidates")
    if not isinstance(candidates, list) or not candidates:
        raise EvidenceError(f"curve-a={curve_a} discovery emitted no candidates")
    for candidate in candidates:
        if not isinstance(candidate, dict):
            raise EvidenceError(f"curve-a={curve_a} discovery candidate is malformed")
        if (
            not isinstance(candidate.get("divisor_indices"), list)
            or not candidate["divisor_indices"]
            or not all(type(index) is int and index >= 0 for index in candidate["divisor_indices"])
            or candidate.get("dimension") != values["dimension"]
        ):
            raise EvidenceError(f"curve-a={curve_a} discovery candidate basis is malformed")
        for field in (
            "divisor_polynomial", "abscissae", "rational_points",
            "signed_frobenius_orbits_before_projection", "projected_signed_frobenius_orbits",
        ):
            if type(candidate.get(field)) is not int or candidate[field] < 0:
                raise EvidenceError(f"curve-a={curve_a} discovery field {field} is invalid")
    if result.get("selected") not in candidates:
        raise EvidenceError(f"curve-a={curve_a} selected candidate is absent")
    if curve_a == frozen["curve"]["a"]:
        matches = [candidate for candidate in candidates if candidate.get("divisor_indices") == frozen["factor_base"]["divisor_indices"]]
        if not matches:
            raise EvidenceError("public discovery omitted the frozen K0 divisor base")
        if values["n"] == frozen["curve"]["n"]:
            candidate = matches[0]
            expected = frozen["factor_base"]
            if (
                candidate.get("divisor_polynomial") != expected["divisor_polynomial"]
                or candidate.get("rational_points") != expected["rational_points"]
                or candidate.get("signed_frobenius_orbits_before_projection") != expected["signed_frobenius_orbits_before_projection"]
                or candidate.get("projected_signed_frobenius_orbits") != expected["projected_signed_frobenius_columns"]
            ):
                raise EvidenceError("public discovery changed the production factor base")


def recursively_forbid_scalar_labels(value: Any) -> None:
    forbidden = {"secret", "target_scalar", "known_scalar", "planted_scalar"}
    if isinstance(value, dict):
        overlap = forbidden.intersection(value)
        if overlap:
            raise EvidenceError(f"forbidden scalar-label keys: {sorted(overlap)}")
        for child in value.values():
            recursively_forbid_scalar_labels(child)
    elif isinstance(value, list):
        for child in value:
            recursively_forbid_scalar_labels(child)


def validate_targets(result: dict[str, Any], profile: str, expected: int) -> list[dict[str, Any]]:
    if result.get("schema") != "koblitz_unknown_scalar_target_panel.v1" or result.get("status") != "complete":
        raise EvidenceError("target panel is not complete")
    if result.get("target_scalar_constructed_or_recorded") is not False or result.get("factor_base_log_labels_constructed_or_recorded") is not False:
        raise EvidenceError("target panel constructed forbidden scalar labels")
    recursively_forbid_scalar_labels(result)
    identity_value = result.get("identity", {})
    if json_blake3(identity_value) != result.get("identity_blake3"):
        raise EvidenceError("target-panel identity hash is invalid")
    if identity_value.get("profile") != profile or len(identity_value.get("targets", [])) != expected:
        raise EvidenceError("target-panel profile or count changed")
    targets = identity_value["targets"]
    ids = [target.get("target_id") for target in targets]
    points = [(target.get("point", {}).get("x"), target.get("point", {}).get("y")) for target in targets]
    if len(set(ids)) != expected or len(set(points)) != expected:
        raise EvidenceError("target panel contains duplicate identities or points")
    return targets


def canonical_residue(value: Any, modulus: int, context: str) -> int:
    if not isinstance(value, str):
        raise EvidenceError(f"{context} must be a canonical decimal string")
    try:
        parsed = int(value)
    except ValueError as error:
        raise EvidenceError(f"{context} is not a decimal integer") from error
    if value != str(parsed) or not 0 <= parsed < modulus:
        raise EvidenceError(f"{context} is not a canonical residue")
    return parsed


def relation_arithmetic(profile: str, frozen: dict[str, Any]) -> tuple[int, int]:
    if profile == "production":
        return int(frozen["curve"]["subgroup_order"]), int(frozen["curve"]["cofactor"])
    if profile == "smoke":
        return 71, 2
    raise EvidenceError("unknown relation arithmetic profile")


def target_column_certificate(
    result: dict[str, Any], profile: str, frozen: dict[str, Any]
) -> dict[str, Any]:
    report = result.get("report", {})
    relations = result.get("relation_matrix", [])
    columns = report.get("matrix_columns")
    if type(columns) is not int or columns < 1 or not isinstance(relations, list):
        raise EvidenceError("IC relation matrix has invalid dimensions")
    modulus, cofactor = relation_arithmetic(profile, frozen)
    orbit_columns = columns - 1
    augmented: list[list[int]] = []
    for index, relation in enumerate(relations):
        if not isinstance(relation, dict):
            raise EvidenceError(f"IC relation {index} is not an object")
        row = relation.get("row")
        if not isinstance(row, list) or len(row) != orbit_columns:
            raise EvidenceError(f"IC relation {index} has the wrong orbit-column count")
        coefficients = [
            canonical_residue(value, modulus, f"IC relation {index} column {column}")
            for column, value in enumerate(row)
        ]
        coefficient_a = canonical_residue(
            relation.get("coefficient_a"), modulus, f"IC relation {index} coefficient_a"
        )
        coefficient_b = canonical_residue(
            relation.get("coefficient_b"), modulus, f"IC relation {index} coefficient_b"
        )
        coefficients.append((-cofactor * coefficient_b) % modulus)
        augmented.append(coefficients + [(cofactor * coefficient_a) % modulus])
    pivot_columns: list[int] = []
    pivot_row = 0
    for column in range(columns):
        selected = next(
            (row for row in range(pivot_row, len(augmented)) if augmented[row][column] % modulus),
            None,
        )
        if selected is None:
            continue
        augmented[pivot_row], augmented[selected] = augmented[selected], augmented[pivot_row]
        try:
            inverse = pow(augmented[pivot_row][column], -1, modulus)
        except ValueError as error:
            raise EvidenceError("IC relation pivot is not invertible") from error
        augmented[pivot_row] = [value * inverse % modulus for value in augmented[pivot_row]]
        for row in range(len(augmented)):
            if row == pivot_row or augmented[row][column] == 0:
                continue
            factor = augmented[row][column]
            augmented[row] = [
                (augmented[row][entry] - factor * augmented[pivot_row][entry]) % modulus
                for entry in range(columns + 1)
            ]
        pivot_columns.append(column)
        pivot_row += 1
        if pivot_row == len(augmented):
            break
    consistent = not any(
        all(value == 0 for value in row[:columns]) and row[columns] != 0 for row in augmented
    )
    rank = len(pivot_columns)
    target_column = columns - 1
    target_pivot = target_column in pivot_columns
    free_columns = [column for column in range(columns) if column not in pivot_columns]
    target_row = pivot_columns.index(target_column) if target_pivot else None
    target_invariant = bool(
        consistent and target_row is not None
        and all(augmented[target_row][column] == 0 for column in free_columns)
    )
    target_scalar = augmented[target_row][columns] if target_invariant else None
    return {
        "rows": len(relations), "columns": columns, "rank": rank,
        "nullity": columns - rank, "consistent": consistent,
        "target_column": target_column, "target_pivot": target_pivot,
        "free_columns": free_columns, "target_invariant": target_invariant,
        "target_scalar": str(target_scalar) if target_scalar is not None else None,
    }


def validate_ic(
    result: dict[str, Any], target: dict[str, Any], profile: str, frozen: dict[str, Any]
) -> bool:
    if result.get("schema") != "koblitz_unknown_scalar_ic_result.v1":
        raise EvidenceError("wrong IC result schema")
    if result.get("status") not in {"complete_verified", "incomplete"}:
        raise EvidenceError("unknown IC terminal status")
    if result.get("profile") != profile or profile not in {"production", "smoke"}:
        raise EvidenceError("IC profile mismatch")
    if result.get("target_id") != target["target_id"] or result.get("target") != target["point"]:
        raise EvidenceError("IC target mismatch")
    if result.get("seed") != target.get("ic_seed"):
        raise EvidenceError("IC seed differs from the target panel")
    if result.get("target_scalar_constructed_or_supplied") is not False or result.get("factor_base_logs_constructed_or_supplied") is not False:
        raise EvidenceError("IC process used forbidden scalar labels")
    report = result.get("report", {})
    attempts = result.get("attempt_records", [])
    matrix = result.get("relation_matrix", [])
    if report.get("trials") != len(attempts) or report.get("relations") != len(matrix):
        raise EvidenceError("IC attempt or matrix inventory is incomplete")
    outcomes = report.get("outcomes", {})
    if not isinstance(outcomes, dict) or sum(outcomes.values()) != len(attempts):
        raise EvidenceError("IC outcome partition is incomplete")
    if [row.get("trial") for row in attempts] != list(range(1, len(attempts) + 1)):
        raise EvidenceError("IC attempt sequence is not exact")
    counted = {name: 0 for name in outcomes}
    for row in attempts:
        disposition = row.get("disposition")
        if disposition not in counted:
            raise EvidenceError("IC attempt has an unknown disposition")
        counted[disposition] += 1
    if counted != outcomes or counted.get("relation_found") != len(matrix):
        raise EvidenceError("IC outcome counts differ from retained evidence")
    if sum(row.get("conflicts", 0) for row in attempts) != report.get("sat_conflicts"):
        raise EvidenceError("IC conflict total differs from attempts")
    if sum(row.get("solver_calls", 0) for row in attempts) != report.get("sat_calls"):
        raise EvidenceError("IC solver-call total differs from attempts")
    if sum(row.get("models", 0) for row in attempts) != report.get("sat_models"):
        raise EvidenceError("IC model total differs from attempts")
    if json_blake3(attempts) != result.get("attempt_records_blake3") or json_blake3(matrix) != result.get("relation_matrix_blake3"):
        raise EvidenceError("IC retained transcript hash is invalid")
    progress = result.get("progress", [])
    if json_blake3(progress) != result.get("progress_blake3"):
        raise EvidenceError("IC progress hash is invalid")
    attempt_events = [row for row in progress if row.get("event") == "relation_attempt_finished"]
    if len(attempt_events) != len(attempts) or [row.get("trial") for row in attempt_events] != list(range(1, len(attempts) + 1)):
        raise EvidenceError("IC progress omits or reorders attempts")
    rank_history = report.get("rank_history", [])
    rank_events = [row for row in progress if row.get("event") == "matrix_rank"]
    if (
        not rank_history or report.get("rank_checks") != len(rank_history)
        or report.get("linear_solve_attempts") != len(rank_history)
        or len(rank_events) != len(rank_history) or report.get("matrix_rows") != len(matrix)
        or type(report.get("matrix_columns")) is not int or type(report.get("orbit_count")) is not int
        or report.get("matrix_columns") != report.get("orbit_count") + 1
        or type(report.get("terminal_matrix_rank")) is not int
        or not 0 <= report["terminal_matrix_rank"] <= min(report["matrix_rows"], report["matrix_columns"])
    ):
        raise EvidenceError("IC rank history or progress is incomplete")
    terminal = rank_history[-1]
    terminal_event = rank_events[-1]
    for field in ("rows", "columns", "rank"):
        report_field = "terminal_matrix_rank" if field == "rank" else f"matrix_{field}"
        if terminal.get(field) != report.get(report_field) or terminal_event.get(field) != report.get(report_field):
            raise EvidenceError("IC terminal rank evidence is inconsistent")
    certificate = target_column_certificate(result, profile, frozen)
    if (
        certificate["rows"] != report.get("matrix_rows")
        or certificate["columns"] != report.get("matrix_columns")
        or certificate["rank"] != report.get("terminal_matrix_rank")
    ):
        raise EvidenceError("IC matrix does not reproduce its rank")
    if result["status"] == "complete_verified":
        if not report.get("recovered_scalar_point_verified") or report.get("recovered_scalar") is None:
            raise EvidenceError("completed IC row lacks project point verification")
        if report.get("direct_relation") or report.get("sat_invalid_models") != 0:
            raise EvidenceError("completed IC row bypassed the matrix or admitted invalid models")
        if report["rank_history"][-1].get("candidate_produced") is not True or report["rank_history"][-1].get("candidate_verified") is not True:
            raise EvidenceError("completed IC row lacks candidate history")
        recovered = canonical_residue(
            report["recovered_scalar"], relation_arithmetic(profile, frozen)[0], "IC recovered scalar"
        )
        if not certificate["consistent"] or not certificate["target_invariant"] or certificate["target_scalar"] != str(recovered):
            raise EvidenceError("completed IC row does not identify its reported target scalar")
        return True
    if report.get("recovered_scalar_point_verified") or report.get("recovered_scalar") is not None:
        raise EvidenceError("incomplete IC row retained a completed scalar claim")
    return False


def validate_rho(result: dict[str, Any], target: dict[str, Any], profile: str) -> bool:
    if result.get("schema") != "koblitz_unknown_scalar_rho_result.v1":
        raise EvidenceError("wrong rho result schema")
    if result.get("status") not in {"complete_verified", "incomplete"}:
        raise EvidenceError("unknown rho terminal status")
    if result.get("profile") != profile or profile not in {"production", "smoke"}:
        raise EvidenceError("rho profile mismatch")
    if result.get("target_id") != target["target_id"] or result.get("target") != target["point"]:
        raise EvidenceError("rho target mismatch")
    if result.get("seed") != target.get("rho_seed"):
        raise EvidenceError("rho seed differs from the target panel")
    if result.get("target_scalar_constructed_or_supplied") is not False:
        raise EvidenceError("rho process received a scalar label")
    report = result.get("report", {})
    charges = result.get("charges", {})
    if report.get("jump_table_rebuilds") != report.get("restarts_attempted"):
        raise EvidenceError("rho restart and jump-table counts differ")
    if charges.get("setup_scalar_multiplications") != 34 * report.get("jump_table_rebuilds", -1):
        raise EvidenceError("rho setup scalar-multiplication ledger does not balance")
    if charges.get("setup_group_additions") != 17 * report.get("jump_table_rebuilds", -1):
        raise EvidenceError("rho setup addition ledger does not balance")
    if charges.get("coefficient_draws") != charges.get("setup_scalar_multiplications"):
        raise EvidenceError("rho coefficient-draw ledger does not balance")
    if charges.get("walk_group_additions") != charges.get("partition_hashes") or charges.get("walk_group_additions") != 3 * report.get("iterations", -1):
        raise EvidenceError("rho walk ledger does not balance")
    if charges.get("canonicalizations") != report.get("restarts_attempted", -1) + charges.get("walk_group_additions", -1):
        raise EvidenceError("rho canonicalization ledger does not balance")
    curve_n = 23 if profile == "production" else 7
    if (
        charges.get("frobenius_maps") != charges.get("negations_examined")
        or charges.get("frobenius_maps", 0) > charges.get("canonicalizations", 0) * curve_n
        or charges.get("frobenius_maps", 0) % curve_n != 0
    ):
        raise EvidenceError("rho automorphism ledger does not balance")
    if charges.get("collisions", 0) < charges.get("failed_collisions", 0):
        raise EvidenceError("rho collision ledger does not balance")
    progress = result.get("progress", [])
    if json_blake3(progress) != result.get("progress_blake3"):
        raise EvidenceError("rho progress hash is invalid")
    if sum(row.get("event") == "rho_restart_started" for row in progress) != report.get("restarts_attempted"):
        raise EvidenceError("rho progress omits restarts")
    if sum(row.get("event") == "rho_jump_table_ready" for row in progress) != report.get("jump_table_rebuilds"):
        raise EvidenceError("rho progress omits jump tables")
    if sum(row.get("event") == "rho_collision" for row in progress) != charges.get("collisions"):
        raise EvidenceError("rho progress omits collisions")
    if not progress or progress[-1].get("event") != "rho_finished":
        raise EvidenceError("rho progress lacks its terminal event")
    timing = result.get("timing_ns", {})
    timing_fields = (
        "target_and_subgroup_validation", "rho_setup", "rho_walk",
        "candidate_verification", "end_to_end",
    )
    if any(type(timing.get(field)) is not int or timing[field] < 0 for field in timing_fields):
        raise EvidenceError("rho timing ledger is incomplete")
    if sum(timing[field] for field in timing_fields[:-1]) > timing["end_to_end"]:
        raise EvidenceError("rho timing stages overlap or exceed end-to-end")
    if result["status"] == "complete_verified":
        if not report.get("recovered_scalar_point_verified") or report.get("recovered_scalar") is None or report.get("exhausted") is not False:
            raise EvidenceError("completed rho row lacks project point verification")
        return True
    if report.get("recovered_scalar_point_verified") or report.get("recovered_scalar") is not None or report.get("exhausted") is not True:
        raise EvidenceError("incomplete rho row retained a completed scalar claim")
    return False


def task_resource(task: dict[str, Any]) -> dict[str, Any]:
    process = task["receipt"]["process"]
    return {
        "returncode": process["returncode"],
        "timed_out": process["timed_out"],
        "orphan_group_terminated": process["orphan_group_terminated"],
        **process["metrics"],
    }


def aggregate_task_resources(tasks: list[dict[str, Any]]) -> dict[str, Any]:
    resources = [task_resource(task) for task in tasks]
    return {
        "processes": len(resources),
        "summed_user_seconds": sum(row["user_seconds"] for row in resources),
        "summed_system_seconds": sum(row["system_seconds"] for row in resources),
        "total_core_seconds": sum(row["total_core_seconds"] for row in resources),
        "summed_process_wall_seconds": sum(row["wall_seconds"] for row in resources),
        "peak_process_rss_bytes": max((row["peak_rss_bytes"] for row in resources), default=0),
        "single_core_elapsed_seconds": None,
        "single_core_seconds_legacy_alias": "total_core_seconds",
        "aggregate_parallel_rss_bytes": None,
    }


def reconstruct_row(
    index: int,
    target: dict[str, Any],
    ic_task: dict[str, Any],
    rho_task: dict[str, Any],
    profile: str,
    frozen: dict[str, Any],
) -> dict[str, Any]:
    ic_error = None
    rho_error = None
    try:
        ic_complete = bool(ic_task["result"] and validate_ic(ic_task["result"], target, profile, frozen))
    except EvidenceError as error:
        ic_complete = False
        ic_error = str(error)
    try:
        rho_complete = bool(rho_task["result"] and validate_rho(rho_task["result"], target, profile))
    except EvidenceError as error:
        rho_complete = False
        rho_error = str(error)
    recovered_match = bool(
        ic_complete and rho_complete
        and ic_task["result"]["report"]["recovered_scalar"] == rho_task["result"]["report"]["recovered_scalar"]
    )
    return {
        "index": index, "target_id": target["target_id"], "target": target["point"],
        "ic": {
            "complete": ic_complete,
            "status": ic_task["result"].get("status") if ic_task["result"] else "process_incomplete",
            "resources": task_resource(ic_task), "result": ic_task["receipt"]["result"],
            "validation_error": ic_error,
        },
        "rho": {
            "complete": rho_complete,
            "status": rho_task["result"].get("status") if rho_task["result"] else "process_incomplete",
            "resources": task_resource(rho_task), "result": rho_task["receipt"]["result"],
            "validation_error": rho_error,
        },
        "recovered_scalars_match": recovered_match,
        "row_complete": ic_complete and rho_complete and recovered_match,
        "incomplete_is_inconclusive": not (ic_complete and rho_complete and recovered_match),
    }


def compute_task_ratios(tasks: list[dict[str, Any]], rows: list[dict[str, Any]]) -> dict[str, Any]:
    all_complete = bool(rows) and all(row["row_complete"] for row in rows)
    ic_core = sum(row["ic"]["resources"]["total_core_seconds"] for row in rows)
    rho_core = sum(row["rho"]["resources"]["total_core_seconds"] for row in rows)
    discovery_core = sum(
        task_resource(task)["total_core_seconds"] for task in tasks
        if "discovery-" in task["receipt"]["name"]
    )
    return {
        "online_ic_over_rho_core": ic_core / rho_core if rho_core else None,
        "setup_charged_ic_over_rho_core": (ic_core + discovery_core) / rho_core if rho_core else None,
        "verdict": "finite_complete_panel" if all_complete else "incomplete_inconclusive",
    }


def _source_tree_listing(records: list[dict[str, Any]]) -> bytes:
    return b"".join(
        f"{record['mode']} blob {record['blob_oid']}\t{record['path']}\n".encode()
        for record in records
    )


def _source_manifest(
    *,
    repo: Path,
    source: dict[str, Any],
    bundle: Path,
    original_root_lock: str,
) -> dict[str, Any]:
    git_state = require_exact(
        source.get("git"),
        {"commit", "branch", "dirty", "porcelain", "committed_tree_sha256"},
        "source Git binding",
    )
    commit = require_hex(git_state["commit"], 40, "source commit")
    records, listing = parse_git_tree(repo, commit)
    if sha256(listing) != git_state["committed_tree_sha256"]:
        raise EvidenceError("source binding committed-tree digest differs from Git")
    tree_oid = run_command(
        ["git", "-C", str(repo), "rev-parse", f"{commit}^{{tree}}"],
        "source root tree",
    ).decode().strip()
    require_hex(tree_oid, 40, "source root tree")
    if compute_tree_oid(records) != tree_oid:
        raise EvidenceError("source root tree cannot be reconstructed from its Git entries")
    commit_bytes = run_command(
        ["git", "-C", str(repo), "cat-file", "commit", commit],
        "source commit object",
    )
    if git_object_oid("commit", commit_bytes) != commit:
        raise EvidenceError("source commit object hash is invalid")
    tree_lines = [line for line in commit_bytes.splitlines() if line.startswith(b"tree ")]
    if tree_lines != [f"tree {tree_oid}".encode()]:
        raise EvidenceError("source commit does not bind the reconstructed root tree")

    completed = materialize_git_tree(repo, records, bundle / "source-tree")
    write_new(bundle / "source-metadata/commit.object", commit_bytes)

    by_path = {record["path"]: record for record in completed}
    direct = source.get("direct_sources")
    if not isinstance(direct, dict) or set(direct) != set(DIRECT_SOURCE_PATHS):
        raise EvidenceError("source binding direct-source closure changed")
    direct_manifest: list[dict[str, Any]] = []
    for relative in DIRECT_SOURCE_PATHS:
        bound = validate_identity_shape(direct[relative], f"direct source {relative}")
        expected_original = str(Path(repo) / relative)
        if bound["path"] != expected_original:
            raise EvidenceError(f"direct source {relative} points outside the typed source root")
        if relative == "Cargo.lock":
            if relative in by_path:
                raise EvidenceError("workspace Cargo.lock must remain ignored, not tracked")
            frozen = by_path.get(FROZEN_LOCK_RELATIVE)
            if frozen is None:
                raise EvidenceError("tracked frozen Cargo.lock is absent from the source tree")
            frozen_bytes = regular_bytes(
                bundle / "source-tree" / FROZEN_LOCK_RELATIVE,
                "tracked frozen Cargo.lock",
            )
            if (len(frozen_bytes), sha256(frozen_bytes)) != (bound["bytes"], bound["sha256"]):
                raise EvidenceError("ignored workspace Cargo.lock differs from the tracked frozen lock")
            write_new(bundle / "source-extra/Cargo.lock", frozen_bytes)
            archive_relative = "source-extra/Cargo.lock"
            tracked = False
            blob_oid = None
            origin = "ignored_workspace_copy_byte_equal_to_tracked_frozen_lock"
        else:
            tracked_record = by_path.get(relative)
            if tracked_record is None:
                raise EvidenceError(f"direct source {relative} is not tracked at the bound commit")
            if (tracked_record["bytes"], tracked_record["sha256"]) != (bound["bytes"], bound["sha256"]):
                raise EvidenceError(f"direct source {relative} differs from its bound commit")
            archive_relative = f"source-tree/{relative}"
            tracked = True
            blob_oid = tracked_record["blob_oid"]
            origin = "bound_commit_tree"
        direct_manifest.append(
            {
                "relative": relative,
                "original_identity": bound,
                "archive_relative": archive_relative,
                "tracked_at_bound_commit": tracked,
                "git_blob_oid": blob_oid,
                "origin": origin,
            }
        )

    dependency_lock = validate_identity_shape(source.get("dependency_lock"), "dependency lock")
    if dependency_lock != direct[FROZEN_LOCK_RELATIVE]:
        raise EvidenceError("source dependency-lock binding changed")
    return {
        "closure_kind": "complete_committed_git_tree_plus_ignored_workspace_lock",
        "commit": commit,
        "commit_object": identity(
            bundle / "source-metadata/commit.object",
            recorded_path="source-metadata/commit.object",
            context="source commit object",
        ),
        "root_tree_oid": tree_oid,
        "committed_tree_sha256": git_state["committed_tree_sha256"],
        "tracked_files": completed,
        "direct_sources": direct_manifest,
        "workspace_cargo_lock": {
            "original_absolute": original_root_lock,
            "archive_relative": "source-extra/Cargo.lock",
            "tracked_at_bound_commit": False,
            "ignored_workspace_file": True,
            "byte_equal_to": FROZEN_LOCK_RELATIVE,
            "identity": identity(
                bundle / "source-extra/Cargo.lock",
                recorded_path="source-extra/Cargo.lock",
                context="ignored workspace Cargo.lock",
            ),
        },
        "external_dependency_sources_archived": False,
        "external_dependency_binding": "tracked frozen Cargo.lock",
    }


def _validate_terminal_run_input(run_root: Path) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    seal, _ = read_json(run_root / "run-seal.json", "original run seal")
    require_exact(
        seal,
        {
            "schema", "status", "profile", "protocol_sha256", "source_commit",
            "summary", "inventory", "inventory_sha256", "expected_inner_command",
            "panel_complete", "scientific_measurement_admitted", "seal_payload_sha256",
        },
        "original run seal",
    )
    if seal["schema"] != RUN_SEAL_SCHEMA or seal["status"] != "outputs_frozen":
        raise EvidenceError("original run is not terminal")
    validate_self_hash(seal, "seal_payload_sha256", "original run seal")
    inventory = validate_inventory_records(seal["inventory"], "original run-seal inventory")
    if canonical_sha256(inventory) != seal["inventory_sha256"]:
        raise EvidenceError("original run-seal inventory hash is invalid")
    actual = tree_inventory(run_root, excluded={"run-seal.json"})
    if actual != inventory:
        raise EvidenceError("original run bytes changed after terminal sealing")
    if seal["scientific_measurement_admitted"] is not False:
        raise EvidenceError("original run widened scientific admission")
    return seal, inventory


def package_bundle(
    *,
    run_root: Path,
    outer_root: Path,
    outer_metrics: Path,
    project_verification_root: Path,
    output: Path,
) -> dict[str, Any]:
    run_root = real_directory(run_root, "original terminal run")
    outer_root = real_directory(outer_root, "original outer receipt root")
    project_verification_root = real_directory(
        project_verification_root, "original project verification root"
    )
    output = safe_new_output(output)
    outer_metrics = Path(outer_metrics)
    if not outer_metrics.is_absolute():
        outer_metrics = Path.cwd() / outer_metrics
    normalized_absolute(str(outer_metrics), "outer metrics")
    try:
        outer_metrics_relative = outer_metrics.relative_to(outer_root).as_posix()
    except ValueError as error:
        raise EvidenceError("outer metrics must be inside the typed outer root") from error
    safe_relative(outer_metrics_relative, "outer metrics relative path")

    seal, run_inventory = _validate_terminal_run_input(run_root)
    if any(record["path"] == "build-target" for record in run_inventory):
        raise EvidenceError("build-target must be a directory, never a sealed file")
    retained = [record for record in run_inventory if not is_build_target_path(record["path"])]
    omitted = [record for record in run_inventory if is_build_target_path(record["path"])]
    if retained + omitted == run_inventory:
        # The source inventory is globally sorted, so interleaving is possible.  The
        # partition itself remains separately sorted and is checked by the verifier.
        pass
    if len(retained) + len(omitted) != len(run_inventory):
        raise EvidenceError("run inventory partition is not exhaustive")
    if any(record["path"] == "build-target" for record in retained):
        raise EvidenceError("build-target root was retained unexpectedly")

    summary, _ = read_json(run_root / "run-summary.json", "original run summary")
    source, _ = read_json(run_root / "inputs/source.json", "original source binding")
    source = require_exact(
        source,
        {"schema", "git", "direct_sources", "dependency_lock", "python", "cargo", "rustc"},
        "original source binding",
    )
    if source["schema"] != "koblitz_unknown_scalar_source_binding.v1":
        raise EvidenceError("original source binding schema changed")
    direct = source.get("direct_sources", {})
    cargo_toml = validate_identity_shape(direct.get("Cargo.toml"), "bound Cargo.toml")
    source_repo = real_directory(Path(cargo_toml["path"]).parent, "bound source repository")
    if cargo_toml["path"] != str(source_repo / "Cargo.toml"):
        raise EvidenceError("bound source repository path is not normalized")
    root_lock_original = validate_identity_shape(
        direct.get("Cargo.lock"), "bound ignored workspace Cargo.lock"
    )["path"]

    outer_inventory = tree_inventory(outer_root)
    project_inventory = tree_inventory(project_verification_root)
    if {record["path"] for record in outer_inventory} != {
        "driver.metrics.json", "driver.stdout", "driver.stderr"
    }:
        raise EvidenceError("original outer tree violates the exact driver receipt grammar")
    if outer_metrics_relative != "driver.metrics.json":
        raise EvidenceError("outer metrics must be driver.metrics.json")
    if {record["path"] for record in project_inventory} != {
        "verification.json", "verification-seal.json"
    }:
        raise EvidenceError("original project-verification tree violates the exact two-file grammar")
    verification_seal, _ = read_json(
        project_verification_root / "verification-seal.json", "original project verification seal"
    )
    if (
        verification_seal.get("schema") != PROJECT_VERIFICATION_SEAL_SCHEMA
        or verification_seal.get("status") != "verification_frozen"
    ):
        raise EvidenceError("original project verification is not terminal")
    validate_self_hash(
        verification_seal, "seal_payload_sha256", "original project verification seal"
    )
    frozen_project_inventory = validate_inventory_records(
        verification_seal.get("inventory"), "project verification inventory"
    )
    if (
        project_inventory
        != sorted(
            [
                *frozen_project_inventory,
                identity(
                    project_verification_root / "verification-seal.json",
                    recorded_path="verification-seal.json",
                    context="project verification seal",
                ),
            ],
            key=lambda item: item["path"],
        )
        or canonical_sha256(frozen_project_inventory) != verification_seal.get("inventory_sha256")
    ):
        raise EvidenceError("original project verification inventory changed after sealing")

    create_new_directory(output, "bundle output")
    copy_tree_bytes(run_root, output / "original/run", retained)
    write_new(
        output / "original/run/run-seal.json",
        regular_bytes(run_root / "run-seal.json", "original run seal"),
    )
    copy_tree_bytes(outer_root, output / "original/outer", outer_inventory)
    copy_tree_bytes(
        project_verification_root,
        output / "original/project-verification",
        project_inventory,
    )

    core_source = Path(__file__).resolve()
    package_source = core_source.with_name("package_koblitz_stage23_terminal_evidence.py")
    verifier_source = core_source.with_name("verify_koblitz_stage23_terminal_evidence.py")
    for source_path in (core_source, package_source, verifier_source):
        write_new(
            output / "verifier" / source_path.name,
            regular_bytes(source_path, f"portable tool source {source_path.name}"),
            0o644,
        )

    source_manifest = _source_manifest(
        repo=source_repo,
        source=source,
        bundle=output,
        original_root_lock=root_lock_original,
    )

    by_run_path = {record["path"]: record for record in run_inventory}
    binary_bindings: list[dict[str, Any]] = []
    claimed_binaries = summary.get("binaries")
    if not isinstance(claimed_binaries, dict):
        raise EvidenceError("run summary binary inventory is malformed")
    for name in sorted(claimed_binaries):
        if name not in BINARY_NAMES:
            raise EvidenceError(f"unexpected retained Stage-23 binary {name}")
        retained_path = f"binaries/{name}"
        omitted_path = f"build-target/release/examples/{name}"
        retained_record = by_run_path.get(retained_path)
        omitted_record = by_run_path.get(omitted_path)
        claimed = validate_identity_shape(claimed_binaries[name], f"summary binary {name}")
        if retained_record is None or omitted_record is None:
            raise EvidenceError(f"binary {name} lacks retained/omitted sealed counterparts")
        if claimed["path"] != str(run_root / retained_path):
            raise EvidenceError(f"summary binary {name} points outside the typed run root")
        if (claimed["bytes"], claimed["sha256"]) != (
            retained_record["bytes"], retained_record["sha256"]
        ):
            raise EvidenceError(f"summary binary {name} differs from the retained sealed file")
        if (retained_record["bytes"], retained_record["sha256"]) != (
            omitted_record["bytes"], omitted_record["sha256"]
        ):
            raise EvidenceError(f"retained binary {name} differs from its omitted fresh build output")
        binary_bindings.append(
            {
                "name": name,
                "retained_run_path": retained_path,
                "retained_identity": retained_record,
                "omitted_build_path": omitted_path,
                "omitted_identity": omitted_record,
                "bytes_equal_by_original_run_seal": True,
            }
        )
    if claimed_binaries and set(claimed_binaries) != set(BINARY_NAMES):
        raise EvidenceError("complete binary inventory must contain both Stage-23 executables")

    roots = [
        {
            "role": role,
            "kind": kind,
            "original_absolute": {
                "run": str(run_root),
                "outer": str(outer_root),
                "project_verification": str(project_verification_root),
                "source_repository": str(source_repo),
            }[role],
            "archive_relative": archive,
        }
        for role, (kind, archive) in ROOT_LAYOUT.items()
    ]
    path_map = {
        "schema": PATH_MAP_SCHEMA,
        "roots": roots,
        "files": [
            {
                "role": "ignored_workspace_lock",
                "kind": "ignored_workspace_dependency_lock",
                "original_absolute": root_lock_original,
                "archive_relative": "source-extra/Cargo.lock",
            },
            {
                "role": "source_commit_object",
                "kind": "git_commit_object",
                "original_absolute": str(source_repo / ".git-objects" / seal["source_commit"]),
                "archive_relative": "source-metadata/commit.object",
            },
            {
                "role": "packager_source",
                "kind": "portable_tool_source",
                "original_absolute": str(package_source),
                "archive_relative": "verifier/package_koblitz_stage23_terminal_evidence.py",
            },
            {
                "role": "core_source",
                "kind": "portable_tool_source",
                "original_absolute": str(core_source),
                "archive_relative": "verifier/koblitz_stage23_terminal_evidence.py",
            },
            {
                "role": "verifier_source",
                "kind": "portable_tool_source",
                "original_absolute": str(verifier_source),
                "archive_relative": "verifier/verify_koblitz_stage23_terminal_evidence.py",
            },
        ],
        "external_tools": [
            {
                "role": role,
                "kind": "external_tool_identity",
                "original_absolute": source[role]["path"],
                "archived": False,
                "identity": source[role],
            }
            for role in ("python", "cargo", "rustc")
        ],
    }

    manifest = {
        "schema": MANIFEST_SCHEMA,
        "status": "compact_terminal_evidence_packaged",
        "profile": seal["profile"],
        "path_map": path_map,
        "run_partition": {
            "source_inventory_sha256": seal["inventory_sha256"],
            "omission_rule": "reject a build-target file and omit exactly build-target descendants",
            "retained": retained,
            "omitted": omitted,
            "run_seal": identity(
                output / "original/run/run-seal.json",
                recorded_path="original/run/run-seal.json",
                context="archived original run seal",
            ),
        },
        "binary_build_cross_bindings": binary_bindings,
        "outer": {
            "inventory": outer_inventory,
            "metrics_relative": outer_metrics_relative,
            "metrics_original_identity": identity(
                outer_metrics,
                recorded_path=str(outer_metrics),
                context="original outer metrics",
            ),
        },
        "project_verification": {
            "inventory": project_inventory,
            "seal_relative": "verification-seal.json",
        },
        "source_closure": source_manifest,
        "portable_tools": {
            "packager": identity(
                output / "verifier/package_koblitz_stage23_terminal_evidence.py",
                recorded_path="verifier/package_koblitz_stage23_terminal_evidence.py",
                context="archived packager",
            ),
            "core": identity(
                output / "verifier/koblitz_stage23_terminal_evidence.py",
                recorded_path="verifier/koblitz_stage23_terminal_evidence.py",
                context="archived portable core",
            ),
            "verifier": identity(
                output / "verifier/verify_koblitz_stage23_terminal_evidence.py",
                recorded_path="verifier/verify_koblitz_stage23_terminal_evidence.py",
                context="archived verifier",
            ),
        },
        "claim_boundary": {
            "project_custody_reconstruction_only": True,
            "independent_mathematical_payload_replay_completed": False,
            "pending_independent_mathematical_payload_replay": PENDING_MATH_REPLAY,
            "scientific_measurement_admitted": False,
            "external_portable_verification_satisfied": False,
            "independent_external_reproduction_satisfied": False,
            "full_cost_gate_passed": False,
            "koblitz_index_calculus_sota": False,
        },
    }
    write_json_new(output / "bundle-manifest.json", manifest)
    bundle_inventory = tree_inventory(output, excluded={"bundle-seal.json"})
    seal_payload = {
        "schema": SEAL_SCHEMA,
        "status": "compact_terminal_evidence_frozen",
        "manifest": identity(
            output / "bundle-manifest.json",
            recorded_path="bundle-manifest.json",
            context="bundle manifest",
        ),
        "inventory": bundle_inventory,
        "inventory_sha256": canonical_sha256(bundle_inventory),
        "source_run_inventory_sha256": seal["inventory_sha256"],
        "scientific_measurement_admitted": False,
        "external_portable_verification_satisfied": False,
        "independent_external_reproduction_satisfied": False,
        "full_cost_gate_passed": False,
        "koblitz_index_calculus_sota": False,
    }
    bundle_seal = dict(seal_payload)
    bundle_seal["seal_payload_sha256"] = canonical_sha256(seal_payload)
    write_json_new(output / "bundle-seal.json", bundle_seal)
    verification = verify_bundle(output)
    return {
        "schema": "koblitz_stage23_terminal_evidence_package_result.v1",
        "status": "package_frozen_and_verified",
        "bundle": str(output),
        "bundle_seal": identity(output / "bundle-seal.json", context="bundle seal"),
        "verification": verification,
    }


def _verify_source_closure(
    bundle: Path,
    manifest: dict[str, Any],
    path_map: ArchiveMap,
    source_binding: dict[str, Any],
    run_seal: dict[str, Any],
) -> None:
    source = require_exact(
        manifest.get("source_closure"),
        {
            "closure_kind", "commit", "commit_object", "root_tree_oid",
            "committed_tree_sha256", "tracked_files", "direct_sources",
            "workspace_cargo_lock", "external_dependency_sources_archived",
            "external_dependency_binding",
        },
        "source closure",
    )
    if source["closure_kind"] != "complete_committed_git_tree_plus_ignored_workspace_lock":
        raise EvidenceError("source closure kind changed")
    commit = require_hex(source["commit"], 40, "source closure commit")
    if commit != run_seal["source_commit"] or commit != source_binding["git"]["commit"]:
        raise EvidenceError("source closure commit differs from the sealed run")
    tree_oid = require_hex(source["root_tree_oid"], 40, "source closure root tree")

    tracked = source["tracked_files"]
    if not isinstance(tracked, list) or not tracked:
        raise EvidenceError("source closure has no tracked files")
    checked: list[dict[str, Any]] = []
    seen: set[str] = set()
    for index, record in enumerate(tracked):
        record = require_exact(
            record, {"path", "mode", "bytes", "blob_oid", "sha256"},
            f"source tracked file {index}",
        )
        relative = safe_relative(record["path"], f"source tracked file {index}.path")
        if relative in seen:
            raise EvidenceError(f"source closure contains duplicate tracked path {relative}")
        seen.add(relative)
        if record["mode"] not in {"100644", "100755"}:
            raise EvidenceError(f"source tracked file {relative} has an unsupported mode")
        if isinstance(record["bytes"], bool) or not isinstance(record["bytes"], int) or record["bytes"] < 0:
            raise EvidenceError(f"source tracked file {relative} has an invalid size")
        require_hex(record["blob_oid"], 40, f"source tracked file {relative} blob")
        require_hex(record["sha256"], 64, f"source tracked file {relative} sha256")
        data = regular_bytes(path_map.archived("source_repository", relative), f"source tracked file {relative}")
        if (
            len(data) != record["bytes"]
            or sha256(data) != record["sha256"]
            or git_object_oid("blob", data) != record["blob_oid"]
        ):
            raise EvidenceError(f"source tracked file {relative} changed")
        checked.append(record)
    if checked != sorted(checked, key=lambda record: record["path"]):
        raise EvidenceError("source tracked-file inventory is not sorted")
    actual_source_paths = [record["path"] for record in tree_inventory(path_map.archived("source_repository"))]
    if actual_source_paths != [record["path"] for record in checked]:
        raise EvidenceError("source tree contains missing or extra files")
    if compute_tree_oid(checked) != tree_oid:
        raise EvidenceError("archived source files do not reconstruct the bound Git tree")
    listing = _source_tree_listing(checked)
    if (
        sha256(listing) != source["committed_tree_sha256"]
        or source["committed_tree_sha256"] != source_binding["git"]["committed_tree_sha256"]
    ):
        raise EvidenceError("archived source listing differs from the source binding")

    commit_identity = validate_identity_shape(source["commit_object"], "source commit object identity")
    if commit_identity["path"] != "source-metadata/commit.object":
        raise EvidenceError("source commit object has the wrong archive path")
    commit_path = path_map.archived_file("source_commit_object")
    commit_bytes = regular_bytes(commit_path, "source commit object")
    if identity(commit_path, recorded_path=commit_identity["path"]) != commit_identity:
        raise EvidenceError("source commit object identity changed")
    if git_object_oid("commit", commit_bytes) != commit:
        raise EvidenceError("archived source commit object hash is invalid")
    if [line for line in commit_bytes.splitlines() if line.startswith(b"tree ")] != [
        f"tree {tree_oid}".encode()
    ]:
        raise EvidenceError("archived source commit does not point to the reconstructed tree")

    direct_manifest = source["direct_sources"]
    if not isinstance(direct_manifest, list) or len(direct_manifest) != len(DIRECT_SOURCE_PATHS):
        raise EvidenceError("source closure direct-source manifest changed")
    direct_by_path: dict[str, dict[str, Any]] = {}
    for entry in direct_manifest:
        entry = require_exact(
            entry,
            {
                "relative", "original_identity", "archive_relative",
                "tracked_at_bound_commit", "git_blob_oid", "origin",
            },
            "source direct-source entry",
        )
        relative = safe_relative(entry["relative"], "source direct-source relative path")
        if relative in direct_by_path:
            raise EvidenceError(f"source closure duplicates direct source {relative}")
        direct_by_path[relative] = entry
    if set(direct_by_path) != set(DIRECT_SOURCE_PATHS):
        raise EvidenceError("source closure direct-source set changed")
    if not isinstance(source_binding.get("direct_sources"), dict) or set(source_binding["direct_sources"]) != set(DIRECT_SOURCE_PATHS):
        raise EvidenceError("archived source binding direct-source set changed")
    tracked_by_path = {record["path"]: record for record in checked}
    for relative in DIRECT_SOURCE_PATHS:
        entry = direct_by_path[relative]
        original = validate_identity_shape(entry["original_identity"], f"direct source {relative}")
        if original != source_binding["direct_sources"][relative]:
            raise EvidenceError(f"source manifest and binding differ for {relative}")
        if original["path"] != path_map.original("source_repository", relative):
            raise EvidenceError(f"direct source {relative} escapes the typed source root")
        if relative == "Cargo.lock":
            if (
                entry["archive_relative"] != "source-extra/Cargo.lock"
                or entry["tracked_at_bound_commit"] is not False
                or entry["git_blob_oid"] is not None
                or entry["origin"] != "ignored_workspace_copy_byte_equal_to_tracked_frozen_lock"
                or relative in tracked_by_path
            ):
                raise EvidenceError("workspace Cargo.lock was incorrectly claimed as tracked")
            archive_path = path_map.archived_file("ignored_workspace_lock")
        else:
            record = tracked_by_path.get(relative)
            if (
                record is None
                or entry["archive_relative"] != f"source-tree/{relative}"
                or entry["tracked_at_bound_commit"] is not True
                or entry["git_blob_oid"] != record["blob_oid"]
                or entry["origin"] != "bound_commit_tree"
            ):
                raise EvidenceError(f"tracked direct source classification changed for {relative}")
            archive_path = path_map.archived("source_repository", relative)
        data = regular_bytes(archive_path, f"archived direct source {relative}")
        if (len(data), sha256(data)) != (original["bytes"], original["sha256"]):
            raise EvidenceError(f"archived direct source {relative} differs from its original identity")

    workspace_lock = require_exact(
        source["workspace_cargo_lock"],
        {
            "original_absolute", "archive_relative", "tracked_at_bound_commit",
            "ignored_workspace_file", "byte_equal_to", "identity",
        },
        "workspace Cargo.lock classification",
    )
    if (
        workspace_lock["original_absolute"] != path_map.files["ignored_workspace_lock"]["original_absolute"]
        or workspace_lock["archive_relative"] != "source-extra/Cargo.lock"
        or workspace_lock["tracked_at_bound_commit"] is not False
        or workspace_lock["ignored_workspace_file"] is not True
        or workspace_lock["byte_equal_to"] != FROZEN_LOCK_RELATIVE
    ):
        raise EvidenceError("workspace Cargo.lock tracking classification changed")
    lock_identity = validate_identity_shape(workspace_lock["identity"], "workspace Cargo.lock identity")
    lock_path = path_map.archived_file("ignored_workspace_lock")
    if identity(lock_path, recorded_path="source-extra/Cargo.lock") != lock_identity:
        raise EvidenceError("workspace Cargo.lock archive identity changed")
    frozen_bytes = regular_bytes(
        path_map.archived("source_repository", FROZEN_LOCK_RELATIVE),
        "tracked frozen Cargo.lock",
    )
    if regular_bytes(lock_path, "ignored workspace Cargo.lock") != frozen_bytes:
        raise EvidenceError("ignored workspace Cargo.lock is not byte-equal to the tracked frozen lock")
    if source_binding["dependency_lock"] != source_binding["direct_sources"][FROZEN_LOCK_RELATIVE]:
        raise EvidenceError("archived dependency-lock binding changed")
    if source["external_dependency_sources_archived"] is not False or source["external_dependency_binding"] != "tracked frozen Cargo.lock":
        raise EvidenceError("source closure overstated external dependency custody")


def _reconstruct_task(
    *,
    path_map: ArchiveMap,
    name: str,
    command: list[str],
    watchdog: float,
    environment: dict[str, str],
    inputs: list[dict[str, Any]],
    meter_identity: dict[str, Any],
    expect_json: bool,
    run_paths: set[str],
) -> dict[str, Any]:
    base = f"tasks/{name}"
    expected_files = {"intent.json", "stdout", "stderr", "metrics.json", "receipt.json"}
    present = {
        path[len(base) + 1 :]
        for path in run_paths
        if path.startswith(base + "/") and "/" not in path[len(base) + 1 :]
    }
    receipt, _ = read_json(path_map.archived("run", f"{base}/receipt.json"), f"{name} receipt")
    receipt = require_exact(
        receipt,
        {
            "schema", "name", "process", "intent", "stdout", "stderr", "metrics",
            "inputs", "result", "parse_error", "terminal_status",
        },
        f"{name} receipt",
    )
    if receipt["schema"] != "koblitz_unknown_scalar_process_receipt.v1" or receipt["name"] != name:
        raise EvidenceError(f"{name} receipt identity changed")
    for field, relative in (
        ("intent", f"{base}/intent.json"),
        ("stdout", f"{base}/stdout"),
        ("stderr", f"{base}/stderr"),
        ("metrics", f"{base}/metrics.json"),
    ):
        path_map.validate_bound_file(receipt[field], "run", relative, f"{name} {field}")

    intent, _ = read_json(path_map.archived("run", f"{base}/intent.json"), f"{name} intent")
    expected_intent = {
        "schema": "koblitz_unknown_scalar_process_intent.v1",
        "name": name,
        "command": command,
        "watchdog_seconds": watchdog,
        "environment": environment,
        "inputs": inputs,
        "meter": meter_identity,
        "meter_exclusive_create": True,
        "stdin": "devnull",
        "close_fds": True,
    }
    if intent != expected_intent or receipt["inputs"] != inputs:
        raise EvidenceError(f"{name} intent or input graph differs from reconstruction")
    raw_metrics, _ = read_json(path_map.archived("run", f"{base}/metrics.json"), f"{name} metrics")
    raw_metrics = validate_metrics(raw_metrics, f"{name} metrics")
    if receipt["process"] != raw_metrics:
        raise EvidenceError(f"{name} receipt differs from raw metrics")
    if raw_metrics["command"] != command or raw_metrics["watchdog_seconds"] != watchdog:
        raise EvidenceError(f"{name} raw command or watchdog differs from reconstruction")
    clean = process_clean(raw_metrics)
    parsed: dict[str, Any] | None = None
    parse_error: str | None = None
    if expect_json and clean:
        try:
            parsed, _ = read_json(path_map.archived("run", f"{base}/stdout"), f"{name} stdout")
        except Exception as error:
            parse_error = str(error)
    result_relative = f"{base}/result.json"
    if parsed is not None:
        expected_files.add("result.json")
        path_map.validate_bound_file(receipt["result"], "run", result_relative, f"{name} result")
        result, _ = read_json(path_map.archived("run", result_relative), f"{name} result")
        if result != parsed:
            raise EvidenceError(f"{name} result differs from parsed stdout")
    elif receipt["result"] is not None or result_relative in run_paths:
        raise EvidenceError(f"{name} retained a result for an incomplete process")
    terminal = "complete" if clean and (parsed is not None or not expect_json) else "incomplete"
    if receipt["parse_error"] != parse_error or receipt["terminal_status"] != terminal:
        raise EvidenceError(f"{name} terminal status is not derivable")
    if present != expected_files:
        raise EvidenceError(f"{name} task contains missing or extra files")
    nested = [path for path in run_paths if path.startswith(base + "/") and path[len(base) + 1 :].count("/")]
    if nested:
        raise EvidenceError(f"{name} task contains nested extra files")
    return {"receipt": receipt, "result": parsed}


def _verify_run_reconstruction(
    *,
    path_map: ArchiveMap,
    manifest: dict[str, Any],
    run_seal: dict[str, Any],
    retained: list[dict[str, Any]],
) -> dict[str, Any]:
    run_paths = {record["path"] for record in retained}
    summary, _ = read_json(path_map.archived("run", "run-summary.json"), "archived run summary")
    summary = require_exact(
        summary,
        {
            "schema", "status", "profile", "evidence_class", "protocol_sha256",
            "source_revision", "host", "critical_failure", "target_panel", "binaries",
            "rows", "completed_rows", "expected_rows", "resources", "ratios",
            "full_cost_gate_passed", "independent_external_reproduction_satisfied",
            "koblitz_index_calculus_sota", "outer_driver_accounting", "claim_boundary",
        },
        "archived run summary",
    )
    if summary["schema"] != RUN_SUMMARY_SCHEMA or summary["profile"] not in {"production", "smoke"}:
        raise EvidenceError("run summary schema or profile changed")
    if summary["profile"] != run_seal["profile"]:
        raise EvidenceError("run summary and seal profiles differ")
    path_map.validate_bound_file(run_seal["summary"], "run", "run-summary.json", "run summary")

    inputs = {path for path in run_paths if path.startswith("inputs/")}
    if inputs != {"inputs/protocol.json", "inputs/source.json", "inputs/host.json"}:
        raise EvidenceError("run input inventory changed")
    protocol, _ = read_json(path_map.archived("run", "inputs/protocol.json"), "archived protocol")
    source_binding, _ = read_json(path_map.archived("run", "inputs/source.json"), "archived source binding")
    source_binding = require_exact(
        source_binding,
        {"schema", "git", "direct_sources", "dependency_lock", "python", "cargo", "rustc"},
        "archived source binding",
    )
    if source_binding["schema"] != "koblitz_unknown_scalar_source_binding.v1":
        raise EvidenceError("source binding schema changed")
    for role in ("python", "cargo", "rustc"):
        bound_tool = validate_identity_shape(source_binding[role], f"source-bound {role}")
        mapped_tool = path_map.external_tools[role]
        if (
            mapped_tool["identity"] != bound_tool
            or mapped_tool["original_absolute"] != bound_tool["path"]
        ):
            raise EvidenceError(f"external tool map differs from source binding for {role}")
    host, _ = read_json(path_map.archived("run", "inputs/host.json"), "archived host binding")
    host = require_exact(
        host,
        {"schema", "platform", "uname", "logical_cpus", "scientific_environment"},
        "archived host binding",
    )
    if host["schema"] != "koblitz_unknown_scalar_host_binding.v1" or host["scientific_environment"] != EXPECTED_ENVIRONMENT:
        raise EvidenceError("host or scientific environment binding changed")
    if host != summary["host"]:
        raise EvidenceError("summary host differs from the archived host binding")
    if source_binding["git"] != summary["source_revision"] or source_binding["git"]["commit"] != run_seal["source_commit"]:
        raise EvidenceError("summary source revision differs from the archived source binding")
    _verify_source_closure(path_map.bundle, manifest, path_map, source_binding, run_seal)

    protocol_identity = identity(
        path_map.archived("run", "inputs/protocol.json"),
        recorded_path=path_map.original("run", "inputs/protocol.json"),
        context="archived protocol",
    )
    source_identity = identity(
        path_map.archived("run", "inputs/source.json"),
        recorded_path=path_map.original("run", "inputs/source.json"),
        context="archived source binding",
    )
    if protocol_identity["sha256"] != summary["protocol_sha256"] or protocol_identity["sha256"] != run_seal["protocol_sha256"]:
        raise EvidenceError("archived protocol digest differs from the run seal")
    meter_identity = source_binding["direct_sources"].get("scripts/process_meter.py")
    validate_identity_shape(meter_identity, "bound process meter")
    frozen = stage23_protocol(protocol, summary["profile"])
    values = stage23_profile_values(summary["profile"], frozen)
    if summary["expected_rows"] != values["targets"]:
        raise EvidenceError("summary expected-row count is not derivable")

    full_names = ["00-build", "01-discovery-a0", "02-discovery-a1", "03-targets"]
    for index in range(1, values["targets"] + 1):
        full_names.extend(
            [f"{2 * index + 2:02d}-row-{index:02d}-ic", f"{2 * index + 3:02d}-row-{index:02d}-rho"]
        )
    task_names = sorted(
        {
            PurePosixPath(path).parts[1]
            for path in run_paths
            if len(PurePosixPath(path).parts) >= 3 and PurePosixPath(path).parts[0] == "tasks"
        }
    )
    if not task_names or task_names != full_names[: len(task_names)]:
        raise EvidenceError("task graph is not an exact Stage-23 execution prefix")
    tasks: list[dict[str, Any]] = []
    task_by_name: dict[str, dict[str, Any]] = {}

    def add_task(name: str, command: list[str], watchdog: float, task_inputs: list[dict[str, Any]], expect_json: bool) -> dict[str, Any]:
        if name not in task_names:
            raise EvidenceError(f"task graph omitted required task {name}")
        task = _reconstruct_task(
            path_map=path_map,
            name=name,
            command=command,
            watchdog=watchdog,
            environment=EXPECTED_ENVIRONMENT,
            inputs=task_inputs,
            meter_identity=meter_identity,
            expect_json=expect_json,
            run_paths=run_paths,
        )
        tasks.append(task)
        task_by_name[name] = task
        return task

    build = add_task(
        "00-build",
        [
            source_binding["cargo"]["path"], "build", "--release", "--locked", "--jobs", "1",
            "--target-dir", path_map.original("run", "build-target"),
            "--example", "koblitz_public_factor_base_discovery",
            "--example", "koblitz_unknown_scalar_panel",
        ],
        float(frozen["execution"]["build_watchdog_seconds"]),
        [source_identity, protocol_identity],
        False,
    )
    build_complete = build["receipt"]["terminal_status"] == "complete"
    binaries = summary["binaries"]
    if not isinstance(binaries, dict):
        raise EvidenceError("summary binaries must be an object")
    binary_paths = {path for path in run_paths if path.startswith("binaries/")}
    if build_complete:
        if set(binaries) != set(BINARY_NAMES) or binary_paths != {f"binaries/{name}" for name in BINARY_NAMES}:
            raise EvidenceError("completed build binary inventory changed")
        for name in BINARY_NAMES:
            path_map.validate_bound_file(binaries[name], "run", f"binaries/{name}", f"binary {name}")
    elif binaries or binary_paths:
        raise EvidenceError("incomplete build retained executables")

    critical_failure: str | None = None
    targets_result: dict[str, Any] | None = None
    targets: list[dict[str, Any]] = []
    rows: list[dict[str, Any]] = []
    if not build_complete:
        critical_failure = "build"
    else:
        discovery_binary = binaries["koblitz_public_factor_base_discovery"]["path"]
        for curve_a in (0, 1):
            name = f"0{curve_a + 1}-discovery-a{curve_a}"
            task = add_task(
                name,
                [discovery_binary, str(values["n"]), str(curve_a), "2", str(values["dimension"])],
                float(frozen["execution"]["discovery_watchdog_seconds"]),
                [protocol_identity],
                True,
            )
            if task["receipt"]["terminal_status"] != "complete":
                critical_failure = f"discovery-a{curve_a}"
                break
            try:
                validate_discovery_result(task["result"], curve_a, values, frozen)
            except EvidenceError:
                critical_failure = f"discovery-a{curve_a}-invalid"
                break
    if critical_failure is None:
        panel_binary = binaries["koblitz_unknown_scalar_panel"]["path"]
        target_task = add_task(
            "03-targets",
            [panel_binary, "targets", summary["profile"]],
            float(frozen["execution"]["target_generation_watchdog_seconds"]),
            [protocol_identity],
            True,
        )
        if target_task["receipt"]["terminal_status"] != "complete":
            critical_failure = "targets"
        else:
            targets_result = target_task["result"]
            try:
                targets = validate_targets(targets_result, summary["profile"], values["targets"])
            except EvidenceError:
                critical_failure = "targets-invalid"
    if critical_failure is None:
        target_identity = task_by_name["03-targets"]["receipt"]["result"]
        panel_binary = binaries["koblitz_unknown_scalar_panel"]["path"]
        for index, target in enumerate(targets, 1):
            common = [target["target_id"], target["point"]["x"], target["point"]["y"]]
            ic_name = f"{2 * index + 2:02d}-row-{index:02d}-ic"
            rho_name = f"{2 * index + 3:02d}-row-{index:02d}-rho"
            ic_task = add_task(
                ic_name,
                [panel_binary, "ic", summary["profile"], *common, str(target["ic_seed"])],
                float(frozen["execution"]["ic_watchdog_seconds"]),
                [target_identity],
                True,
            )
            rho_task = add_task(
                rho_name,
                [panel_binary, "rho", summary["profile"], *common, str(target["rho_seed"])],
                float(frozen["execution"]["rho_watchdog_seconds"]),
                [target_identity],
                True,
            )
            try:
                rows.append(reconstruct_row(index, target, ic_task, rho_task, summary["profile"], frozen))
            except EvidenceError as error:
                raise EvidenceError(f"row {index} reconstruction failed: {error}") from error

    if task_names != [task["receipt"]["name"] for task in tasks]:
        raise EvidenceError("task graph continued after failure or omitted a required task")
    expected_run_paths = {
        "run-summary.json",
        "inputs/protocol.json",
        "inputs/source.json",
        "inputs/host.json",
        *(f"binaries/{name}" for name in binaries),
    }
    for task in tasks:
        name = task["receipt"]["name"]
        expected_run_paths.update(
            {
                f"tasks/{name}/intent.json",
                f"tasks/{name}/stdout",
                f"tasks/{name}/stderr",
                f"tasks/{name}/metrics.json",
                f"tasks/{name}/receipt.json",
            }
        )
        if task["receipt"]["result"] is not None:
            expected_run_paths.add(f"tasks/{name}/result.json")
    if run_paths != expected_run_paths:
        missing = sorted(expected_run_paths - run_paths)[:5]
        extra = sorted(run_paths - expected_run_paths)[:5]
        raise EvidenceError(
            f"retained run violates the exact file grammar (missing={missing}, extra={extra})"
        )
    resources = aggregate_task_resources(tasks)
    ratios = compute_task_ratios(tasks, rows)
    all_complete = len(rows) == values["targets"] and all(row["row_complete"] for row in rows)
    status = "complete_verified_panel" if all_complete else "incomplete_panel"
    if summary["critical_failure"] != critical_failure:
        raise EvidenceError("summary critical failure is not derivable")
    if summary["target_panel"] != targets_result or summary["rows"] != rows:
        raise EvidenceError("summary target panel or rows differ from raw task evidence")
    if summary["resources"] != resources or summary["ratios"] != ratios:
        raise EvidenceError("summary resources or ratios differ from raw task evidence")
    if (
        summary["status"] != status
        or summary["completed_rows"] != sum(row["row_complete"] for row in rows)
        or summary["evidence_class"] != values["evidence_class"]
        or summary["claim_boundary"] != frozen["claim_boundary"]
        or summary["outer_driver_accounting"] != "required_before_admission"
        or summary["full_cost_gate_passed"] is not False
        or summary["independent_external_reproduction_satisfied"] is not False
        or summary["koblitz_index_calculus_sota"] is not False
        or run_seal["panel_complete"] is not all_complete
        or run_seal["scientific_measurement_admitted"] is not False
    ):
        raise EvidenceError("run panel state or claim boundary is not derivable")

    base_inner = [
        source_binding["python"]["path"],
        path_map.original("source_repository", "scripts/run_koblitz_unknown_scalar_panel.py"),
        "run", "--profile", summary["profile"], "--output", path_map.original("run"),
        "--meter", path_map.original("source_repository", "scripts/process_meter.py"),
    ]
    allowed = [base_inner]
    if summary["profile"] == "smoke":
        allowed.append([*base_inner, "--allow-dirty"])
    if run_seal["expected_inner_command"] not in allowed:
        raise EvidenceError("sealed inner command is not reconstructible from typed paths")
    if summary["profile"] == "production" and summary["source_revision"]["dirty"]:
        raise EvidenceError("production source binding is dirty")
    return {
        "summary": summary,
        "protocol": frozen,
        "source_binding": source_binding,
        "tasks": tasks,
        "rows": rows,
        "resources": resources,
        "ratios": ratios,
        "all_complete": all_complete,
    }


def _verify_original_project_verification(
    *,
    path_map: ArchiveMap,
    manifest: dict[str, Any],
    run_seal: dict[str, Any],
    reconstructed: dict[str, Any],
) -> dict[str, Any]:
    project = require_exact(
        manifest.get("project_verification"),
        {"inventory", "seal_relative"},
        "project-verification manifest",
    )
    if project["seal_relative"] != "verification-seal.json":
        raise EvidenceError("project-verification seal path changed")
    inventory = validate_inventory_records(project["inventory"], "project-verification manifest inventory")
    if {record["path"] for record in inventory} != {"verification.json", "verification-seal.json"}:
        raise EvidenceError("project-verification tree violates the exact two-file grammar")
    if tree_inventory(path_map.archived("project_verification")) != inventory:
        raise EvidenceError("project-verification archive contains missing or extra files")
    seal, _ = read_json(
        path_map.archived("project_verification", "verification-seal.json"),
        "archived project-verification seal",
    )
    seal = require_exact(
        seal,
        {"schema", "status", "verification", "inventory", "inventory_sha256", "seal_payload_sha256"},
        "archived project-verification seal",
    )
    if seal["schema"] != PROJECT_VERIFICATION_SEAL_SCHEMA or seal["status"] != "verification_frozen":
        raise EvidenceError("archived project verification is not terminal")
    validate_self_hash(seal, "seal_payload_sha256", "archived project-verification seal")
    frozen_inventory = validate_inventory_records(seal["inventory"], "archived project-verification frozen inventory")
    if frozen_inventory != [record for record in inventory if record["path"] != "verification-seal.json"]:
        raise EvidenceError("archived project-verification frozen inventory changed")
    if canonical_sha256(frozen_inventory) != seal["inventory_sha256"]:
        raise EvidenceError("archived project-verification inventory hash changed")
    verification_identity = validate_identity_shape(seal["verification"], "project verification identity")
    if verification_identity["path"] != path_map.original("project_verification", "verification.json"):
        raise EvidenceError("project verification identity points outside its typed root")
    actual_verification_identity = identity(
        path_map.archived("project_verification", "verification.json"),
        recorded_path=verification_identity["path"],
        context="archived project verification",
    )
    if actual_verification_identity != verification_identity:
        raise EvidenceError("archived project verification bytes changed")
    verification, _ = read_json(
        path_map.archived("project_verification", "verification.json"),
        "archived project verification",
    )
    verification = require_exact(
        verification,
        {
            "schema", "status", "production", "run_root", "run_seal",
            "run_inventory_sha256", "outer_metrics", "outer_driver", "panel_status",
            "completed_rows", "expected_rows", "resources",
            "scientific_measurement_admitted", "measurement_admission_status",
            "pending_independent_payload_replay", "external_portable_verification_satisfied",
            "external_portable_verification_blocker", "independent_external_reproduction_satisfied",
            "full_cost_gate_passed", "koblitz_index_calculus_sota", "claim_boundary",
        },
        "archived project verification",
    )
    if (
        verification["schema"] != PROJECT_VERIFICATION_SCHEMA
        or verification["status"] != "run_custody_and_derivable_relationships_verified"
        or verification["run_root"] != path_map.original("run")
        or verification["run_inventory_sha256"] != run_seal["inventory_sha256"]
    ):
        raise EvidenceError("archived project verification run binding changed")
    path_map.validate_bound_file(verification["run_seal"], "run", "run-seal.json", "project-bound run seal")

    outer_manifest = require_exact(
        manifest.get("outer"),
        {"inventory", "metrics_relative", "metrics_original_identity"},
        "outer manifest",
    )
    outer_inventory = validate_inventory_records(outer_manifest["inventory"], "outer inventory")
    if {record["path"] for record in outer_inventory} != {
        "driver.metrics.json", "driver.stdout", "driver.stderr"
    }:
        raise EvidenceError("outer tree violates the exact driver receipt grammar")
    if tree_inventory(path_map.archived("outer")) != outer_inventory:
        raise EvidenceError("outer archive contains missing or extra files")
    metrics_relative = safe_relative(outer_manifest["metrics_relative"], "outer metrics path")
    if metrics_relative != "driver.metrics.json":
        raise EvidenceError("outer metrics path changed from driver.metrics.json")
    metrics_identity = validate_identity_shape(
        outer_manifest["metrics_original_identity"], "outer metrics original identity"
    )
    if metrics_identity["path"] != path_map.original("outer", metrics_relative):
        raise EvidenceError("outer metrics identity escapes the typed outer root")
    actual_metrics_identity = identity(
        path_map.archived("outer", metrics_relative),
        recorded_path=metrics_identity["path"],
        context="archived outer metrics",
    )
    if actual_metrics_identity != metrics_identity or verification["outer_metrics"] != metrics_identity:
        raise EvidenceError("outer metrics bytes or project binding changed")
    outer, _ = read_json(path_map.archived("outer", metrics_relative), "archived outer metrics")
    outer = validate_metrics(outer, "archived outer metrics")
    if verification["outer_driver"] != outer or outer["command"] != run_seal["expected_inner_command"]:
        raise EvidenceError("outer receipt does not enclose the exact sealed inner command")
    expected_watchdog = float(reconstructed["protocol"]["execution"]["whole_driver_watchdog_seconds"])
    if outer["watchdog_seconds"] != expected_watchdog or not process_clean(outer):
        raise EvidenceError("outer receipt is not a clean protocol-bounded terminal process")
    child = reconstructed["resources"]
    if outer["metrics"]["total_core_seconds"] + 1e-9 < child["total_core_seconds"]:
        raise EvidenceError("outer CPU does not enclose charged child CPU")
    if outer["metrics"]["wall_seconds"] + 1e-9 < child["summed_process_wall_seconds"]:
        raise EvidenceError("outer wall time does not enclose sequential child wall time")
    if outer["metrics"]["peak_rss_bytes"] < child["peak_process_rss_bytes"]:
        raise EvidenceError("outer RSS does not enclose child RSS")

    summary = reconstructed["summary"]
    expected_pending = [
        "domain-separated target generation and subgroup membership",
        "public factor-base discovery and materialization",
        "every relation-attempt target, decomposition witness, and relation row",
        "relation-matrix rank and modular linear algebra",
        "IC and rho recovered-scalar point verification",
        "signed-Frobenius rho operation ledger",
    ]
    if (
        verification["production"] is not (summary["profile"] == "production")
        or verification["panel_status"] != summary["status"]
        or verification["completed_rows"] != summary["completed_rows"]
        or verification["expected_rows"] != summary["expected_rows"]
        or verification["resources"] != {"charged_children": child, "whole_driver": outer["metrics"]}
        or verification["pending_independent_payload_replay"] != expected_pending
        or verification["measurement_admission_status"]
        != ("pending_independent_mathematical_payload_replay" if summary["profile"] == "production" else "operational_smoke_only")
        or verification["claim_boundary"] != summary["claim_boundary"]
        or verification["external_portable_verification_blocker"] != ORIGINAL_PORTABLE_BLOCKER
        or verification["scientific_measurement_admitted"] is not False
        or verification["external_portable_verification_satisfied"] is not False
        or verification["independent_external_reproduction_satisfied"] is not False
        or verification["full_cost_gate_passed"] is not False
        or verification["koblitz_index_calculus_sota"] is not False
    ):
        raise EvidenceError("project verification cannot be reconstructed from raw evidence")
    return {"verification": verification, "outer": outer}


def _expected_bundle_paths(
    manifest: dict[str, Any],
    run_retained: list[dict[str, Any]],
    outer_inventory: list[dict[str, Any]],
    project_inventory: list[dict[str, Any]],
    source_records: list[dict[str, Any]],
) -> set[str]:
    paths = {
        "bundle-manifest.json",
        "original/run/run-seal.json",
        "source-extra/Cargo.lock",
        "source-metadata/commit.object",
        "verifier/package_koblitz_stage23_terminal_evidence.py",
        "verifier/koblitz_stage23_terminal_evidence.py",
        "verifier/verify_koblitz_stage23_terminal_evidence.py",
    }
    paths.update(f"original/run/{record['path']}" for record in run_retained)
    paths.update(f"original/outer/{record['path']}" for record in outer_inventory)
    paths.update(f"original/project-verification/{record['path']}" for record in project_inventory)
    paths.update(f"source-tree/{record['path']}" for record in source_records)
    return paths


def verify_bundle(bundle: Path) -> dict[str, Any]:
    bundle = real_directory(bundle, "Stage-23 terminal-evidence bundle")
    initial_full_inventory = tree_inventory(bundle)
    seal, _ = read_json(bundle / "bundle-seal.json", "bundle seal")
    seal = require_exact(
        seal,
        {
            "schema", "status", "manifest", "inventory", "inventory_sha256",
            "source_run_inventory_sha256", "scientific_measurement_admitted",
            "external_portable_verification_satisfied",
            "independent_external_reproduction_satisfied", "full_cost_gate_passed",
            "koblitz_index_calculus_sota", "seal_payload_sha256",
        },
        "bundle seal",
    )
    if seal["schema"] != SEAL_SCHEMA or seal["status"] != "compact_terminal_evidence_frozen":
        raise EvidenceError("bundle seal is not terminal")
    validate_self_hash(seal, "seal_payload_sha256", "bundle seal")
    for field in (
        "scientific_measurement_admitted", "external_portable_verification_satisfied",
        "independent_external_reproduction_satisfied", "full_cost_gate_passed",
        "koblitz_index_calculus_sota",
    ):
        if seal[field] is not False:
            raise EvidenceError(f"bundle seal widened {field}")
    frozen_inventory = validate_inventory_records(seal["inventory"], "bundle frozen inventory")
    actual_inventory = [
        record for record in initial_full_inventory if record["path"] != "bundle-seal.json"
    ]
    if actual_inventory != frozen_inventory or canonical_sha256(actual_inventory) != seal["inventory_sha256"]:
        raise EvidenceError("bundle inventory changed after sealing")
    manifest_identity = validate_identity_shape(seal["manifest"], "bundle manifest identity")
    if manifest_identity["path"] != "bundle-manifest.json" or identity(
        bundle / "bundle-manifest.json", recorded_path="bundle-manifest.json"
    ) != manifest_identity:
        raise EvidenceError("bundle manifest identity changed")
    manifest, _ = read_json(bundle / "bundle-manifest.json", "bundle manifest")
    manifest = require_exact(
        manifest,
        {
            "schema", "status", "profile", "path_map", "run_partition",
            "binary_build_cross_bindings", "outer", "project_verification",
            "source_closure", "portable_tools", "claim_boundary",
        },
        "bundle manifest",
    )
    if manifest["schema"] != MANIFEST_SCHEMA or manifest["status"] != "compact_terminal_evidence_packaged":
        raise EvidenceError("bundle manifest schema or status changed")
    path_map = ArchiveMap.from_manifest(bundle, manifest["path_map"])

    partition = require_exact(
        manifest["run_partition"],
        {"source_inventory_sha256", "omission_rule", "retained", "omitted", "run_seal"},
        "run partition",
    )
    if partition["omission_rule"] != "reject a build-target file and omit exactly build-target descendants":
        raise EvidenceError("run omission rule changed")
    retained = validate_inventory_records(partition["retained"], "run retained partition")
    omitted = validate_inventory_records(partition["omitted"], "run omitted partition")
    run_seal_identity = validate_identity_shape(partition["run_seal"], "archived run-seal identity")
    if run_seal_identity["path"] != "original/run/run-seal.json" or identity(
        path_map.archived("run", "run-seal.json"),
        recorded_path="original/run/run-seal.json",
    ) != run_seal_identity:
        raise EvidenceError("archived run-seal identity changed")
    run_seal, _ = read_json(path_map.archived("run", "run-seal.json"), "archived original run seal")
    run_seal = require_exact(
        run_seal,
        {
            "schema", "status", "profile", "protocol_sha256", "source_commit",
            "summary", "inventory", "inventory_sha256", "expected_inner_command",
            "panel_complete", "scientific_measurement_admitted", "seal_payload_sha256",
        },
        "archived original run seal",
    )
    if run_seal["schema"] != RUN_SEAL_SCHEMA or run_seal["status"] != "outputs_frozen":
        raise EvidenceError("archived original run is not terminal")
    if manifest["profile"] != run_seal["profile"]:
        raise EvidenceError("manifest profile differs from the original run seal")
    validate_self_hash(run_seal, "seal_payload_sha256", "archived original run seal")
    source_inventory = validate_inventory_records(run_seal["inventory"], "archived original run inventory")
    if any(record["path"] == "build-target" for record in source_inventory):
        raise EvidenceError("original run seal contains a file at reserved build-target path")
    derived_retained = [record for record in source_inventory if not is_build_target_path(record["path"])]
    derived_omitted = [record for record in source_inventory if is_build_target_path(record["path"])]
    if retained != derived_retained or omitted != derived_omitted:
        raise EvidenceError("retained/omitted partition is not derived from the original run seal")
    if len(retained) + len(omitted) != len(source_inventory):
        raise EvidenceError("retained/omitted partition is not exhaustive")
    if (
        canonical_sha256(source_inventory) != run_seal["inventory_sha256"]
        or partition["source_inventory_sha256"] != run_seal["inventory_sha256"]
        or seal["source_run_inventory_sha256"] != run_seal["inventory_sha256"]
    ):
        raise EvidenceError("source run inventory hash binding changed")
    actual_run = tree_inventory(path_map.archived("run"), excluded={"run-seal.json"})
    if actual_run != retained:
        raise EvidenceError("compact run archive contains missing, omitted, or extra files")
    if any(is_build_target_path(record["path"]) for record in actual_run):
        raise EvidenceError("compact run archive retained build-target bytes")

    summary, _ = read_json(path_map.archived("run", "run-summary.json"), "archived run summary")
    bindings = manifest["binary_build_cross_bindings"]
    if not isinstance(bindings, list):
        raise EvidenceError("binary/build cross-bindings must be a list")
    expected_bindings: list[dict[str, Any]] = []
    inventory_by_path = {record["path"]: record for record in source_inventory}
    summary_binaries = summary.get("binaries", {})
    if not isinstance(summary_binaries, dict):
        raise EvidenceError("summary binary inventory is malformed")
    for name in sorted(summary_binaries):
        if name not in BINARY_NAMES:
            raise EvidenceError(f"unexpected summary binary {name}")
        retained_path = f"binaries/{name}"
        omitted_path = f"build-target/release/examples/{name}"
        retained_record = inventory_by_path.get(retained_path)
        omitted_record = inventory_by_path.get(omitted_path)
        if retained_record is None or omitted_record is None:
            raise EvidenceError(f"binary {name} lacks a sealed fresh-build counterpart")
        if (retained_record["bytes"], retained_record["sha256"]) != (
            omitted_record["bytes"], omitted_record["sha256"]
        ):
            raise EvidenceError(f"binary {name} differs from its omitted fresh-build output")
        claimed = validate_identity_shape(summary_binaries[name], f"summary binary {name}")
        if claimed["path"] != path_map.original("run", retained_path) or (
            claimed["bytes"], claimed["sha256"]
        ) != (retained_record["bytes"], retained_record["sha256"]):
            raise EvidenceError(f"summary binary {name} binding changed")
        expected_bindings.append(
            {
                "name": name,
                "retained_run_path": retained_path,
                "retained_identity": retained_record,
                "omitted_build_path": omitted_path,
                "omitted_identity": omitted_record,
                "bytes_equal_by_original_run_seal": True,
            }
        )
    if bindings != expected_bindings:
        raise EvidenceError("binary/build cross-bindings are incomplete or forged")

    portable_tools = require_exact(
        manifest["portable_tools"], {"packager", "core", "verifier"}, "portable tools"
    )
    for key, relative in (
        ("packager", "verifier/package_koblitz_stage23_terminal_evidence.py"),
        ("core", "verifier/koblitz_stage23_terminal_evidence.py"),
        ("verifier", "verifier/verify_koblitz_stage23_terminal_evidence.py"),
    ):
        tool_identity = validate_identity_shape(portable_tools[key], f"portable {key}")
        if tool_identity["path"] != relative or identity(bundle / relative, recorded_path=relative) != tool_identity:
            raise EvidenceError(f"portable {key} source identity changed")

    reconstructed = _verify_run_reconstruction(
        path_map=path_map,
        manifest=manifest,
        run_seal=run_seal,
        retained=retained,
    )
    if manifest["profile"] != reconstructed["summary"]["profile"]:
        raise EvidenceError("manifest profile differs from reconstructed run profile")
    project = _verify_original_project_verification(
        path_map=path_map,
        manifest=manifest,
        run_seal=run_seal,
        reconstructed=reconstructed,
    )

    claim = require_exact(
        manifest["claim_boundary"],
        {
            "project_custody_reconstruction_only",
            "independent_mathematical_payload_replay_completed",
            "pending_independent_mathematical_payload_replay",
            "scientific_measurement_admitted", "external_portable_verification_satisfied",
            "independent_external_reproduction_satisfied", "full_cost_gate_passed",
            "koblitz_index_calculus_sota",
        },
        "bundle claim boundary",
    )
    if (
        claim["project_custody_reconstruction_only"] is not True
        or claim["independent_mathematical_payload_replay_completed"] is not False
        or claim["pending_independent_mathematical_payload_replay"] != PENDING_MATH_REPLAY
    ):
        raise EvidenceError("bundle claim boundary hides or changes pending mathematical replay")
    for field in (
        "scientific_measurement_admitted", "external_portable_verification_satisfied",
        "independent_external_reproduction_satisfied", "full_cost_gate_passed",
        "koblitz_index_calculus_sota",
    ):
        if claim[field] is not False:
            raise EvidenceError(f"bundle manifest widened {field}")

    outer_inventory = validate_inventory_records(manifest["outer"]["inventory"], "outer inventory")
    project_inventory = validate_inventory_records(
        manifest["project_verification"]["inventory"], "project-verification inventory"
    )
    source_records = manifest["source_closure"]["tracked_files"]
    expected_paths = _expected_bundle_paths(
        manifest, retained, outer_inventory, project_inventory, source_records
    )
    actual_paths = {record["path"] for record in actual_inventory}
    if actual_paths != expected_paths:
        missing = sorted(expected_paths - actual_paths)[:5]
        extra = sorted(actual_paths - expected_paths)[:5]
        raise EvidenceError(f"bundle has missing or extra files (missing={missing}, extra={extra})")
    result = {
        "schema": VERIFY_SCHEMA,
        "status": "compact_terminal_evidence_structurally_verified",
        "profile": reconstructed["summary"]["profile"],
        "panel_status": reconstructed["summary"]["status"],
        "task_count": len(reconstructed["tasks"]),
        "row_count": len(reconstructed["rows"]),
        "completed_rows": reconstructed["summary"]["completed_rows"],
        "resources": reconstructed["resources"],
        "ratios": reconstructed["ratios"],
        "source_commit": run_seal["source_commit"],
        "source_tree_oid": manifest["source_closure"]["root_tree_oid"],
        "run_inventory_sha256": run_seal["inventory_sha256"],
        "omitted_build_files": len(omitted),
        "binary_build_cross_bindings": len(bindings),
        "original_project_verification_preserved": bool(project),
        "bundled_verifier_sources_executed": False,
        "typed_archive_path_map_verified": True,
        "exact_task_graph_commands_resources_rows_and_ratios_reconstructed": True,
        "independent_mathematical_payload_replay_completed": False,
        "pending_independent_mathematical_payload_replay": PENDING_MATH_REPLAY,
        "scientific_measurement_admitted": False,
        "external_portable_verification_satisfied": False,
        "independent_external_reproduction_satisfied": False,
        "full_cost_gate_passed": False,
        "koblitz_index_calculus_sota": False,
    }
    final_full_inventory = tree_inventory(bundle)
    if final_full_inventory != initial_full_inventory:
        raise EvidenceError("bundle changed between authenticated inventory and verifier return")
    final_authenticated_inventory = [
        record for record in final_full_inventory if record["path"] != "bundle-seal.json"
    ]
    if final_authenticated_inventory != actual_inventory:
        raise EvidenceError("authenticated bundle inventory changed before verifier return")
    return result
