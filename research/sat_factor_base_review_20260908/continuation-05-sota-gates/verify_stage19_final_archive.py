#!/usr/bin/env python3
"""Replay the pinned completed Stage 19 archive offline on a different host.

The live runner is unchanged. Only historical command/authorization paths and
the recorded CA identity are adapted in an isolated copy of its exact verifier.
No caller-supplied mapping, current-host CA, or network operation is accepted.
"""

from __future__ import annotations

import argparse
from contextlib import contextmanager
import hashlib
import importlib.util
import json
from pathlib import Path, PurePosixPath
import socket
import stat
import subprocess
import tempfile
from unittest.mock import patch


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
ARCHIVE_COMMIT = "b8cc5d684dbaf4144d2d9fdc6d53da567142ce45"
ARCHIVE_TREE = "e6a6d04f51cd47b9ced85238a6f59c4f0a22f795"
FIRST_COMMIT = "5a09482aea0ce57239b9a06a7cf93c03121730a7"
EXECUTION_COMMIT = "a1001c18524e1e84d83f3e55498b9a26242c1a13"
GATE_PATH = PurePosixPath("research/sat_factor_base_review_20260908/continuation-05-sota-gates")
ARTIFACT = GATE_PATH / "stage-19-magma-calculator-panel-amendment-03-20260910"
EXPECTED = GATE_PATH / "stage-19-amendment-03-final-verification.json"
CHILD = GATE_PATH / "post_stage19_amendment03_magma_calculator_request.py"
VERIFIER = GATE_PATH / "verify_stage19_amendment03_magma_calculator_panel.py"
HISTORICAL_ROOT = PurePosixPath("/Volumes/SSD990/crypto-koblitz-magma-calculator-panel-20260909")
HISTORICAL_PYTHON = "/Library/Frameworks/Python.framework/Versions/3.13/bin/python3.13"
RECORDED_CA = {
    "path": "/etc/ssl/cert.pem", "bytes": 333483,
    "sha256": "9dae8d76e55cb08991f2b672d58999ea15560d910759c16b544f843bdffbb994",
}


class ArchiveError(RuntimeError):
    pass


def git(repo: Path, *args: str) -> bytes:
    result = subprocess.run(["git", *args], cwd=repo, capture_output=True, check=False)
    if result.returncode:
        raise ArchiveError(result.stderr.decode(errors="replace").strip())
    return result.stdout


def regular_bytes(root: Path, relative: str) -> bytes:
    rel = PurePosixPath(relative)
    if rel.is_absolute() or ".." in rel.parts:
        raise ArchiveError("archive path escapes its checkout")
    path = root
    for part in rel.parts:
        path /= part
        if path.is_symlink():
            raise ArchiveError(f"archive path is a symlink: {relative}")
    info = path.stat()
    if not stat.S_ISREG(info.st_mode) or info.st_nlink != 1:
        raise ArchiveError(f"archive path is not a single-link file: {relative}")
    return path.read_bytes()


def verify_current_archive(repo: Path) -> dict:
    """Authenticate the full current result against an immutable Git commit."""
    if git(repo, "rev-parse", ARCHIVE_COMMIT + "^{tree}").decode().strip() != ARCHIVE_TREE:
        raise ArchiveError("archive tree identity changed")
    entries = git(repo, "ls-tree", "-r", ARCHIVE_COMMIT, "--", str(ARTIFACT), str(EXPECTED))
    rows = []
    expected_names = set()
    for line in entries.decode().splitlines():
        metadata, relative = line.split("\t", 1)
        mode, kind, oid = metadata.split()
        if mode != "100644" or kind != "blob":
            raise ArchiveError("unexpected archive Git entry")
        data = regular_bytes(repo, relative)
        if data != git(repo, "cat-file", "blob", oid):
            raise ArchiveError(f"current archive differs from frozen Git blob: {relative}")
        rows.append({"path": relative, "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()})
        expected_names.add(relative)
    actual_names = {
        path.relative_to(repo).as_posix()
        for path in (repo / ARTIFACT).rglob("*") if not path.is_dir() or path.is_symlink()
    } | {str(EXPECTED)}
    if len(rows) != 96 or actual_names != expected_names:
        raise ArchiveError("completed archive file inventory changed")
    encoded = (json.dumps(rows, indent=2, sort_keys=True) + "\n").encode()
    return {"files": len(rows), "sha256": hashlib.sha256(encoded).hexdigest()}


@contextmanager
def detached_archive(repo: Path):
    with tempfile.TemporaryDirectory(prefix="stage19-offline-archive-") as temporary:
        destination = Path(temporary) / "checkout"
        git(repo, "worktree", "add", "--detach", str(destination), ARCHIVE_COMMIT)
        try:
            yield destination
        finally:
            git(repo, "worktree", "remove", "--force", str(destination))


def historical_command(task: dict) -> list[str]:
    ordinal = task["ordinal"]
    if type(ordinal) is not int or not 1 <= ordinal <= 9:
        raise ArchiveError("unexpected historical task ordinal")
    expected_stem = f"{ordinal:02d}-{task['id'].replace('/', '-')}"
    if task["named_input"]["path"] != f"inputs/{expected_stem}.magma":
        raise ArchiveError("historical input name does not match the task")
    return [
        HISTORICAL_PYTHON, str(HISTORICAL_ROOT / CHILD), "--attempt",
        str(HISTORICAL_ROOT / ARTIFACT / "attempts" / expected_stem), "--execute-child",
    ]


def deny_network(*_args, **_kwargs):
    raise ArchiveError("network access is forbidden during historical verification")


def verify_historical_commands(snapshot: Path) -> None:
    artifact = snapshot / ARTIFACT
    plan = json.loads((artifact / "plan.json").read_text())
    for task in plan["tasks"]:
        command = historical_command(task)
        attempt = artifact / "attempts" / Path(task["named_input"]["path"]).stem
        metrics = json.loads((attempt / "transport-metrics.json").read_text())
        authorization = json.loads((attempt / "launch-authorization.consumed.json").read_text())
        if metrics["command"] != command or authorization["child"]["argv"] != command[1:]:
            raise ArchiveError("recorded command does not match the fixed historical host context")


def replay(snapshot: Path) -> dict:
    # The snapshot is created from the immutable archive commit before importing
    # verifier code; never import mutable current-checkout implementation here.
    if git(snapshot, "rev-parse", "HEAD").decode().strip() != ARCHIVE_COMMIT:
        raise ArchiveError("replay is not at the pinned archive commit")
    if git(snapshot, "status", "--porcelain", "--untracked-files=all"):
        raise ArchiveError("historical replay checkout is dirty")
    verify_historical_commands(snapshot)
    with patch.object(socket.socket, "connect", deny_network), patch.object(socket.socket, "connect_ex", deny_network), patch.object(socket, "create_connection", deny_network):
        spec = importlib.util.spec_from_file_location("stage19_pinned_archive_verifier", snapshot / VERIFIER)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        artifact = snapshot / ARTIFACT
        run = json.loads((artifact / "run.json").read_text())
        binding = run["execution_binding"]
        if binding["first_preexecution_commit"] != FIRST_COMMIT or binding["execution_commit"] != EXECUTION_COMMIT:
            raise ArchiveError("historical pre-execution commits changed")
        if module.CHILD.expected_ca_bundle_record() != RECORDED_CA or binding["ca_bundle"] != RECORDED_CA:
            raise ArchiveError("historical CA record changed")
        original_authorization = module.CHILD.expected_launch_authorization

        def transport_command(observed_artifact, task):
            if observed_artifact != artifact:
                raise ArchiveError("unexpected replay artifact")
            return historical_command(task)

        def authorization(attempt, plan, task, start, artifact=artifact):
            expected = original_authorization(attempt, plan, task, start, artifact=artifact)
            command = transport_command(artifact, task)
            if attempt != artifact / "attempts" / Path(task["named_input"]["path"]).stem:
                raise ArchiveError("unexpected replay attempt path")
            expected["child"]["argv"] = command[1:]
            return expected

        # These are interpretation adapters for recorded host metadata. Neither
        # touches the live transport, mutates receipts, nor skips the verifier's
        # full binding, source, authorization, chronology and response checks.
        module.expected_transport_command = transport_command
        module.CHILD.expected_launch_authorization = authorization
        module.CHILD.ca_bundle_record = lambda: dict(RECORDED_CA)
        summary = module.verify(artifact)
        rendered = module.canonical_bytes(summary)
        if rendered != (snapshot / EXPECTED).read_bytes():
            raise ArchiveError("historical final summary did not regenerate byte-for-byte")
        return summary


def verify(repo: Path = REPO) -> dict:
    verify_current_archive(repo)
    with detached_archive(repo) as snapshot:
        return replay(snapshot)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--expected", type=Path, default=REPO / EXPECTED)
    args = parser.parse_args()
    try:
        summary = verify()
        rendered = (json.dumps(summary, indent=2, sort_keys=True) + "\n").encode()
        if args.expected.read_bytes() != rendered:
            raise ArchiveError("expected final summary differs from authenticated replay")
        print(rendered.decode(), end="")
    except (ArchiveError, OSError, ValueError) as error:
        parser.exit(2, f"error: {error}\n")


if __name__ == "__main__":
    main()
