"""Pack a complete or interrupted campaign with an explicit read-status seal.

The upload action must not traverse a raw million-file tree. This command
runs in a separate `if: always()` step after measurement and never edits
measured receipts. It preserves the compressed stream and tar diagnostics
even when a concurrent read prevents a trustworthy archive claim.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import time


def pack(source, archive, *, cargo_lock=None, tar_factory=None):
    source, archive = Path(source).resolve(), Path(archive).resolve()
    if archive.is_relative_to(source):
        raise ValueError('archive cannot live inside the measured result tree')
    existed = source.is_dir()
    source.mkdir(parents=True, exist_ok=True)
    if cargo_lock is not None:
        lock = Path(cargo_lock).resolve()
        if lock.is_file():
            destination = source / 'workflow-Cargo.lock'
            if not destination.exists():
                shutil.copy2(lock, destination)
    capture = dict(schema_version=1, source_existed_before_pack=existed,
                   captured_unix=time.time(),
                   state_present=(source/'tournament/state.json').is_file(),
                   summary_present=(source/'summary.json').is_file(),
                   gate_present=(source/'family-gate.json').is_file())
    (source/'capture.json').write_text(json.dumps(capture, sort_keys=True)+'\n')
    archive.parent.mkdir(parents=True, exist_ok=True)
    if archive.exists():
        raise FileExistsError(archive)
    # A timed-out producer may still be writing when tar walks the million-file
    # tree. GNU tar exits 1 for changed files. Retain that compressed stream as
    # explicitly unverified evidence so the always-run upload can preserve it;
    # never silently treat it as a complete snapshot.
    tar_stderr = archive.with_name(archive.name+'.tar.stderr')
    tar_factory = subprocess.Popen if tar_factory is None else tar_factory
    with tar_stderr.open('xb') as errors:
        tar = tar_factory(
            ['tar', '-C', str(source.parent), '-cf', '-', source.name],
            stdout=subprocess.PIPE, stderr=errors)
        try:
            compressed = subprocess.run(['zstd', '-T2', '-3', '-o', str(archive)],
                                        stdin=tar.stdout, check=False)
        finally:
            tar.stdout.close()
        tar_exit_code = tar.wait()
    tar_stderr_bytes = tar_stderr.read_bytes()
    archived = archive.is_file()
    if archived:
        with archive.open('rb') as stream:
            archive_sha256 = hashlib.file_digest(stream, 'sha256').hexdigest()
    else:
        archive_sha256 = None
    status = ('SOURCE_MISSING' if not existed
              else 'COMPRESSION_FAILURE' if compressed.returncode != 0 or not archived
              else 'ARCHIVE_READ_UNVERIFIED' if tar_exit_code != 0 or tar_stderr_bytes
              else 'ARCHIVE_READ_COMPLETE')
    manifest = dict(schema_version=2, archive=archive.name,
                    sha256=archive_sha256,
                    byte_count=archive.stat().st_size if archived else None,
                    capture=capture, pack_status=status,
                    tar_exit_code=tar_exit_code,
                    compression_exit_code=compressed.returncode,
                    tar_stderr_file=tar_stderr.name,
                    tar_stderr_sha256=hashlib.sha256(tar_stderr_bytes).hexdigest(),
                    tar_stderr_bytes=len(tar_stderr_bytes),
                    independent_audit_required=True)
    archive.with_name(archive.name+'.manifest.json').write_text(
        json.dumps(manifest, sort_keys=True, indent=2)+'\n')
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--archive', type=Path, required=True)
    parser.add_argument('--cargo-lock', type=Path)
    args = parser.parse_args()
    manifest = pack(args.source, args.archive, cargo_lock=args.cargo_lock)
    print(json.dumps(manifest), flush=True)
    if manifest['pack_status'] != 'ARCHIVE_READ_COMPLETE':
        raise SystemExit(1)


if __name__ == '__main__':
    main()
