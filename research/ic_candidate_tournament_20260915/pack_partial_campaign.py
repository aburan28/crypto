"""Pack a complete or interrupted campaign as one upload-safe artifact.

The upload action traversed the first campaign's raw million-file tree and
failed after the job timed out. This command runs in a separate `if: always()`
step with time reserved after measurement. It never edits measured receipts.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import time


def pack(source, archive, *, cargo_lock=None):
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
    tar = subprocess.Popen(['tar', '-C', str(source.parent), '-cf', '-', source.name],
                           stdout=subprocess.PIPE)
    try:
        compressed = subprocess.run(['zstd', '-T2', '-3', '-o', str(archive)],
                                    stdin=tar.stdout, check=False)
    finally:
        tar.stdout.close()
    if tar.wait() != 0 or compressed.returncode != 0:
        archive.unlink(missing_ok=True)
        raise RuntimeError('partial campaign archive failed')
    with archive.open('rb') as stream:
        digest = hashlib.file_digest(stream, 'sha256').hexdigest()
    manifest = dict(schema_version=1, archive=archive.name, sha256=digest,
                    byte_count=archive.stat().st_size, capture=capture)
    archive.with_name(archive.name+'.manifest.json').write_text(
        json.dumps(manifest, sort_keys=True, indent=2)+'\n')
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--archive', type=Path, required=True)
    parser.add_argument('--cargo-lock', type=Path)
    args = parser.parse_args()
    print(json.dumps(pack(args.source, args.archive, cargo_lock=args.cargo_lock)),
          flush=True)


if __name__ == '__main__':
    main()
