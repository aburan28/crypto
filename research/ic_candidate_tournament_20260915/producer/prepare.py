#!/usr/bin/env python3
"""Build reviewable derivatives of hash-checked archived optimized sources.

This only materializes code. It does not label a source qualified, measure a
candidate, or overwrite the original archive/previous output directory.
"""
import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import shutil
import subprocess
import sys
import tarfile

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from identity import sha256, write_immutable
from oracle import require

SOURCES = {
    'scaled': ('round-0023', 'source_candidates/scaled',
               '55154f73c35b1f55240b35fd1a5e8e2488c6df41444b5114949e41741c0c3db1'),
    'both': ('round-0020', 'source_candidates/both',
             '563eb460f29d9ef09a2adbde4770a466d3dc16566f5f08f1e6d8238325b104d4'),
}


def filehash(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def extract_source(name, destination):
    """Verify the whole archive, extract only the selected sealed source tree."""
    round_name, relative, _ = SOURCES[name]
    entries = json.loads((HERE.parent/'evidence/manifest.json').read_text())['archives']
    entry = next(e for e in entries if e['file'] == round_name+'.tar.zst')
    archive = HERE.parent/'evidence'/entry['file']
    require(filehash(archive) == entry['sha256'], 'archive digest mismatch')
    prefix = f'runs/{round_name}/{relative}/'
    target_prefix = prefix+'source/'
    destination.mkdir(parents=True, exist_ok=False)
    process = subprocess.Popen(['zstd', '-dq', '-c', str(archive)], stdout=subprocess.PIPE)
    try:
        with tarfile.open(fileobj=process.stdout, mode='r|') as tar:
            for item in tar:
                if not (item.name.startswith(target_prefix) or item.name == prefix+'source-manifest.json'):
                    continue
                require(item.isfile(), 'non-file source archive entry')
                relative_path = PurePosixPath(item.name.removeprefix(prefix))
                require(not relative_path.is_absolute() and '..' not in relative_path.parts,
                        'unsafe source archive path')
                target = destination/relative_path
                target.parent.mkdir(parents=True, exist_ok=True)
                with target.open('xb') as stream:
                    stream.write(tar.extractfile(item).read())
        require(process.wait() == 0, 'archive decompressor failed')
    finally:
        process.stdout.close()
        if process.poll() is None:
            process.terminate()
        process.wait()


def verify_source(root, expected):
    manifest = json.loads((root/'source-manifest.json').read_text())
    require(sha256(manifest) == expected, 'unexpected archived source manifest')
    source = root/'source'
    files = {str(p.relative_to(source)) for p in source.rglob('*') if p.is_file()}
    require(files == set(manifest), 'source inventory differs from manifest')
    for relative, expected_digest in manifest.items():
        path = source/relative
        require(not path.is_symlink() and path.resolve().is_relative_to(source.resolve()),
                'unsafe source path')
        require(filehash(path) == expected_digest, 'source digest differs: '+relative)
    return manifest


def prepare(name, output, restored=None, *, instrument=True):
    require(name in ('scaled', 'pairinv', 'both'), 'unknown archived source')
    base_name = 'scaled' if name == 'pairinv' else name
    round_name, relative, expected = SOURCES[base_name]
    output = output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    if restored is None:
        origin = output/'archive-input'
        extract_source(base_name, origin)
    else:
        origin = restored.resolve()/'runs'/round_name/relative
    manifest = verify_source(origin, expected)
    source = output/'source'
    shutil.copytree(origin/'source', source)
    patches = []
    if name == 'pairinv':
        patches.append(HERE.parent/'campaign_20260916/round24-pairinv.patch')
    if instrument:
        if name == 'both':
            patches.append(HERE/'both-test-convention.patch')
        patches.append(HERE/'exclusive-phases.patch')
    for patch in patches:
        # Check every hunk before writing. Some system patch implementations
        # return success after printing that an out-of-range hunk was ignored.
        subprocess.run(['git', 'apply', '--check', str(patch)], cwd=source, check=True)
        subprocess.run(['git', 'apply', str(patch)], cwd=source, check=True)
    if instrument:
        shutil.copy2(HERE/'ic_phase.rs', source/'src/cryptanalysis/ic_phase.rs')
        module = source/'src/cryptanalysis/mod.rs'
        with module.open('a') as stream:
            stream.write('\npub mod ic_phase;\n')
    derived = {str(p.relative_to(source)): filehash(p) for p in sorted(source.rglob('*')) if p.is_file()}
    write_immutable(output/'source-manifest.json', derived)
    receipt = {'status': 'SOURCE_MATERIALIZED_NOT_QUALIFIED', 'reference': name,
               'instrumented': instrument, 'archive_source_sha256': sha256(manifest),
               'source_manifest_sha256': sha256(derived),
               'patches': [{'name': p.name, 'sha256': filehash(p)} for p in patches],
               'phase_module_sha256': filehash(HERE/'ic_phase.rs') if instrument else None,
               'candidate_id': None, 'total_operations': None}
    write_immutable(output/'preparation.json', receipt)
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference', choices=('scaled', 'pairinv', 'both'), required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--restored-root', type=Path)
    parser.add_argument('--original', action='store_true')
    args = parser.parse_args()
    print(json.dumps(prepare(args.reference, args.out, args.restored_root, instrument=not args.original)))


if __name__ == '__main__':
    main()
