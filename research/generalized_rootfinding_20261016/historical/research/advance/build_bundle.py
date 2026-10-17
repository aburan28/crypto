"""Package the frozen release plus experimental source and complete evidence."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import zipfile

PROJECT = Path(__file__).resolve().parents[2]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    if output.exists():
        raise SystemExit('Refusing to overwrite an existing source archive')
    previous = json.loads((PROJECT / 'BUNDLE_MANIFEST.json').read_text())
    files = set(previous['files']) | {'BUNDLE_MANIFEST.json'}
    for name, metadata in previous['files'].items():
        if digest(PROJECT / name) != metadata['sha256']:
            raise SystemExit(f'Frozen version0.2 file changed: {name}')
    for directory in ('research/continuation', 'research/advance'):
        for path in (PROJECT / directory).rglob('*'):
            if not path.is_file():
                continue
            parts = path.relative_to(PROJECT).parts
            if any(p == '__pycache__' or p == 'target' or p.startswith('deps') for p in parts):
                continue
            if path.suffix in ('.pyc', '.so') or path.name.endswith('.tmp'):
                continue
            files.add(path.relative_to(PROJECT).as_posix())
    manifest_name = 'research/advance/BUNDLE_MANIFEST.json'
    files.discard(manifest_name)
    manifest = {
        'schema': 'generalized_rootfinding_advanced_source_bundle_v1',
        'created_at_utc': datetime.now(timezone.utc).isoformat(),
        'status': 'EXPERIMENTAL',
        'performance_release_gate': 'INCOMPLETE',
        'read_first': 'research/advance/ADVANCE.md',
        'measured_results': 'research/advance/performance/RESULTS.md',
        'correctness_cases': 89,
        'synthetic_harness_tests': 37,
        'native_adapter_tests': 2,
        'exact_symbolic_verified_inputs': 30,
        'previous_release_files_unchanged': len(previous['files']),
        'excluded': ['vendor checkout', 'compiled executables and target directories',
                     'optional native SAT wheels/dependencies', 'Python caches'],
        'external_source': {'repository': 'https://github.com/aburan28/crypto',
                            'revision': '3d737486869d1ca77e14ad56bdaf376916eecce7'},
        'hash_scope': 'All listed source/evidence files; this manifest excludes itself',
        'files': {name: {'bytes': (PROJECT / name).stat().st_size,
                         'sha256': digest(PROJECT / name)} for name in sorted(files)},
    }
    (PROJECT / manifest_name).write_text(json.dumps(manifest, indent=2) + '\n')
    files.add(manifest_name)
    with zipfile.ZipFile(output, 'w', compression=zipfile.ZIP_DEFLATED,
                         compresslevel=6) as archive:
        for name in sorted(files):
            archive.write(PROJECT / name, 'generalized_roots/' + name)
    with zipfile.ZipFile(output) as archive:
        if archive.testzip() is not None:
            raise SystemExit('Archive CRC verification failed')
        for name, metadata in manifest['files'].items():
            data = archive.read('generalized_roots/' + name)
            if hashlib.sha256(data).hexdigest() != metadata['sha256']:
                raise SystemExit(f'Archive source hash differs: {name}')
    receipt = {'schema': 'advanced_source_archive_verification_v1',
               'archive': str(output), 'bytes': output.stat().st_size,
               'sha256': digest(output), 'files': len(files),
               'indexed_files': len(manifest['files']), 'crc_verified': True,
               'all_indexed_hashes_verified': True, 'release_gate': 'INCOMPLETE'}
    output.with_suffix('.verification.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt))


if __name__ == '__main__':
    main()
