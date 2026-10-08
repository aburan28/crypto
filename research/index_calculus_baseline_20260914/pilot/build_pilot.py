"""Fetch pinned third-party files into a cache and build WDSat outside the repo."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parent


def verified_sources(cache, fetch=False):
    for item in json.loads((ROOT / 'source_manifest.json').read_text()):
        path = cache / item['path']
        if not path.exists() and fetch:
            with urlopen(item['url'], timeout=60) as response:
                data = response.read()
            if hashlib.sha256(data).hexdigest() != item['sha256']:
                raise ValueError(f'download hash mismatch: {path}')
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(data)
        if hashlib.sha256(path.read_bytes()).hexdigest() != item['sha256']:
            raise ValueError(f'source hash mismatch: {path}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache-dir', required=True, type=Path)
    parser.add_argument('--build-dir', required=True, type=Path,
                        help='Must not already exist; preserves previous builds.')
    parser.add_argument('--offline', action='store_true')
    args = parser.parse_args()
    if args.build_dir.exists():
        parser.error('--build-dir must not already exist')
    cache = args.cache_dir.resolve()
    verified_sources(cache, fetch=not args.offline)
    dest = args.build_dir.resolve()
    dest.mkdir(parents=True, exist_ok=False)
    cmd = ['gcc', '-O3', '-Wall'] + [str(p) for p in sorted(
        (cache / 'vendor/WDSat/src').glob('*.c'))] + ['-lm', '-o', str(dest / 'wdsat_solver')]
    result = subprocess.run(cmd, capture_output=True, text=True)
    (dest / 'build.stdout').write_text(result.stdout)
    (dest / 'build.stderr').write_text(result.stderr)
    metadata = {'exit_code': result.returncode, 'compiler_command': cmd,
                'compiler': subprocess.check_output(['gcc', '--version'], text=True).splitlines()[0],
                'source_manifest_sha256': hashlib.sha256((ROOT / 'source_manifest.json').read_bytes()).hexdigest(),
                'binary_sha256': hashlib.sha256((dest / 'wdsat_solver').read_bytes()).hexdigest()
                    if result.returncode == 0 else None}
    (dest / 'build.json').write_text(json.dumps(metadata, indent=2) + '\n')
    print(json.dumps(metadata))
    raise SystemExit(result.returncode)


if __name__ == '__main__':
    main()
