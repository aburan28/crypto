#!/usr/bin/env python3
"""Build two Linux CaDiCaL binaries from pinned source; no measured producer."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

from source_manifest import scan

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def logged(command: list[str], *, cwd: Path, prefix: Path,
           env: dict[str, str] | None = None, timeout: int = 600) -> None:
    with prefix.with_suffix('.stdout.txt').open('wb') as stdout:
        with prefix.with_suffix('.stderr.txt').open('wb') as stderr:
            subprocess.run(command, cwd=cwd, env=env, stdout=stdout,
                           stderr=stderr, timeout=timeout, check=True)


def tool(command: list[str]) -> str:
    return subprocess.check_output(command, text=True, stderr=subprocess.STDOUT,
                                   timeout=20).strip()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--probe-receipt', type=Path, required=True)
    args = parser.parse_args()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    prepare_path = HERE / 'PREPARE.json'
    prepare = json.loads(prepare_path.read_text())
    expected = prepare['cadical']
    receipt = {'schema': 'k0-dag-dimacs-linux-cadical-build-v1',
               'utc': datetime.now(timezone.utc).isoformat(),
               'decision': 'PREPARE_FAILURE',
               'prepare_sha256': sha(prepare_path),
               'attempt1_receipt_sha256': prepare['attempt1_receipt_sha256'],
               'runner_image_os': os.getenv('ImageOS'),
               'runner_image_version': os.getenv('ImageVersion'),
               'platform': platform.platform(),
               'machine': platform.machine(),
               'python_executable': sys.executable,
               'python_version': sys.version,
               'checkout_head': subprocess.check_output(
                   ['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
               'builds': []}
    try:
        if sys.platform != 'linux' or platform.machine() != 'x86_64':
            raise RuntimeError('NOT_ADMITTED: target must be Linux x86_64')
        probe = json.loads(args.probe_receipt.read_text())
        if (probe['decision'] != 'PASS_HARMLESS_CAPS_ONLY' or
                [row['cap_bytes'] for row in probe['caps']] !=
                prepare['target']['hard_rlimit_as_caps_bytes'] or
                not receipt['runner_image_version'] or
                probe['runner_image_version'] != receipt['runner_image_version'] or
                probe['checkout_head'] != receipt['checkout_head'] or
                probe['python_executable'] != receipt['python_executable']):
            raise RuntimeError('NOT_ADMITTED: matching hard-RLIMIT_AS probe absent')
        receipt['probe_receipt_sha256'] = sha(args.probe_receipt)
        receipt['toolchain'] = {
            'uname': tool(['uname', '-a']),
            'gcc': tool(['/usr/bin/gcc', '--version']),
            'gxx': tool(['/usr/bin/g++', '--version']),
            'ar': tool(['/usr/bin/ar', '--version']),
            'make': tool(['/usr/bin/make', '--version']),
            'ldd': tool(['ldd', '--version']),
            'os_release': Path('/etc/os-release').read_text(),
        }
        build_env = os.environ.copy()
        build_env.update({'LC_ALL': 'C', 'TZ': 'UTC',
                          'SOURCE_DATE_EPOCH': str(expected['source_date_epoch']),
                          'CC': '/usr/bin/gcc', 'CXX': '/usr/bin/g++',
                          'AR': '/usr/bin/ar'})
        for label in ('a', 'b'):
            source = out / f'source-{label}'
            logged(['git', 'clone', '--depth', '1', '--branch', expected['tag'],
                    '--single-branch', expected['source_repository'], str(source)],
                   cwd=out, prefix=out / f'{label}-clone')
            got = scan(source)
            for key in ('commit', 'tree', 'tracked_files', 'source_manifest_sha256'):
                if got[key] != expected[key]:
                    raise RuntimeError(f'NOT_ADMITTED: source {key} mismatch')
            if got['version'] != '3.0.1':
                raise RuntimeError('NOT_ADMITTED: source VERSION mismatch')
            logged(expected['configure_command'], cwd=source,
                   prefix=out / f'{label}-configure', env=build_env)
            logged(expected['build_command'], cwd=source,
                   prefix=out / f'{label}-make', env=build_env)
            binary = source / 'build/cadical'
            version = tool([str(binary), '--version'])
            if '3.0.1' not in version:
                raise RuntimeError('NOT_ADMITTED: built CaDiCaL version mismatch')
            receipt['builds'].append({'label': label, 'source': got,
                                      'binary_sha256': sha(binary),
                                      'binary_bytes': binary.stat().st_size,
                                      'version_output': version})
        if receipt['builds'][0]['binary_sha256'] != receipt['builds'][1]['binary_sha256']:
            raise RuntimeError('NOT_ADMITTED: independent Linux builds differ')
        shutil.copy2(out / 'source-a/build/cadical', out / 'cadical')
        receipt['linux_binary_sha256'] = sha(out / 'cadical')
        receipt['linux_binary_bytes'] = (out / 'cadical').stat().st_size
        receipt['binary_ldd'] = tool(['ldd', str(out / 'cadical')])
        receipt['decision'] = 'PASS_BUILD_ONLY_NO_MEASURED_CHILD'
    except Exception as exc:
        receipt['error'] = f'{type(exc).__name__}: {exc}'
    (out / 'PREPARE_RECEIPT.json').write_text(json.dumps(receipt, sort_keys=True,
                                                       indent=2) + '\n')
    print(json.dumps({'decision': receipt['decision'],
                      'builds': len(receipt['builds'])}, sort_keys=True))
    return 0 if receipt['decision'] == 'PASS_BUILD_ONLY_NO_MEASURED_CHILD' else 1


if __name__ == '__main__':
    raise SystemExit(main())
