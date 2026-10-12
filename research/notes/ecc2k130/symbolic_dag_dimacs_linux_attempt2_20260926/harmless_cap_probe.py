#!/usr/bin/env python3
"""Linux-only RLIMIT controls; never imports or launches the measured producer."""
from __future__ import annotations

import argparse
import json
import os
import platform
import resource
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
PROBE = '''import json, mmap, resource, sys
cap = int(sys.argv[1])
soft, hard = resource.getrlimit(resource.RLIMIT_AS)
assert (soft, hard) == (cap, cap), (soft, hard)
try:
    mapping = mmap.mmap(-1, cap + 4096)
except (OSError, MemoryError) as exc:
    print(json.dumps({'limit': [soft, hard], 'over_cap_allocation': 'REJECTED',
                      'exception': type(exc).__name__}, sort_keys=True))
else:
    mapping.close()
    raise SystemExit('NOT_ADMITTED: over-cap mapping succeeded')
'''


def limits(cap: int):
    def apply() -> None:
        resource.setrlimit(resource.RLIMIT_AS, (cap, cap))
    return apply


def probe(cap: int) -> dict:
    harmless = subprocess.run(['/usr/bin/true'], preexec_fn=limits(cap),
                              capture_output=True, timeout=15, check=True)
    assert harmless.stdout == b'' and harmless.stderr == b''
    child = subprocess.run([sys.executable, '-c', PROBE, str(cap)],
                           preexec_fn=limits(cap), capture_output=True,
                           text=True, timeout=15, check=True)
    observed = json.loads(child.stdout)
    assert observed['limit'] == [cap, cap]
    assert observed['over_cap_allocation'] == 'REJECTED'
    return {'cap_bytes': cap, 'true_exit_code': harmless.returncode,
            'python_exit_code': child.returncode, 'observation': observed}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    out = args.out.resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.exists():
        raise SystemExit('NOT_ADMITTED: probe receipt path already exists')
    prepare = json.loads((HERE / 'PREPARE.json').read_text())
    receipt = {'schema': 'k0-dag-dimacs-linux-rlimit-probe-v1',
               'utc': datetime.now(timezone.utc).isoformat(),
               'python': sys.version, 'python_executable': sys.executable,
               'platform': platform.platform(), 'machine': platform.machine(),
               'runner_image_os': os.getenv('ImageOS'),
               'runner_image_version': os.getenv('ImageVersion'),
               'checkout_head': subprocess.check_output(
                   ['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
               'caps': [], 'decision': 'PREPARE_FAILURE'}
    try:
        if sys.platform != 'linux' or platform.machine() != 'x86_64':
            raise RuntimeError('NOT_ADMITTED: target must be Linux x86_64')
        release = dict(line.split('=', 1) for line in
                       Path('/etc/os-release').read_text().splitlines()
                       if '=' in line)
        if (release.get('ID', '').strip('"') != 'ubuntu' or
                release.get('VERSION_ID', '').strip('"') != '24.04' or
                not receipt['runner_image_os'] or
                not receipt['runner_image_version']):
            raise RuntimeError('NOT_ADMITTED: exact Ubuntu 24.04 Actions image unavailable')
        for cap in prepare['target']['hard_rlimit_as_caps_bytes']:
            receipt['caps'].append(probe(cap))
        receipt['decision'] = 'PASS_HARMLESS_CAPS_ONLY'
    except Exception as exc:
        receipt['error'] = f'{type(exc).__name__}: {exc}'
    out.write_text(json.dumps(receipt, sort_keys=True, indent=2) + '\n')
    print(json.dumps({'decision': receipt['decision'], 'caps': len(receipt['caps'])},
                     sort_keys=True))
    return 0 if receipt['decision'] == 'PASS_HARMLESS_CAPS_ONLY' else 1


if __name__ == '__main__':
    raise SystemExit(main())
