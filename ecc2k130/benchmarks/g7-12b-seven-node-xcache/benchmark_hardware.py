"""Capture benchmark hardware locally and reject comparisons across host types.

No cloud API or instance credentials are queried. A rate belongs to one GPU
on one provider/instance type; visible GPU count is not an aggregate rate.
"""
import csv
from datetime import datetime, timezone
import hashlib
import json
import os
import re
from pathlib import Path
import subprocess


def platform_identity(environ=None, dmi_root=Path('/sys/class/dmi/id')):
    env = os.environ if environ is None else environ
    provider = env.get('ECC_BENCH_PROVIDER')
    instance = env.get('ECC_BENCH_INSTANCE_TYPE')
    def read(name):
        try:
            return (Path(dmi_root)/name).read_text().strip()
        except OSError:
            return ''
    vendor, product = read('sys_vendor'), read('product_name')
    if provider is None and vendor == 'Amazon EC2':
        provider = 'aws'
    if provider == 'aws' and instance is None and vendor == 'Amazon EC2':
        instance = product if re.fullmatch(r'[a-z0-9-]+\.[a-z0-9]+', product) else None
    return dict(provider=provider or 'unknown', instance_type=instance,
                hardware_type=env.get('ECC_BENCH_HARDWARE_TYPE'),
                identity_source='environment' if env.get('ECC_BENCH_PROVIDER') or env.get('ECC_BENCH_INSTANCE_TYPE') else 'local DMI')


def capture_hardware():
    identity = platform_identity()
    # The client uses CUDA device zero. With a visibility mask that is the
    # first named physical GPU; nvidia-smi itself does not honor that mask.
    selector = os.environ.get('CUDA_VISIBLE_DEVICES', '0').split(',')[0].strip()
    if not selector or selector == '-1':
        raise ValueError('GPU benchmark has no visible CUDA device')
    fields = ['name', 'uuid', 'driver_version', 'power.limit', 'memory.total', 'compute_cap']
    p = subprocess.run(['nvidia-smi', '-i', selector, '--query-gpu='+','.join(fields),
                        '--format=csv,noheader,nounits'], capture_output=True, text=True, check=True)
    rows = list(csv.reader(p.stdout.strip().splitlines()))
    if len(rows) != 1 or len(rows[0]) != len(fields):
        raise ValueError('Expected exactly one identified benchmark GPU')
    name, uuid, driver, power, memory, capability = [v.strip() for v in rows[0]]
    def number(text, kind=float):
        try:
            return kind(text)
        except ValueError:
            return None
    identity.update(captured_at_utc=datetime.now(timezone.utc).isoformat(), gpu_name=name, gpu_uuid=uuid, driver_version=driver,
                    power_limit_watts=number(power), memory_mib=number(memory, int),
                    compute_capability=capability, gpu_selector=selector,
                    measured_gpu_count=1, measurement_scope='single_gpu',
                    cpu_affinity_count=len(os.sched_getaffinity(0)) if hasattr(os, 'sched_getaffinity') else os.cpu_count())
    return identity


def hardware_key(hardware):
    """Unknown instance types never match known types or each other implicitly."""
    provider = hardware.get('provider')
    host_type = hardware.get('instance_type') or hardware.get('hardware_type')
    if not provider or provider == 'unknown' or not host_type:
        raise ValueError('Benchmark needs an explicit provider and instance or hardware type')
    fields = ('provider', 'instance_type', 'hardware_type', 'gpu_name', 'driver_version',
              'power_limit_watts', 'memory_mib', 'compute_capability', 'measured_gpu_count',
              'measurement_scope', 'cpu_affinity_count')
    if not hardware.get('gpu_name') or hardware.get('measured_gpu_count') != 1:
        raise ValueError('This benchmark comparison requires one identified GPU')
    return hashlib.sha256(json.dumps({k:hardware.get(k) for k in fields}, sort_keys=True).encode()).hexdigest()


def check_destination(path, hardware):
    """A reused label must not overwrite a different hardware type's value."""
    path = Path(path)
    if not path.exists():
        return
    previous = json.loads(path.read_text())
    if 'hardware' not in previous:
        raise ValueError('Existing benchmark has no hardware identity; use a new label or directory')
    if hardware_key(previous['hardware']) != hardware_key(hardware):
        raise ValueError('Benchmark label belongs to different hardware; use a new label or directory')


def require_matched_hardware(rows):
    if not rows:
        raise ValueError('No benchmark rows')
    keys = {hardware_key(row['hardware']) for row in rows}
    if len(keys) != 1:
        raise ValueError('Cannot combine benchmarks from different hardware or instance types')
    return next(iter(keys))


if __name__ == '__main__':
    print(json.dumps(capture_hardware(), indent=2))
