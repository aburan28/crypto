"""Bounded, serial, fresh-process measurement; never overwrite a run."""
import argparse
import ctypes
import hashlib
import json
import os
from pathlib import Path
import platform
import random
import resource
import statistics
import subprocess
import sys
import time
import tracemalloc

from backends import compute, semantic_result, structure, verify

ROOT = Path(__file__).resolve().parent


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def read_optional(path):
    try:
        return Path(path).read_text().strip()
    except OSError as error:
        return str(error)


def affinity():
    allowed = sorted(os.sched_getaffinity(0))
    os.sched_setaffinity(0, {allowed[0]})
    status = {'cpu': allowed[0], 'affinity': sorted(os.sched_getaffinity(0)),
              'exclusive_cpu': False}
    try:
        lib = ctypes.CDLL('libnuma.so.1', use_errno=True)
        mask = ctypes.c_ulong(1)
        rc = lib.set_mempolicy(2, ctypes.byref(mask), ctypes.c_ulong(1))
        status['membind_node0'] = 'success' if rc == 0 else os.strerror(ctypes.get_errno())
    except (OSError, AttributeError) as error:
        status['membind_node0'] = str(error)
    return status


def worker(case_index, arm, traced):
    pin = affinity()
    case = json.loads((ROOT / 'inputs.json').read_text())[case_index]
    backend, strategy = arm.split('/')
    runs, fingerprints = [], []
    peak = None
    for _ in range(1 if traced else 3):
        if traced:
            tracemalloc.start()
        start = time.perf_counter()
        result = compute(case['n'], case['generators'], backend=backend, strategy=strategy)
        computed = time.perf_counter()
        if traced:
            _, peak = tracemalloc.get_traced_memory()
            tracemalloc.stop()
        before_verify = time.perf_counter()
        verification = verify(case['n'], case['generators'], result)
        verified = time.perf_counter()
        fingerprints.append(digest(semantic_result(result)))
        runs.append({'compute_s': computed - start, 'verify_s': verified - before_verify,
                     'compute_plus_verify_s': computed - start + verified - before_verify})
    if len(set(fingerprints)) != 1:
        raise AssertionError('nondeterministic output/counters')
    return {'arm': arm, 'case': case['name'], 'runs': runs, 'affinity': pin,
            'semantic_sha256': fingerprints[0], 'stats': result['stats'],
            'verified': verification, 'peak_compute_python_bytes': peak,
            'process_peak_rss_kib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}


def system_snapshot():
    info = subprocess.run(['lscpu', '-J'], text=True, capture_output=True)
    processes = subprocess.run(['ps', '-eo', 'pid,pcpu,comm', '--sort=-pcpu'],
                               text=True, capture_output=True)
    return {'python': sys.version, 'platform': platform.platform(),
            'cpu': json.loads(info.stdout) if info.returncode == 0 else info.stderr,
            'allowed_affinity': sorted(os.sched_getaffinity(0)),
            'loadavg': os.getloadavg(), 'top_processes': processes.stdout.splitlines()[:15],
            'meminfo': read_optional('/proc/meminfo'),
            'cpu_quota': read_optional('/sys/fs/cgroup/cpu.max'),
            'memory_limit': read_optional('/sys/fs/cgroup/memory.max'),
            'numa_online': read_optional('/sys/devices/system/node/online'),
            'ddr_generation': None,
            'hardware_scope': 'virtualized Linux x86-64; no exclusivity or physical topology guarantee'}


def median_interval(values):
    rng = random.Random(912)
    medians = sorted(statistics.median(rng.choices(values, k=len(values))) for _ in range(10000))
    return [medians[249], medians[9749]]


def main(output):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    cases = json.loads((ROOT / 'inputs.json').read_text())
    report = {'before': system_snapshot(), 'cases': [], 'failures': [],
              'source_sha256': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                                for p in sorted(ROOT.rglob('*'))
                                if p.is_file() and p.suffix in ('.py', '.json', '.md')
                                and 'runs' not in p.parts},
              'end_to_end_cryptographic_cost': None, 'end_to_end_cryptographic_speedup': None}
    (output / 'host_before.json').write_text(json.dumps(report['before'], indent=2)+'\n')
    serial = 0
    def sample(index, arm, traced=False):
        nonlocal serial
        command = [sys.executable, str(ROOT / 'measure.py'), '--worker', str(index), '--arm', arm]
        if traced:
            command.append('--trace')
        serial += 1
        file = output / f'sample-{serial:04d}.json'
        try:
            run = subprocess.run(command, capture_output=True, text=True, timeout=30,
                                 env={**os.environ, 'PYTHONHASHSEED': '0', 'OMP_NUM_THREADS': '1'})
            if run.returncode:
                raise RuntimeError(run.stderr)
            value = json.loads(run.stdout)
        except Exception as error:
            failure = {'command': command, 'error': str(error)}
            file.write_text(json.dumps(failure, indent=2)+'\n')
            report['failures'].append(failure)
            (output / 'results.json').write_text(json.dumps(report, indent=2)+'\n')
            raise
        file.write_text(json.dumps(value, indent=2)+'\n')
        value['receipt'] = file.name
        return value
    def mean(sample, metric):
        return statistics.mean(r[metric] for r in sample['runs'])
    for index, case in enumerate(cases):
        entry = {**case, 'structure': structure(case['n'], case['generators']),
                 'aa': [], 'schedules': {}}
        expected = {}
        def confirm(s):
            strategy = s['arm'].split('/')[1]
            old = expected.setdefault(strategy, s['semantic_sha256'])
            if old != s['semantic_sha256']:
                raise AssertionError('backend output or algebraic-work mismatch')
        for _ in range(5):
            a, b = sample(index, 'sparse/frontier'), sample(index, 'sparse/frontier')
            confirm(a); confirm(b)
            entry['aa'].append({'a': a['receipt'], 'b': b['receipt'],
                                'b_over_a': mean(b, 'compute_s') / mean(a, 'compute_s')})
        for strategy in ('exhaustive', 'frontier'):
            pairs, samples = [], {'sparse': [], 'packed': []}
            for round_index in range(7):
                pair = {}
                for backend in (('sparse', 'packed') if round_index % 2 == 0 else ('packed', 'sparse')):
                    s = sample(index, f'{backend}/{strategy}')
                    confirm(s)
                    pair[backend] = s
                    samples[backend].append(s)
                pairs.append({metric: mean(pair['packed'], metric) / mean(pair['sparse'], metric)
                              for metric in ('compute_s', 'compute_plus_verify_s')})
            summary = {'pairs': pairs, 'arms': {}}
            for metric in ('compute_s', 'compute_plus_verify_s'):
                ratios = [p[metric] for p in pairs]
                summary[metric + '_packed_over_sparse'] = {
                    'median': statistics.median(ratios), 'bootstrap_95pct': median_interval(ratios)}
            for backend in ('sparse', 'packed'):
                arm = samples[backend]
                memory = sample(index, f'{backend}/{strategy}', traced=True)
                confirm(memory)
                summary['arms'][backend] = {
                    'receipts': [s['receipt'] for s in arm], 'traced_receipt': memory['receipt'],
                    'stats': arm[0]['stats'], 'semantic_sha256': arm[0]['semantic_sha256'],
                    'compute_s': {'median': statistics.median(mean(s, 'compute_s') for s in arm),
                                  'minimum': min(mean(s, 'compute_s') for s in arm)},
                    'compute_plus_verify_s': {'median': statistics.median(mean(s, 'compute_plus_verify_s') for s in arm),
                                             'minimum': min(mean(s, 'compute_plus_verify_s') for s in arm)},
                    'peak_compute_python_bytes': memory['peak_compute_python_bytes'],
                    'process_peak_rss_kib': [s['process_peak_rss_kib'] for s in arm]}
            entry['schedules'][strategy] = summary
        report['cases'].append(entry)
        (output / 'results.json').write_text(json.dumps(report, indent=2)+'\n')
        print(case['name'], {s: round(v['compute_s_packed_over_sparse']['median'], 3)
                             for s, v in entry['schedules'].items()}, flush=True)
    report['after'] = system_snapshot()
    (output / 'results.json').write_text(json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--worker', type=int)
    parser.add_argument('--arm')
    parser.add_argument('--trace', action='store_true')
    parser.add_argument('--output', default='runs/replay')
    args = parser.parse_args()
    if args.worker is not None:
        print(json.dumps(worker(args.worker, args.arm, args.trace)))
    else:
        main(args.output)
